#pragma once

// Stores information relating to the reserve configuration.
#include <limits>
#include <cmath>
#include <random>
#include <utility>
#include <array>

#include "common.hpp"
#include "pu.hpp"
#include "zones.hpp"

namespace marzone
{
  using namespace std;

  // Tracks all computed changes on each move so they can be applied without recalculation.
  // Changes stored in pairs, usually first term = index of spec or zonespec, and second term is the value of change.
  typedef struct schange : scost
  {
    vector<pair<int,double>> specListChangeTarget; // list of spindex,target pairs.
    vector<int> specListChangeOcc; // index corresponds with above target vector.
    vector<pair<int, double>> zoneTargetChange;
    vector<int> zoneOccChange; // ind corresponds to zoneTargetChange.first
    vector<pair<int, double>> speciesClumpChange;
  }
  schange;

  class Reserve
  {
  public:
    Reserve() // Empty constructor.
    {
      InitializeObjective();
    }

    Reserve(Species &spec, int zoneCount, int clumptype, int id = 0) : clumptype(clumptype), id(id)
    {
      InitializeZoneSpec(spec.spno, zoneCount);
      speciesAmounts.resize(spec.spno, {}); // init amounts struct
      if (spec.aggexist)
        speciesClump.resize(spec.spno);
      
      InitializeObjective();
    }

    // Clone existing reserve.
    Reserve(const Reserve& r) {
      solution = r.solution;
      zoneSpec = r.zoneSpec;
      speciesAmounts = r.speciesAmounts;
      speciesClump = r.speciesClump; 
      objective = r.objective;
      id = r.id;
      clumptype = r.clumptype;
    }

    // Replace existing reserve with given reserve.
    void Assign(const Reserve& r) {
      solution = r.solution;
      zoneSpec = r.zoneSpec;
      speciesAmounts = r.speciesAmounts;
      speciesClump = r.speciesClump;
      objective = r.objective;
      id = r.id;
      clumptype = r.clumptype;
    }

    // Initializes a change structure to re-use for change calculations
    schange InitializeChange(Species& spec, Zones& zones) {
      schange change = {}; // init and pre-allocate
      change.specListChangeTarget.reserve(spec.spno);
      change.specListChangeOcc.reserve(spec.spno);

      if (zones.zoneTarget.size())
      {
        change.zoneTargetChange.reserve(spec.spno * zones.zoneCount);
        change.zoneOccChange.reserve(spec.spno * zones.zoneCount);
      }

      if (spec.aggexist)
        change.speciesClumpChange.reserve(spec.spno);

      return change;
    }

    void InitializeSolution(uint64_t puno)
    {
      solution.assign(puno, 0);
    }

    void RandomiseSolution(Pu &pu, mt19937 &rngEngine, int zoneCount)
    {
      uniform_int_distribution<int> random_dist(0, numeric_limits<int>::max());

      for (int i = 0; i < pu.puno; i++)
      {
        spustuff &puTerm = pu.puList[i];
        if (puTerm.fPULock == 1)
        {
          solution[i] = puTerm.iPULock;
          continue; // pu lock takes precedence above all.
        }

        solution[i] = 0; // default setting is 0

        if (puTerm.status > 0) // starting pu state detected
          solution[i] = puTerm.status;

        if (puTerm.numZones > 0)
        { // enforce puzone
          int zoneInd = random_dist(rngEngine) % puTerm.numZones;
          solution[i] = pu.puZone[i][zoneInd]; // solution array stores the zero-indexed zone.
        }
        else if (puTerm.numZones == 0)
        {
          // set to any random zone
          solution[i] = random_dist(rngEngine) % zoneCount;
        }
      }
    }

    // Recomputes and updates all speciesAmounts and zoneSpec, including clumps, from scratch.
    void ComputeSpeciesAmounts(Pu &pu, Species &spec, Zones &zones)
    {
      int isp, ism, iZoneSpecIndex;
      double rContribAmount;

      // Set species amounts to 0 and set target2 species.
      for (int isp = 0; isp < speciesAmounts.size(); isp++)
      {
        speciesAmounts[isp].amount = 0;
        speciesAmounts[isp].occurrence = 0;
        // Clear map clumps for this species
        if (spec.specList[isp].target2)
          ClearClumps(isp); // for target2 species, clear existing pu.
      }

      for (int ipu = 0; ipu < solution.size(); ipu++)
      {
        if (pu.puList[ipu].richness)
        {
          if (solution[ipu] >= 0)
          {
            for (int i = 0; i < pu.puList[ipu].richness; i++)
            {
              ism = pu.puList[ipu].offset + i;
              isp = pu.puvspr[ism].spindex;
              rContribAmount = pu.puvspr[ism].amount*zones.GetZoneContrib(ipu, pu.puno, isp, solution[ipu]);
              if (spec.specList[isp].target2 == 0)
              {
                if (pu.puvspr[ism].amount) 
                {
                  iZoneSpecIndex = (isp * zones.zoneCount) + solution[ipu];
                  zoneSpec[iZoneSpecIndex].amount += pu.puvspr[ism].amount;
                  zoneSpec[iZoneSpecIndex].occurrence += (pu.puvspr[ism].amount > 0);
                }

                if (rContribAmount)
                {
                  speciesAmounts[isp].amount += rContribAmount;
                  speciesAmounts[isp].occurrence += (rContribAmount > 0);
                }
              }
              else {
                // Add amount to species clumping.
                if (rContribAmount)
                  AddNewPuToClump(ipu, isp, rContribAmount);
              }
            }
          }
        }
      }

      //Compute species amounts for target2 species.
      for (int isp = 0; isp < speciesAmounts.size(); isp++)
      {
        if (spec.specList[isp].target2)
          SpeciesAmounts4(isp, spec.specList[isp].target2); // for target2 species, recompute
      }
    }

    // Counts missing species and proportion of lowest species
    int CountMissing(Species& spec, Zones& zones, double misslevel, double& rMinimumProportionMet) {
      double rProportionMet;
      int specCount = 0, iArrayIndex;

      rMinimumProportionMet = 1;
      for (int i = 0; i < spec.spno; i++)
      {
        rProportionMet = 1;

        if (spec.specList[i].target > 0)
        {
          if (speciesAmounts[i].amount < spec.specList[i].target)
          {
            rProportionMet = speciesAmounts[i].amount / spec.specList[i].target;

            if (rProportionMet < rMinimumProportionMet)
              rMinimumProportionMet = rProportionMet;
          }

          if (rProportionMet < misslevel)
          {
            specCount++;
            continue;
          }
        }

        if (spec.specList[i].targetocc > 0)
        {
          if (speciesAmounts[i].occurrence < spec.specList[i].targetocc)
          {
            rProportionMet = (double)speciesAmounts[i].occurrence / (double)spec.specList[i].targetocc;

            if (rProportionMet < rMinimumProportionMet)
              rMinimumProportionMet = rProportionMet;
          }

          if (rProportionMet < misslevel)
          {
            specCount++;
            continue;
          }
        }

        // only necessary if zoneTarget is set
        if (zones.zoneTarget.size())
          for (int j = 0; j < zones.zoneCount; j++)
          {
            iArrayIndex = (i * zones.zoneCount) + j;
            if (zones.zoneTarget[i][j].target > 0)
            {
              if (zoneSpec[iArrayIndex].amount < zones.zoneTarget[i][j].target)
              {
                rProportionMet = zoneSpec[iArrayIndex].amount / zones.zoneTarget[i][j].target;

                if (rProportionMet < rMinimumProportionMet)
                  rMinimumProportionMet = rProportionMet;

                if (rProportionMet < misslevel)
                {
                  specCount++;
                }
              }
            }
            if (zones.zoneTarget[i][j].occurrence > 0)
            {
              if (zoneSpec[iArrayIndex].occurrence < zones.zoneTarget[i][j].occurrence)
              {
                rProportionMet = zoneSpec[iArrayIndex].occurrence / zones.zoneTarget[i][j].occurrence;

                if (rProportionMet < rMinimumProportionMet)
                  rMinimumProportionMet = rProportionMet;

                if (rProportionMet < misslevel)
                {
                  specCount++;
                }
              }
            }
          }

        if (spec.specList[i].sepdistance && speciesAmounts[i].separation < 3)
        {
          specCount++; /* count species if not met separation and not already counted */
        }
      }

      return specCount;
    }

    // returns a formatted string based on zone distribution of the current solution
    string CountPuZones(Zones& zones) {
      vector<int> zoneDist(zones.zoneCount, 0);
      stringstream sLine;

      for (int i=0; i < solution.size(); i++) {
        zoneDist[solution[i]] += 1;
      }

      for (int i=0;i<zones.zoneCount;i++)
      {
        sLine << " " << zones.zoneNames[i].name << " " << zoneDist[i];
      }

      return sLine.str();
    }

    /* * * * * Set Penalties for a given Type 4 Species ***/
    /* Returns 1 if the species is a 'bad species' and -1 if it is a 'good species' */
    // Used for initial penalty calculation. and updates spec with the penalty.
    // NOTE: currently does no support sepnum/sepdistance feature.
    int ComputePenaltyType4(sspecies& spec, vector<penaltyTerm>& costTerms, int isp, double target, double target2, int targetocc) {
      int sepCount;

      /* Check first to see if I've already satisfied targets for this species with the current clump setup*/
      SpeciesAmounts4(isp, target2); // recompute amounts for this species.
      /* NOTE - reenable later if we have time to migrate the separation code.
      if (spec.specList[isp].sepnum > 0)
      {
        sepCount = CountSeparation2(isp);
        speciesAmounts[isp].separation = sepCount == -1 ? spec.specList[isp].sepnum : sepCount;
      }
      */

      if ((speciesAmounts[isp].amount >= target) 
        && (speciesAmounts[isp].occurrence >= targetocc)) {
          // Targets met for this spec can exit
          return -1;
      }

      // Go down costTerms until an amount is less than target2, or target is met
      bool targetMet = false;
      for (penaltyTerm &p : costTerms)
      {
        if (p.amount > target2)
        {
          speciesAmounts[isp].amount += p.amount;
          speciesAmounts[isp].occurrence++;
          spec.penalty += p.cost;
        }

        // Check if targets met
        if (speciesAmounts[isp].amount >= target && speciesAmounts[isp].occurrence >= targetocc)
        {
          targetMet = true;
          break;
        }
      }

      if (targetMet)
        return -1;
      else
      {
        double scaleFactor = 1;
        // scale up the penalty
        if (speciesAmounts[isp].amount == 0)
          speciesAmounts[isp].amount = numeric_limits<double>::epsilon(); 
        if (speciesAmounts[isp].amount < target)
          scaleFactor = target/speciesAmounts[isp].amount;
        if (speciesAmounts[isp].occurrence && speciesAmounts[isp].occurrence < targetocc)
          scaleFactor += (double) targetocc/ (double) speciesAmounts[isp].occurrence;
        spec.penalty *= scaleFactor;

        return 1;
      }
    }

    // * * * * **** Change in penalty for moving single PU between zones ******
    // we need to know the penalty with ipu in its existing zone and the penalty with ipu in its proposed zone.
    // change in penalty = penalty in proposed configuration - penalty in existing configuration
    double ComputeChangePenalty(Pu& pu, Zones& zones, Species& spec, schange& change, int ipu, int iPreZone, int iPostZone) {
      int ism, isp, iArrayIndex, iNewOccurrence, iCurrentShortfall, iNewShortfall;
      double rShortFraction, rNewShortFraction, rOldShortfall, rNewShortfall, rNewAmount, rSumDeltaPenalty = 0, rDeltaPenalty,
            rCurrentContribAmount, rNewContribAmount, t_diff, o_diff;

      change.shortfall = 0;
      array<int, 2> zonesToTarget{iPreZone, iPostZone};

      // clean up the vectors in the change object.
      change.specListChangeTarget.clear();
      change.specListChangeOcc.clear();
      change.speciesClumpChange.clear();
      change.zoneOccChange.clear();
      change.zoneTargetChange.clear();

      if (pu.puList[ipu].richness)
      {
        for (int i = 0; i < pu.puList[ipu].richness; i++)
        {
          ism = pu.puList[ipu].offset + i;
          isp = pu.puvspr[ism].spindex;
          if (pu.puvspr[ism].amount)
          {
            rOldShortfall = 0, rShortFraction = 0, iCurrentShortfall = 0;

            // init variables tracking shortfall in proposed zone
            rNewShortfall = 0, rNewShortFraction = 0, iNewShortfall = 0;

            // shortfall with respect to overall targets
            rCurrentContribAmount = pu.puvspr[ism].amount*zones.GetZoneContrib(ipu, pu.puno, isp, iPreZone);
            rNewContribAmount = pu.puvspr[ism].amount*zones.GetZoneContrib(ipu, pu.puno, isp, iPostZone);

            // Apply any target2 requirements on these contrib amounts.
            if (spec.specList[isp].target2) {
              change.speciesClumpChange.push_back(pair<int,double>(isp, rNewContribAmount));
              rCurrentContribAmount = rCurrentContribAmount > spec.specList[isp].target2 ? rCurrentContribAmount : PartialPen4(spec.specList[isp].target2, rCurrentContribAmount);
              rNewContribAmount = rNewContribAmount > spec.specList[isp].target2 ? rNewContribAmount : PartialPen4(spec.specList[isp].target2, rNewContribAmount);
            } 

            // Set change vectors
            t_diff = rNewContribAmount - rCurrentContribAmount;
            o_diff = (rNewContribAmount > 0) - (rCurrentContribAmount > 0);

            if (t_diff != 0) // only need to execute below if there is change.
            {
              change.specListChangeTarget.push_back(pair<int,double>(isp, t_diff));
              change.specListChangeOcc.push_back(o_diff);

              rNewAmount = speciesAmounts[isp].amount + t_diff;
              iNewOccurrence = speciesAmounts[isp].occurrence + o_diff;

              if (spec.specList[isp].target > 0)
              {
                if (spec.specList[isp].target > speciesAmounts[isp].amount)
                {
                  rOldShortfall = spec.specList[isp].target - speciesAmounts[isp].amount;
                  rShortFraction += rOldShortfall / spec.specList[isp].target;
                  iCurrentShortfall++;
                }

                if (spec.specList[isp].target > rNewAmount)
                {
                  rNewShortfall = spec.specList[isp].target - rNewAmount;
                  rNewShortFraction += rNewShortfall / spec.specList[isp].target;
                  iNewShortfall++;
                }

                /* 
              #ifdef DEBUG_PEW_CHANGE_PEN
                            sprintf(debugline, "%i,%i,%i,%i,%i,%i,%i,%g,%g,%g,%g,%g,%g,%g,%g,%i,%i,0\n", iIteration, ipu, isp, pu[ipu].id, spec[isp].name, R[ipu], iZone, spec[isp].target, SM[ism].amount, spec[isp].amount, rNewAmount, change.shortfall, rNewShortfall, rShortFraction, rNewShortFraction, iCurrentShortfall, iNewShortfall);
                            AppendDebugFile("debug_MarZone_PewChangePen.csv", debugline, fnames);
                            AppendDebugTraceFile("iteration,ipu,isp,puid,spid,Zone,newZone,Target,PUAmount,Amount,newAmount,Shortfall,newShortfall,rSF,rNSF,iCSF,iNSF,zone\n");
                            AppendDebugTraceFile(debugline);
              #endif
              */
              }

              if (spec.specList[isp].targetocc > 0)
              {
                if (spec.specList[isp].targetocc > speciesAmounts[isp].occurrence)
                {
                  rOldShortfall += spec.specList[isp].targetocc - speciesAmounts[isp].occurrence;
                  rShortFraction += (spec.specList[isp].targetocc - speciesAmounts[isp].occurrence) / spec.specList[isp].targetocc;
                  iCurrentShortfall++;
                }

                if (spec.specList[isp].targetocc > iNewOccurrence)
                {
                  rNewShortfall += spec.specList[isp].targetocc - iNewOccurrence;
                  rNewShortFraction += (spec.specList[isp].targetocc - iNewOccurrence) / spec.specList[isp].targetocc;
                  iNewShortfall++;
                }
              }
            }

            // compute existing & proposed shortfall for this feature across any relevant zone targets
            // only necessary if zoneTarget is set
            if (zones.zoneTarget.size())
              for (int k: zonesToTarget)
              {
                iArrayIndex = (isp * zones.zoneCount) + k;
                double currZoneAmount = zoneSpec[iArrayIndex].amount;
                double currZoneOcc = zoneSpec[iArrayIndex].occurrence;

                // compute amount of feature if change is made
                if (k == iPreZone) // zone is existing zone, reduce zone amount by amount at site
                {
                  change.zoneTargetChange.push_back(pair<int, double>(iArrayIndex, -pu.puvspr[ism].amount));
                  change.zoneOccChange.push_back(-1); // reduce occ by 1
                  rNewAmount = currZoneAmount - pu.puvspr[ism].amount;
                  iNewOccurrence = currZoneOcc - 1;
                }
                else // zone is proposed zone, increase zone amount by amount at site
                {
                  change.zoneTargetChange.push_back(pair<int, double>(iArrayIndex, pu.puvspr[ism].amount));
                  change.zoneOccChange.push_back(1); // add occ by one.
                  rNewAmount = currZoneAmount + pu.puvspr[ism].amount;
                  iNewOccurrence = currZoneOcc + 1;
                }

                // do we have areal zone target?
                double currZoneTarget = zones.zoneTarget[isp][k].target;
                if (currZoneTarget > 0)
                {
                  // compute existing shortfall
                  if (currZoneTarget > currZoneAmount)
                  {
                    rOldShortfall += currZoneTarget - currZoneAmount;
                    rShortFraction += (currZoneTarget - currZoneAmount) / currZoneTarget;
                    iCurrentShortfall++;
                  }

                  // compute proposed shortfall
                  if (currZoneTarget > rNewAmount)
                  {
                    rNewShortfall += currZoneTarget - rNewAmount;
                    rNewShortFraction += (currZoneTarget - rNewAmount) / currZoneTarget;
                    iNewShortfall++;
                  }

                  /*
                    #ifdef DEBUG_PEW_CHANGE_PEN
                                      sprintf(debugline, "%i,%i,%i,%i,%i,%i,%i,%g,%g,%g,%g,%g,%g,%g,%g,%i,%i,1\n", iIteration, ipu, isp, pu[ipu].id, spec[isp].name, R[ipu], iZone, zones.zoneTarget[i][j].target, SM[ism].amount, ZoneSpec[iArrayIndex].amount, rNewAmount, change.shortfall, rNewShortfall, rShortFraction, rNewShortFraction, iCurrentShortfall, iNewShortfall);
                                      AppendDebugFile("debug_MarZone_PewChangePen.csv", debugline, fnames);
                                      AppendDebugTraceFile("iteration,ipu,isp,puid,spid,Zone,newZone,ZoneTarget,PUAmount,Amount,newAmount,Shortfall,newShortfall,rSF,rNSF,iCSF,iNSF,zone\n");
                                      AppendDebugTraceFile(debugline);
                    #endif
                    */
                }

                // do we have occurrence zone target?
                double currZoneTargetOcc = zones.zoneTarget[isp][k].occurrence;
                if (currZoneTargetOcc > 0)
                {
                  // compute existing shortfall
                  if (currZoneTargetOcc > currZoneOcc)
                  {
                    rOldShortfall += currZoneTargetOcc - currZoneOcc;
                    rShortFraction += (currZoneTargetOcc - currZoneOcc) / currZoneTargetOcc;
                    iCurrentShortfall++;
                  }

                  // compute proposed shortfall
                  if (currZoneTargetOcc > iNewOccurrence)
                  {
                    rNewShortfall += currZoneTargetOcc - iNewOccurrence;
                    rNewShortFraction += (currZoneTargetOcc - iNewOccurrence) / currZoneTargetOcc;
                    iNewShortfall++;
                  }
                }
              }

            if (rNewShortFraction - rShortFraction != 0)
            {
              rDeltaPenalty = spec.specList[isp].penalty * spec.specList[isp].spf * (rNewShortFraction - rShortFraction);
  /* renable if needed
  #ifdef DEBUG_PEW_CHANGE_PEN
              // rDeltaPenalty,spec[isp].penalty,spec[isp].spf,rNewShortFraction,rShortFraction
              AppendDebugTraceFile("rDeltaPenalty,spec.penalty,spec.spf,rNewShortFraction,rShortFraction\n");
              sprintf(debugline, "%g,%g,%g,%g,%g\n", rDeltaPenalty, spec[isp].penalty, spec[isp].spf, rNewShortFraction, rShortFraction);
              AppendDebugTraceFile(debugline);
  #endif
  */

              rSumDeltaPenalty += rDeltaPenalty;
              change.shortfall += rNewShortfall - rOldShortfall;
            }
          } // Only worry about PUs where species occurs

/* renable if needed
#ifdef DEBUG_PEW_CHANGE_PEN
          sprintf(debugline, "rSumDeltaPenalty %g\n", rSumDeltaPenalty);
          AppendDebugTraceFile(debugline);
#endif

#ifdef DEBUG_CHANGE_PEN
          sprintf(debugline, "%i,%g,%g,%g,%g,%g,%g\n",
                  spec[isp].name, spec[isp].target, _ZoneTarget[(isp * iZoneCount) + R[ipu] - 1].target,
                  _ZoneTarget[(isp * iZoneCount) + iZone - 1].target, rOldShortfall, rNewShortfall, rDeltaPenalty);
          AppendDebugTraceFile(debugline);
#endif
*/
        }
      }

      /* renable if needed
        #ifdef DEBUG_CHANGE_PEN
       sprintf(debugline,"ChangePen end rSumDeltaPenalty %g\n",rSumDeltaPenalty);
       AppendDebugTraceFile(debugline);
       #endif
      */

      return rSumDeltaPenalty;
    }

    // Computes the change in value of a pu from preZone to postZone
    void CheckChangeValue(schange& change, int puindex, int iPreZone, int iPostZone, Pu& pu, Zones& zones, Species& spec, 
    double costthresh, double blm, double tpf1 = 0, double tpf2 = 0, double timeprop = 1) {
      int imode = 1;
      double threshpen = 0;

      // Change cost only if there's zone cost adjustment
      if (zones.availableZoneCost)
      {
        auto& puCostArray = pu.puList[puindex].costBreakdown;
        change.cost = zones.AggregateTotalCostByPuAndZone(iPostZone, puCostArray) - zones.AggregateTotalCostByPuAndZone(iPreZone, puCostArray);
      }

      // Connection cost
      if (pu.connectionsEntered)
        change.connection = zones.ConnectionCost2(pu, puindex, imode, solution, iPostZone, blm) 
            - zones.ConnectionCost2(pu, puindex, imode, solution, iPreZone, blm);

      change.penalty = ComputeChangePenalty(pu, zones, spec, change, puindex, iPreZone, iPostZone);

      if (costthresh)
      { // Threshold Penalty for costs
        if (objective.cost + objective.connection <= costthresh)
        {
          if (change.cost + change.connection + objective.cost + objective.connection <= costthresh)
            threshpen = 0;
          else
            threshpen = (change.cost + change.connection +
                         objective.cost + objective.connection - costthresh) *
                        ThresholdPenalty(tpf1, tpf2, timeprop);
        }
        else
        {
          if (change.cost + change.connection + objective.cost + objective.connection <= costthresh)
            threshpen = (objective.cost + objective.connection - costthresh) * ThresholdPenalty(tpf1, tpf2, timeprop);
          else
            threshpen = (change.cost + change.connection) * ThresholdPenalty(tpf1, tpf2, timeprop);
        }
      }

      change.threshpen = threshpen;
      change.total = change.cost + change.connection + change.penalty + change.threshpen;

      /*
      #ifdef DEBUGCHECKCHANGE
      if (iDebugMode)
      {
          sprintf(debugline,"%i,%i,%i,%g,%g,%g,%g,%g\n",ipu,pu[ipu].id,R[ipu],change->total,change->cost,change->connection,change->penalty,change->threshpen);
          AppendDebugFile("debug_MarZone_CheckChange.csv",debugline,fnames);
      }
      #endif

      #ifdef DEBUGCHECKCHANGE
      if (iDebugMode)
          AppendDebugTraceFile("CheckChange end\n");
      #endif
      */
      //return change;
    }

    // Applies changes described. Switches puindex to zone, and applies scost change.
    void ApplyChange(int ipu, int iZone, schange& change, Pu& pu, Zones& zones, Species& spec) {
      int i,ism,isp, iPreviousZone, iArrayIndexDestination;
      double rAmount, rCurrentContribAmount, rNewContribAmount;

      objective.cost += change.cost;
      objective.connection += change.connection;
      objective.penalty += change.penalty;
      objective.shortfall += change.shortfall;

      // iterate over amount/target changes and any clumping
      for (int i = 0; i < change.specListChangeTarget.size(); i++) {
        isp = change.specListChangeTarget[i].first; // contains spindex

        if (change.specListChangeTarget[i].second != 0) { //only need to worry if there's an actual change
            // Apply amounts as usual
            speciesAmounts[isp].amount += change.specListChangeTarget[i].second;
            speciesAmounts[isp].occurrence += change.specListChangeOcc[i];
        }
        /* feature does not work - disable for now.
        if (spec[isp].sepnum > 0) // Count separation but only if it is possible that it has changed
          if ((imode == 1 && spec[isp].separation < spec[isp].sepnum) || (imode == -1 && spec[isp].separation > 1))
            spec[isp].separation = CountSeparation2(isp, 0, NULL, puno, R, pu, SM, spec, 0);
        */
      }

      // Iterate through clump changes (if any)
      for (int i = 0; i < change.speciesClumpChange.size(); i++) {
        isp = change.specListChangeTarget[i].first; // contains spindex
        if (change.specListChangeTarget[i].second) {
          // add or modify existing pu clump
          AddNewPuToClump(ipu, isp, change.specListChangeTarget[i].second);
        }
        else {
          // Contrib is now 0, remove
          RemovePuFromClump(ipu, isp);
        }
      }

      // Iterate through zone changes
      int ind;
      for (int i = 0; i < change.zoneTargetChange.size(); i++) {
        ind = change.zoneTargetChange[i].first;
        zoneSpec[ind].amount += change.zoneTargetChange[i].second;
        zoneSpec[ind].occurrence += change.zoneOccChange[i];
      }

      solution[ipu] = iZone; // change the zone set in solution
      objective.total = objective.cost + objective.connection + objective.penalty;
    }

    // Writes the current solution to a file.
    void WriteSolution(string filename, Pu &pu, Zones& zones, int imode)
    {
      ofstream myfile;
      myfile.open(filename);

      if (imode > 1)
        myfile << "planning_unit,zone\n";

      for (int i = 0; i < pu.puno; i++)
      {
        if (imode > 1)
          myfile << pu.puList[i].id << "," << zones.IndexToId(solution[i]) << "\n";
      }

      myfile.close();
    }

    void WriteZoneConnectivitySum(string filename, Pu &pu, Zones &zones, int imode)
    {
      vector<vector<double>> ZCS = zones.InitializeZoneMatrix();
      //string delim = imode > 1 ? "," : "    ";

      ComputeZoneConnectivitySum(pu, ZCS); // populate ZCS

      ofstream myfile;
      myfile.open(filename);

      if (imode > 1)
        myfile << "\"Zone_Connectivity_Sum\"";
      for (int i = 0; i < zones.zoneCount; i++)
        myfile << ",\"" << zones.IndexToName(i) << "\"";
      myfile << "\n";

      // write a data row for each zone
      for (int i = 0; i < zones.zoneCount; i++)
      {
        myfile << "\"" << zones.IndexToName(i) << "\"";
        for (int j = 0; j < zones.zoneCount; j++)
          myfile << "," << ZCS[i][j];
        myfile << "\n";
      }

      myfile.close();
    }

    // Given zones and pu, and using current solution.
    void ComputeZoneConnectivitySum(Pu &pu, vector<vector<double>> &ZCS)
    {
      if (pu.connectionsEntered) // only applies if pu has connections
      {
        for (int i = 0; i < pu.puno; i++)
        {
          // fixed cost for ipu is between this zone and itself
          ZCS[solution[i]][solution[i]] += pu.connections[i].fixedcost;

          // traverse connections for this ipu
          for (sneighbour &p : pu.connections[i].first)
          {
            if (p.nbr > i) // avoid double counting connnections
            {
              if (solution[i] != solution[p.nbr]) // ignore internal connections within a zone
              {
                // connections are symmetric
                ZCS[solution[i]][solution[p.nbr]] += p.cost;
                ZCS[solution[p.nbr]][solution[i]] += p.cost;
              }
            }
          }
        }
      }
    }

    /* * * * * Greedy Species Penalty * * * * * * * * * * * */
    double GreedyPen(int ipu, Pu& pu, Species& spec)
    {
      double famount = 0.0, fold, newamount, newpen;

      // Get all spec contribs for the given ipu
      vector<double> specValuesForPu = pu.RtnAmountAllSpecAtPu(ipu, spec.spno);

      for (int i = 0; i < spec.spno; i++)
      {
        fold = (spec.specList[i].target - speciesAmounts[i].amount);
        if (fold > 0)
        {
          if (spec.specList[i].target2)
          {
            // Convert contribution down to tagret4 penalty version.
            double newSpecValue = PartialPen4(spec.specList[i].target2, specValuesForPu[i]);
            newpen = spec.specList[i].target - speciesAmounts[i].amount - newSpecValue;
          }
          else
          {
            newpen = spec.specList[i].target - speciesAmounts[i].amount - specValuesForPu[i];
          }
          newamount = (newpen < 0) ? 0 : newpen;
          famount += (newamount - fold) * spec.specList[i].spf;
        } /* Add new penalty if species isn't already in the system */
      }
      return (famount); /* Negative means decrease in amount missing */
    } /** Greedy Species Penalty **/

    /**** Partial Penalty for type 4 species ***/
    double PartialPen4(double target2, double amount)
    {
      if (amount >= target2)
        return (amount); /* I'm not a partial penalty */
      else
        switch (clumptype)
        {
        case 0:
          return (0.0); /* default step function */
        case 1:
          return (amount / 2.0); /* nicer step function */
        case 2:
          if (target2)
            return (amount / target2 * amount);
        default:
          return (0.0);
        }
    } /* Partial Penalty for type 4 species */

    // * * * * ****** Value of a Zonation System * * * *
    // Note! This call is expensive as it does a full recompute of the reserve value, based on the current solution.
    void EvaluateObjectiveValue(Pu &pu, Species &spec, Zones &zones, double blm)
    {
      // Evaluate the existing solution and update the scost values
      objective.cost = 0, objective.penalty = 0, objective.connection = 0, objective.shortfall = 0;
      double rShortfall, iShortfall, rShortFraction, rTotalShortfall = 0;
      int specZoneIndex, iMissingFeatures = 0;

      ComputeSpeciesAmounts(pu, spec, zones);

      for (int i = 0; i < spec.spno; i++) {
        rShortfall = 0;
        iShortfall = 0;
        rShortFraction = 0;

        // shortfall with respect to overall targets
        if (spec.specList[i].target > 0)
        {
          if (spec.specList[i].target > speciesAmounts[i].amount)
          {
            rShortfall += spec.specList[i].target - speciesAmounts[i].amount;
            rShortFraction += (spec.specList[i].target - speciesAmounts[i].amount) / spec.specList[i].target;
            iShortfall++;

            objective.shortfall += spec.specList[i].target - speciesAmounts[i].amount;
          }
        }

        if (spec.specList[i].targetocc > 0)
        {
          if (spec.specList[i].targetocc > speciesAmounts[i].occurrence)
          {
            rShortfall += spec.specList[i].targetocc - speciesAmounts[i].occurrence;
            rShortFraction += (spec.specList[i].targetocc - speciesAmounts[i].occurrence) / spec.specList[i].targetocc;
            iShortfall++;

            objective.shortfall += spec.specList[i].targetocc - speciesAmounts[i].occurrence;
          }
        }

        // shortfall with respect to zone targets (removing pu from existing zone)
        // loop through all zones to compute target achievement (only necessary if zoneTarget set)
        if (zones.zoneTarget.size())
          for (int k = 0; k < zones.zoneCount; k++)
          {
            specZoneIndex = (i * zones.zoneCount) + k;
            if (zones.zoneTarget[i][k].target > 0)
              if (zones.zoneTarget[i][k].target > zoneSpec[specZoneIndex].amount)
              {
                  rShortfall += zones.zoneTarget[i][k].target - zoneSpec[specZoneIndex].amount;
                  rShortFraction += (zones.zoneTarget[i][k].target - zoneSpec[specZoneIndex].amount) / zones.zoneTarget[i][k].target;
                  iShortfall++;

                  objective.shortfall += zones.zoneTarget[i][k].target - zoneSpec[specZoneIndex].amount;
              }

            if (zones.zoneTarget[i][k].occurrence > 0)
              if (zones.zoneTarget[i][k].occurrence > zoneSpec[specZoneIndex].occurrence)
              {
                  rShortfall += zones.zoneTarget[i][k].occurrence - zoneSpec[specZoneIndex].occurrence;
                  rShortFraction += (zones.zoneTarget[i][k].occurrence - zoneSpec[specZoneIndex].occurrence) / zones.zoneTarget[i][k].occurrence;
                  iShortfall++;

                  objective.shortfall += zones.zoneTarget[i][k].occurrence - zoneSpec[specZoneIndex].occurrence;
              }
          }

        rTotalShortfall += rShortfall;
        iMissingFeatures += iShortfall;

        objective.penalty += rShortFraction * spec.specList[i].penalty * spec.specList[i].spf;

        /*
          if (fDebugPenaltyNegative)
          {
              fprintf(DebugFile,"%i,%i,%g,%g,%g,%g\n",i,spec[i].name,rCurrentShortfall,spec[i].penalty,spec[i].spf,reserve->penalty);
              //"i,SPID,shortfall,spec_penalty,spf,reserve_penalty"
          }

          #ifdef DEBUG_ZONATION_COST
          sprintf(debugbuffer,"ZonationCost spid %i targ %g reserved %g Shortfall %g rShort %g missing features %i spec.pen %g spec.spf %g\n"
                              ,spec[i].name,spec[i].target,spec[i].amount,rShortfall,rCurrentShortfall,iShortfall,spec[i].penalty,spec[i].spf);
          AppendDebugTraceFile(debugbuffer);
          #endif
          */
      }

      for (int j = 0; j < pu.puno; j++)
      {
        objective.cost += zones.AggregateTotalCostByPuAndZone(solution[j], pu.puList[j].costBreakdown);
        if (pu.connectionsEntered)
        {
          objective.connection += zones.ConnectionCost2Linear(pu, j, 1, solution, blm); // confirm imode 1
        }
      }

      objective.total = objective.cost + objective.connection + objective.penalty;
    }

    // Writes amounts related to species amounts in this particular reserve to a file.
    // See zones.WriteZoneTargetHeaders for header structure.
    void WriteSpeciesAmounts(string filename, Species &spec, Zones &zones, int imode, double misslevel)
    {
      ofstream myfile;
      myfile.open(filename, ofstream::app); // append mode
      double rMPM, rTestMPM, rTarget, rAmountHeld;
      int rOccurrenceTarget, rOccurrencesHeld;
      string temp, d = imode > 1 ? "," : "\t";

      for (int isp = 0; isp < spec.spno; isp++)
      {
        rMPM = 1;

        // Write species status
        myfile << spec.specList[isp].name << d << spec.specList[isp].sname << d << spec.specList[isp].target;

        if (imode > 1)
          myfile << d << spec.specList[isp].totalarea;
        else
          myfile << d << speciesAmounts[isp].amount;

        myfile << d << speciesAmounts[isp].amount << d << spec.specList[isp].targetocc << d << speciesAmounts[isp].occurrence;

        temp = ReturnStringIfTargetMet(spec.specList[isp].target, speciesAmounts[isp].amount,
                                       spec.specList[isp].targetocc, speciesAmounts[isp].occurrence, rMPM, misslevel);

        if (imode > 1)
        {
          if (spec.specList[isp].sepnum)
          {
            if (speciesAmounts[isp].separation / spec.specList[isp].sepnum < misslevel)
              temp = "no";
          }
        }
        myfile << d << temp;

        int iZoneArrayIndex;
        if (zones.zoneTarget.size())
          for (int i = 0; i < zones.zoneCount; i++)
          {
            iZoneArrayIndex = (isp * zones.zoneCount) + i;
            rTarget = zones.zoneTarget[isp][i].target;
            rAmountHeld = zoneSpec[iZoneArrayIndex].amount;
            rOccurrenceTarget = zones.zoneTarget[isp][i].occurrence;
            rOccurrencesHeld = zoneSpec[iZoneArrayIndex].occurrence;

            myfile << d << rTarget << d << rAmountHeld << d << zones.zoneContribValues[iZoneArrayIndex] << d << rOccurrenceTarget << d << rOccurrencesHeld;

            temp = ReturnStringIfTargetMet(rTarget, rAmountHeld, rOccurrenceTarget, rOccurrencesHeld, rMPM, misslevel);
            myfile << d << temp;
          }

        myfile << d << rMPM << "\n";
      }

      myfile.close();
    }

    // Not threadsafe.
    void AppendSolutionsMatrix(string filename, int zoneCount, int delimType, int iIncludeHeaders) {
      stringstream text;
      string sDelimiter = delimType == 3 ? "," : "    ";

      for (int i = 1; i <= zoneCount; i++)
      {
        if (iIncludeHeaders == 1)
        {
          text << "Z" << i << "S" << id << sDelimiter; 
        }

        for (int j = 0; j < solution.size(); j++)
        {
          if (j > 0)
          text << sDelimiter;

          if (solution[j] == i-1)
            text << "1";
          else
            text << "0";
        }

        text << "\n";
      }

      // open file for as little time as possible.
      ofstream myfile;
      myfile.open(filename, ofstream::app);
      myfile << text.str();
      myfile.close();
    }

    // Not threadsafe.
    void AppendSolutionsMatrixZone(string filename, int iZone, int delimType, int iIncludeHeaders)
    {
      stringstream text;
      string sDelimiter = delimType == 3 ? "," : "    ";

      if (iIncludeHeaders == 1)
      {
        text << "S" << id << sDelimiter;
      }

      for (int j = 0; j < solution.size(); j++)
      {
        if (j > 0)
          text << sDelimiter;

        if (solution[j] == iZone)
          text << "1";
        else
          text << "0";
      }

      text << "\n";

      ofstream myfile;
      myfile.open(filename, ofstream::app);
      myfile << text.str();
      myfile.close();
    }

    // Returns a formatted string for zone count based on this reserve configuration
    void CountPuZones2(Zones& zones, int imode, string& sCounts) {
      stringstream s2; 
      string d = imode > 1 ? "," : "    ";

      vector<int> zCounts(zones.zoneCount, 0);
      for (int i = 0; i < solution.size(); i++) {
        zCounts[solution[i]]++;
      }

      for (int i = 0; i < zones.zoneCount; i++) {
        s2 << d << zCounts[i];
      }

      sCounts = s2.str();
    }

    void CostPuZones(Pu &pu, Zones &zones, string &sCount, int imode)
    {
      int i;
      string d = imode > 1 ? "," : "    ";
      double rZoneCost;
      stringstream sCounts;
      vector<double> zCosts(zones.zoneCount, 0.0);

      for (int i = 0; i < pu.puno; i++)
      {
        rZoneCost = zones.AggregateTotalCostByPuAndZone(solution[i], pu.GetCostBreakdown(i));
        zCosts[solution[i]] += rZoneCost;
      }

      for (int i = 0; i < zones.zoneCount; i++)
      {
        sCounts << d << zCosts[i];
      }

      sCount = sCounts.str();
    }

    vector<zonespecstruct> zoneSpec;

    // Species data and amounts
    vector<reservespecies> speciesAmounts;
    
    // spec -> map of pu contributing to it (for species4 only!). Amount contains the puvspr value with zonecontrib applied.
    vector<map<int,sclumps>> speciesClump; 

    // Pu configuration (assignment)
    vector<int> solution; // contains the zone index in the corresponding to a supplied Pu

    // Reserve metadata/objective values
    scost objective;
    int id;
    int clumptype;

  private:
    // This is a modified form of count separation where the user can specify any
    // maximum separation distance rather than just assuming a sep distance of three
    /* NOTE - separation does not seem to work. If there's time we can implement it later.
    int CountSeparation2(int isp, int sepdistance) {
      double targetdist = sepdistance*sepdistance;
      if (targetdist == 0)
        return -1; //Shortcut if sep not apply to this species
    }
    */

    // For a given species with target2 requirements, sum up its clumping pu amounts
    // and set into speciesAmounts.
    void SpeciesAmounts4(int isp, double target2) {
      double totalAmount = 0, ftemp;
      int totalOcc = 0;
      for (auto& [ipu, clump]: speciesClump[isp]) {
        ftemp = PartialPen4(target2, clump.amount);
        totalAmount += ftemp;
        totalOcc += clump.occs*(ftemp>0);
      }
    }

    /* Clumping functionality functions */
    // Sets a clump number for a species in this reserve
    void SetClumpSpecAtPu(int spindex, int puindex, int iSetClump)
    {
      sclumps& current = RtnClumpSpecAtPu(spindex, puindex);
      if (current.clumpid == 0)
        speciesClump[spindex][puindex] = {iSetClump, 0, 0};
      else 
        current.clumpid = iSetClump;
    }

    // overwrites an entire sclump.
    void SetClumpSpecAtPu(int spindex, int puindex, sclumps& current)
    {
        speciesClump[spindex][puindex] = current;
    }

    // define a "null" clump
    sclumps nullClump = {0,0,0};
    sclumps& RtnClumpSpecAtPu(int spindex, int puindex) {
      auto it = speciesClump[spindex].find(puindex);
      if (it != speciesClump[spindex].end())
        return it->second;
      return nullClump;
    }

    // Add new pu to clump of a species.
    void AddNewPuToClump(Pu& pu, int ipu, int isp) {
      double rAmount = pu.RtnAmountSpecAtPu(ipu, isp);
      AddNewPuToClump(ipu, isp, rAmount);
    }

    // Adds or updates existing clump with amount
    void AddNewPuToClump(int ipu, int isp, double amount) {
      sclumps& clump = RtnClumpSpecAtPu(isp, ipu); // keep reference so it can be changed later.
      if (clump.clumpid) {
        // For some reason, species already has this pu as its clump.
        if (amount == clump.amount)
          return; // nothing to do here, pu is already included
      }

      // Now set up the new clump. Set clump id to ipu.
      clump.clumpid = ipu;
      clump.amount = amount;
      clump.occs = (amount > 0);
    }

    // Removes a pu from clump in isp, does nothing if clump does not exist.
    void RemovePuFromClump(int ipu, int isp) {
      auto it = speciesClump[isp].find(ipu);
      if (it != speciesClump[isp].end())
        speciesClump[isp].erase(it);
    }

    void ClearClumps(int isp) {
      speciesClump[isp].clear();
    }

    // Sets all amounts and occurences to 0
    void InitializeZoneSpec(uint64_t spno, uint64_t zoneCount)
    {
      zoneSpec.assign(spno * zoneCount, {0,0});
    }

    string ReturnStringIfTargetMet(double targetArea, double area, int targetOcc, int occ, double &rMPM, double &misslevel)
    {
      string temp = "yes";
      double rTestMPM;

      if (targetArea)
      {
        if (area / targetArea < misslevel)
          temp = "no";

        rTestMPM = area / targetArea;
        if (rTestMPM < rMPM)
          rMPM = rTestMPM;
      }
      if (targetOcc)
      {
        if (occ / targetOcc < misslevel)
          temp = "no";

        rTestMPM = occ / targetOcc;
        if (rTestMPM < rMPM)
          rMPM = rTestMPM;
      }
      return temp;
    }

    double ThresholdPenalty(double tpf1,double tpf2,double timeprop) {
        if (tpf2 < 0)
          return(tpf1);
        else if (tpf1 == 0)
          return 0;

       return(tpf1*exp(tpf2*timeprop));
    }

    // initializes the objective struct and anything else related to it
    void InitializeObjective() {
      objective = {};
    }
  };

} // namespace marzone