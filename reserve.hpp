#pragma once

// Stores information relating to the reserve configuration.
#include <limits>
#include <random>

#include "marzone.hpp"
#include "pu.hpp"
#include "zones.hpp"

namespace marzone
{
  using namespace std;

  class Reserve
  {
  public:
    Reserve(Species &spec, Zones &zones)
    {
      InitializeZoneSpec(spec.spno, zones.zoneCount);
      speciesAmounts.resize(spec.spno, {}); // init amounts struct
    }

    // Sets all amounts and occurences to 0
    void InitializeZoneSpec(uint64_t spno, uint64_t zoneCount)
    {
      zoneSpec.resize(spno * zoneCount);

      for (int j = 0; j < spno; j++)
      {
        for (int i = 0; i < zoneCount; i++)
        {
          zoneSpec[(j * zoneCount) + i].amount = 0;
          zoneSpec[(j * zoneCount) + i].occurrence = 0;
        }
      }
    }

    void InitializeSolution(uint64_t puno)
    {
      solution.assign(puno, 0);
    }

    void RandomiseSolution(Pu &pu, mt19937 &rngEngine)
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
          solution[i] = pu.puZone[i][zoneInd] - 1; //TODO - change this depending on later design, but for now it is the zero-indexed zone.
        }
      }
    }

    void ComputeSpeciesAmounts(Pu &pu, Species &spec, Zones &zones)
    {
      int isp, ism, iZoneSpecIndex;
      double rContribAmount;

      // Set species amounts to 0 and set target2 species.
      for (int isp = 0; isp < speciesAmounts.size(); isp++)
      {
        speciesAmounts[isp].amount = 0;
        speciesAmounts[isp].occurrence = 0;
        if (spec.specList[isp].target2)
          SpeciesAmounts4(isp, spec, clumptype); //TODO
      }

      for (int ipu = 0; ipu < solution.size(); ipu++)
      {
        if (pu.puList[ipu].richness)
        {
          if (solution[ipu] != 0)
          {
            for (int i = 0; i < pu.puList[ipu].richness; i++)
            {
              ism = pu.puList[ipu].offset + i;
              isp = pu.puvspr[ism].spindex;
              if (spec.specList[isp].target2 == 0)
              {
                rContribAmount = 0.0; // Todo
                if (rContribAmount)
                {
                  speciesAmounts[isp].amount += rContribAmount;
                  speciesAmounts[isp].occurrence += (rContribAmount > 0);

                  iZoneSpecIndex = (isp * zones.zoneCount) + (solution[ipu] - 1);
                  zoneSpec[iZoneSpecIndex].amount += pu.puvspr[ism].amount;
                  zoneSpec[iZoneSpecIndex].occurrence++;
                }
              }
            }
          }
        }
      }
    }

    // * * * * **** Change in penalty for moving single PU between zones ******
    // we need to know the penalty with ipu in its existing zone and the penalty with ipu in its proposed zone.
    // change in penalty = penalty in proposed configuration - penalty in existing configuration
    double ComputeChangePenalty(Pu& pu, Zones& zones, Species& spec, double& rShortfall, int ipu, int iPreZone, int iPostZone) {
      int ism, isp, iArrayIndex, iNewOccurrence, iCurrentShortfall, iNewShortfall;
      double rShortFraction, rNewShortFraction, rOldShortfall, rNewShortfall, rNewAmount, rSumDeltaPenalty = 0, rDeltaPenalty,
            rCurrentContribAmount, rNewContribAmount;

      rShortfall = 0;
      vector<int> zonesToTarget{iPreZone, iPostZone};

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
            rNewContribAmount = pu.puvspr[ism].amount*zones.GetZoneContrib(ipu, pu.puno, isp, iPreZone);

            rNewAmount = speciesAmounts[isp].amount + rNewContribAmount - rCurrentContribAmount;
            iNewOccurrence = speciesAmounts[isp].occurrence + (rNewContribAmount > 0) - (rCurrentContribAmount > 0);

            if (spec.specList[isp].target > 0)
            {
              if (spec.specList[isp].target > speciesAmounts[isp].amount)
              {
                rOldShortfall += spec.specList[isp].target - speciesAmounts[isp].amount;
                rShortFraction += (spec.specList[isp].target - speciesAmounts[isp].amount) / spec.specList[isp].target;
                iCurrentShortfall++;
              }

              if (spec.specList[isp].target > rNewAmount)
              {
                rNewShortfall += spec.specList[isp].target - rNewAmount;
                rNewShortFraction += (spec.specList[isp].target - rNewAmount) / spec.specList[isp].target;
                iNewShortfall++;
              }

              /* TODO - enable
              #ifdef DEBUG_PEW_CHANGE_PEN
                            sprintf(debugline, "%i,%i,%i,%i,%i,%i,%i,%g,%g,%g,%g,%g,%g,%g,%g,%i,%i,0\n", iIteration, ipu, isp, pu[ipu].id, spec[isp].name, R[ipu], iZone, spec[isp].target, SM[ism].amount, spec[isp].amount, rNewAmount, rShortfall, rNewShortfall, rShortFraction, rNewShortFraction, iCurrentShortfall, iNewShortfall);
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

            // compute existing & proposed shortfall for this feature across any relevant zone targets
            for (int k: zonesToTarget)
            {
              iArrayIndex = (isp * zones.zoneCount) + k;
              double currZoneAmount = zoneSpec[iArrayIndex].amount;
              double currZoneOcc = zoneSpec[iArrayIndex].occurrence;

              // compute amount of feature if change is made
              if (k == iPreZone) // zone is existing zone, reduce zone amount by amount at site
              {
                rNewAmount = currZoneAmount - pu.puvspr[ism].amount;
                iNewOccurrence = currZoneOcc - 1;
              }
              else // zone is proposed zone, increase zone amount by amount at site
              {
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
                                    sprintf(debugline, "%i,%i,%i,%i,%i,%i,%i,%g,%g,%g,%g,%g,%g,%g,%g,%i,%i,1\n", iIteration, ipu, isp, pu[ipu].id, spec[isp].name, R[ipu], iZone, _ZoneTarget[iArrayIndex].target, SM[ism].amount, ZoneSpec[iArrayIndex].amount, rNewAmount, rShortfall, rNewShortfall, rShortFraction, rNewShortFraction, iCurrentShortfall, iNewShortfall);
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

            rDeltaPenalty = spec.specList[isp].penalty * spec.specList[isp].spf * (rNewShortFraction - rShortFraction);
/* TODO - renable
#ifdef DEBUG_PEW_CHANGE_PEN
            // rDeltaPenalty,spec[isp].penalty,spec[isp].spf,rNewShortFraction,rShortFraction
            AppendDebugTraceFile("rDeltaPenalty,spec.penalty,spec.spf,rNewShortFraction,rShortFraction\n");
            sprintf(debugline, "%g,%g,%g,%g,%g\n", rDeltaPenalty, spec[isp].penalty, spec[isp].spf, rNewShortFraction, rShortFraction);
            AppendDebugTraceFile(debugline);
#endif
*/

            rSumDeltaPenalty += rDeltaPenalty;
            rShortfall += rNewShortfall - rOldShortfall;

          } // Only worry about PUs where species occurs

/* TODO - renable
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

      /* TODO - renable
        #ifdef DEBUG_CHANGE_PEN
       sprintf(debugline,"ChangePen end rSumDeltaPenalty %g\n",rSumDeltaPenalty);
       AppendDebugTraceFile(debugline);
       #endif
      */

      return rSumDeltaPenalty;
    }

    // Computes the change in value of a pu from preZone to postZone
    scost CheckChangeValue(int puindex, int iPreZone, int iPostZone, Pu& pu, Zones& zones, Species& spec, double costthresh) {
      scost change = {};
      int imode = 1;
      double threshpen = 0;

      // Change cost
      auto& puCostArray = pu.puList[puindex].costBreakdown;
      change.cost = zones.AggregateTotalCostByPuAndZone(iPostZone, puCostArray) - zones.AggregateTotalCostByPuAndZone(iPreZone, puCostArray);

      // Connection cost
      change.connection = pu.ConnectionCost2(zones, puindex, imode, solution, iPostZone) 
          - pu.ConnectionCost2(zones, puindex, imode, solution, iPreZone);

      /* Todo - see if we need to remove
      if (threshtype == 1)
      {
        tchangeconnection = change->connection;
        tresconnection = reserve->connection;
        change->connection = 0;
        reserve->connection = 0;
      }
      */

      double rShortfall;
      change.penalty = ComputeChangePenalty(pu, zones, spec, rShortfall, puindex, iPreZone, iPostZone);

      if (costthresh) { // Threshold Penalty for costs

      }
    }

    void WriteSolution(string filename, Pu &pu, int imode)
    {
      ofstream myfile;
      myfile.open(filename);

      if (imode > 1)
        myfile << "planning_unit,zone\n";

      for (int i = 0; i < pu.puno; i++)
      {
        if (imode > 1)
          myfile << pu.puList[i].id << "," << solution[i] + 1 << "\n";
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

    // * * * * ****** Value of a Zonation System * * * *
    void EvaluateObjectiveValue(Pu &pu, Species &spec, Zones &zones)
    {
      // Evaluate the existing solution and update the scost values
      objective.cost = 0, objective.penalty = 0, objective.connection = 0, objective.shortfall = 0;
      double rShortfall, iShortfall, rShortFraction, rTotalShortfall;
      int specZoneIndex, iMissingFeatures;

      //if (aggexist)
      //    SetSpeciesClumps(puno,R,spec,pu,SM,connections,clumptype);

      ComputeSpeciesAmounts(pu, spec, zones);

      for (int i = 0; i < spec.spno; i++) {
        rShortfall = 0;
        iShortfall = 0;
        rShortFraction = 0;

        // shortfall with respect to overall targets
        if (spec.specList[i].target > 0)
          if (spec.specList[i].target > speciesAmounts[i].amount)
          {
            rShortfall += spec.specList[i].target - speciesAmounts[i].amount;
            rShortFraction += (spec.specList[i].target - speciesAmounts[i].amount) / spec.specList[i].target;
            iShortfall++;

            objective.shortfall += spec.specList[i].target - speciesAmounts[i].amount;
          }

        if (spec.specList[i].targetocc > 0)
          if (spec.specList[i].targetocc > speciesAmounts[i].occurrence)
          {
            rShortfall += spec.specList[i].targetocc - speciesAmounts[i].occurrence;
            rShortFraction += (spec.specList[i].targetocc - speciesAmounts[i].occurrence) / spec.specList[i].targetocc;
            iShortfall++;

            objective.shortfall += spec.specList[i].targetocc - speciesAmounts[i].occurrence;
          }

        // shortfall with respect to zone targets (removing pu from existing zone)
        // loop through all zones to compute target achievement
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

        /* TODO - enable and move out
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
        objective.connection += pu.ConnectionCost2Linear(zones, j, 1, solution); // TODO - confirm imode 1
      }

      objective.total = objective.cost + objective.connection + objective.penalty;
    }

    // Writes amounts related to species amounts in this particular reserve to a file.
    void WriteSpeciesAmounts(string filename, Species &spec, Zones &zones, int imode, double misslevel)
    {
      ofstream myfile;
      myfile.open(filename);
      double rMPM, rTestMPM, rTarget, rAmountHeld;
      int rOccurrenceTarget, rOccurrencesHeld;
      string temp, d = imode > 1 ? "," : "\t";

      for (int isp = 0; isp < spec.spno; isp++)
      {
        rMPM = 1;

        // Write species statis
        myfile << spec.specList[isp].name << d << spec.specList[isp].sname << d << spec.specList[isp].target;

        if (imode > 1)
          myfile << d << spec.specList[isp].totalarea;
        else
          myfile << d << speciesAmounts[isp].amount;

        myfile << speciesAmounts[isp].amount << d << spec.specList[isp].targetocc << d << speciesAmounts[isp].occurrence;

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

    vector<zonespecstruct> zoneSpec;

    // Species data and amounts
    vector<reservespecies> speciesAmounts;

    // Pu configuration (assignment)
    vector<int> solution; // contains the zone index in the corresponding to a supplied Pu

    // Reserve metadata/objective values
    scost objective;

  private:
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
  };

} // namespace marzone