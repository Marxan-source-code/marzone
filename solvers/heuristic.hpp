#pragma once

#include <limits>
#include <random>

#include "../marzone.hpp"
#include "../pu.hpp"
#include "../reserve.hpp"
#include "../zones.hpp"

namespace marzone {
using namespace std;

/* * * * ***** Main Heuristic Engine * * * * * * * * ****/
class Heuristic {
    public:
    Heuristic(mt19937& rngEngine, int heuristicMode) : rngEngine(rngEngine), heuristicMode(heuristicMode)
    {
    }

    void RunHeuristic(Reserve& r, Species& spec, Pu& pu, Zones& zones, double tpf1, double tpf2, double costthresh, double blm) {
        /**** Irreplacability ****/

        if (heuristicMode >= 6 && heuristicMode <= 7)
        {
            rare.resize(spec.spno, 0);
            SetAbundance(pu, spec);
        }

        if (heuristicMode >= 2 && heuristicMode <= 5) /* Rareness Setups */
        {
            rare.resize(spec.spno, 0);
            SetRareness(r, pu, spec);
        }

        schange change = r.InitializeChange(spec, zones);

        int bestpu = 0, iZone = 0, bestZone = 0, iPreviousR, iArrayIndex;
        double bestscore = 0, currscore;
        uniform_int_distribution<int> randomDist(0, numeric_limits<int>::max());
        uniform_real_distribution<double> rand1(0.0, 1.0);

        do
        {
            bestpu = 0;
            bestscore = 0;

            for (int i = 0; i < pu.puno; i++)
            {
                // choose a non-available zone to score this available site for
                // we are changing from the available zone to the non-available zone
                iPreviousR = r.solution[i];
                iZone = pu.RtnValidZoneForPu(i, iPreviousR, randomDist, rngEngine, zones.zoneCount);

                /* Set the score for the given Planning Unit */
                currscore = 1; /* null if no other mode set */
                if (heuristicMode == 0)
                    currscore = GreedyScore(i, r, pu, spec, zones, blm);

                if (heuristicMode == 1)
                {
                    r.CheckChangeValue(change, i, iPreviousR, iZone, pu, zones, spec, costthresh, blm, tpf1, tpf2);
                    currscore = change.total;
                }

                if (heuristicMode == 2 || heuristicMode == 3 || heuristicMode == 4 || heuristicMode == 5)
                {
                    currscore = RareScoreByMode(i, r, pu, spec, zones, blm);
                }

                if (heuristicMode == 6)
                    currscore = -ProdIrr(i, r, pu, spec, zones);
                if (heuristicMode == 7)
                    currscore = -SumIrr(i, r, pu, spec, zones);

                currscore *=(double) rand1(rngEngine)*0.001 + 1.0;
                if (!costthresh || pu.puList[i].cost + r.objective.cost <= costthresh)
                    if (currscore < bestscore)
                    {
                        bestpu = i;
                        bestscore = currscore;
                        bestZone = iZone;
                    } /** is this better (ie negative) than bestscore? **/
            }  /** I've looked through each pu to find best **/

            if (bestscore)
            {
                r.CheckChangeValue(change, bestpu, r.solution[bestpu], bestZone, pu, zones, spec, costthresh, blm, tpf1, tpf2);
                r.ApplyChange(bestpu, bestZone, change, pu, zones, spec);

                /* Different Heuristics might have different penalty effects */
                /* Old fashioned penalty and missing counting */
                r.objective.missing = 0;
                for (int i = 0; i < spec.spno; i++)
                {
                    for (int j = 0; j < zones.zoneCount; j++)
                    {
                        iArrayIndex = (i * zones.zoneCount) + j;
                        if (r.zoneSpec[iArrayIndex].amount < zones.zoneTarget[i][j].target)
                            r.objective.missing++;
                    }
                    if (r.speciesAmounts[i].amount < spec.specList[i].target)
                        r.objective.missing++;

                    /** Species missing **/
                } /** checking to see who I am missing **/

                /*
                if (bestscore)
                    ShowGenProgInfo("P.U. %i score %.6f Cost %.1f Connection %.1f Missing %i Amount %.1f \n",
                                    bestpu, bestscore, reserve->cost, reserve->connection, reserve->missing,
                                    reserve->penalty);
                */
            }     /** Add Pu as long as I've found one **/

        } while (bestscore);

        r.objective.total = r.objective.cost + r.objective.connection + r.objective.penalty;
    }

    vector<double> rare;

    private:
    /***** Product Irreplaceability for a single site ****/
    double ProdIrr(int ipu, Reserve& r, Pu& pu, Species& spec, Zones& zones)
    {
        int i, ism, isp;
        double product = 1;

        if (pu.puList[ipu].richness)
            for (i = 0; i < pu.puList[ipu].richness; i++)
            {
                ism = pu.puList[ipu].offset + i;
                isp = pu.puvspr[ism].spindex;
                if (pu.puvspr[ism].amount && (spec.specList[isp].target - r.speciesAmounts[isp].amount) > 0)
                    product *= (1 - Irreplaceability(ipu, isp, r, pu, spec));
            }

        return (1 - product);
    } /* Product Irreplaceability */

    /***** Sum Irreplaceability for a single site *****/
    double SumIrr(int ipu, Reserve& r, Pu& pu, Species& spec, Zones& zones)
    {
        int ism, isp;
        double sum = 0;

        if (pu.puList[ipu].richness)
            for (int i = 0; i < pu.puList[ipu].richness; i++)
            {
                ism = pu.puList[ipu].offset + i;
                isp = pu.puvspr[ism].spindex;
                if (pu.puvspr[ism].amount && (spec.specList[isp].target - r.speciesAmounts[isp].amount) > 0)
                    sum += (Irreplaceability(ipu, isp, r, pu, spec));
            }

        return (sum);
    } /* Sum Irreplaceability */

    /***** Irreplaceability For site for species *****/
    double Irreplaceability(int ipu, int isp, Reserve& r, Pu& pu, Species& spec)
    {
        double buffer, effamount;

        buffer = rare[isp] < spec.specList[isp].target ? 0 : rare[isp] - spec.specList[isp].target;
        if (r.speciesAmounts[isp].amount > spec.specList[isp].target)
            return (0);
        effamount = pu.RtnAmountSpecAtPu(ipu, isp);
        return (buffer < effamount ? 1 : effamount / buffer);
    }

    /* * * * * Greedy Score an alternative to the normal objective function *****/
    double GreedyScore(int ipu, Reserve& r, Pu& pu, Species& spec, Zones& zones, double blm)
    {
        double currpen, currcost, currscore;

        currpen = r.GreedyPen(ipu, pu, spec);
        currcost = pu.puList[ipu].cost + zones.ConnectionCost2(pu, ipu, 1, r.solution, r.solution[ipu], blm);
        if (currcost <= 0)
        {
            currscore = -1.0 / delta;
        } /* otherwise this 'free pu' will have +score */
        else
        {
            currscore = currpen / currcost;
        }

        return (currscore);
    } /* Score for a planning unit based upon greedy algorithm */

    // Supports heuristicMode 2,3,4,5
    double RareScoreByMode(int ipu, Reserve& r, Pu& pu, Species& spec, Zones& zones, double blm)
    {
        int ism, isp;
        int rareno = heuristicMode != 4 ? -1 : 0; // for mode 4, rareno starts from 0
        double rarest = 0;
        double rarescore = 0;

        for (int i = 0; i < pu.puList[ipu].richness; i++)
        {
            ism = pu.puList[ipu].offset + i;
            isp = pu.puvspr[ism].spindex;
            if (pu.puvspr[ism].amount && (spec.specList[isp].target > r.speciesAmounts[isp].amount))
            {
                if (heuristicMode == 2) //Max Rare Score Heuristic. PU scores based on rarest beast on PU
                {
                    if (1.0 / rare[isp] < rarest || rareno < 0)
                    {
                        rareno = isp;
                        rarest = rare[isp];
                    } /* Determine which is the rarest species */
                }
                else if (heuristicMode == 3) //Best Rarity Score. Determines each species rare score
                {
                    rarescore = RareScore(ipu, rareno, r, pu, spec, zones, blm) / rare[isp];
                    if (rarescore > rarest || rareno < 0)
                    {
                        rarest = rarescore;
                        rareno = isp;
                    }
                }
                else // Average(4) & Sum(5) Rare Score. Rare Score for each scoring species/number scoring species
                {
                    rarescore += RareScore(ipu, rareno, r, pu, spec, zones, blm) / rare[isp];
                    rareno++;
                }
            }
        }

        if (heuristicMode == 2)
        {
            if (rareno > -1)
                rarescore = RareScore(ipu, rareno, r, pu, spec, zones, blm) / rarest;
            else
                rarescore = 1.0 / delta;
        }
        else if (heuristicMode == 4) {
            return(rarescore/rareno);
        }

        return (rarescore);
    } /* Max Rare Score */

    /**** RareScore The score for a particular conservation value on a particular PU */
    double RareScore(int ipu, int isp, Reserve& r, Pu& pu, Species& spec, Zones& zones, double blm)
    {
        double currpen, currcost, currscore;
        double fold, newamount;

        fold = (spec.specList[isp].target - r.speciesAmounts[isp].amount);
        if (fold > 0)
        {
            if (spec.specList[isp].target2)
            {
                newamount = spec.specList[isp].target - r.speciesAmounts[isp].amount - r.PartialPen4(spec.specList[isp].target2, pu.RtnAmountSpecAtPu(ipu,isp));
            }
            else
                newamount = spec.specList[isp].target - r.speciesAmounts[isp].amount - pu.RtnAmountSpecAtPu(ipu,isp);

            currpen = newamount - fold;
        } /* Add new penalty if species isn't already in the system */

        currcost = pu.puList[ipu].cost + zones.ConnectionCost2(pu, ipu, 1, r.solution, r.solution[ipu], blm);
        if (currcost <= 0)
        {
            currscore = -1.0 / delta;
        } /* otherwise this 'free pu' will have +score */
        else
        {
            currscore = currpen / currcost;
        }

        return (currscore);
    }

    /****** Set Abundances ******/
    // Which are just the sum of all species in pu
    void SetAbundance(Pu& pu, Species& spec)
    {
        rare = pu.TotalSpeciesAmount(spec.spno);
    } /* Set Abundance */

    /* * * * *** Rarity Settup. Sets up rare score for each species ******/
    /**** score is total abundance / smallest species abundance * * * * */
    void SetRareness(Reserve& r, Pu& pu, Species& spec) {
        double smallest = 0;
        vector<double> fcount(spec.spno, 0);

        fcount = pu.TotalSpeciesAmountByAvailable(spec.spno, r.solution);

        for (int isp = 0; isp < spec.spno; isp++)
        {
            if (smallest == 0 || (fcount[isp] < smallest && fcount[isp] > 0))
                smallest = fcount[isp];
            rare[isp] = fcount[isp];
        }

        /*
        if (smallest == 0)
            ShowErrorMessage("Serious Error in calculating Rarenesses. No species detected.\n");
        */

        for (int isp = 0; isp < spec.spno; isp++)
            rare[isp] /= smallest;
    }

    mt19937 &rngEngine;
    int heuristicMode;

};

} // namespace marzone