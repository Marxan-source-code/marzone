#pragma once

/*
    Population annealing is an extension of the existing simulated annealing solver.
*/

#include <limits>
#include <random>

#include "../common.hpp"
#include "../pu.hpp"
#include "../reserve.hpp"
#include "../zones.hpp"
#include "../logger.hpp"
#include "simulated_annealing.hpp"

namespace marzone {
using namespace std;

class PopulationAnnealing {
    public:
    PopulationAnnealing(sanneal& anneal, mt19937& rngEngine, int id, LoggerBase& logger) 
    : rngEngine(rngEngine), id(id), populationSize(50), numSweeps(0)
    {
        settings = anneal;
        settings.type = 2; // only type 2 thermal annealing supported - override others.
    }

    void Run(Reserve& r, Species& spec, Pu& pu, Zones& zones, double tpf1, double tpf2, double costthresh, double blm) {
        uniform_int_distribution<int> randomDist(0, numeric_limits<int>::max());
        uniform_real_distribution<double> floatRange(0.0, 1.0);
        
        vector<Reserve> population;
        vector<schange> changePop;
        InitPopulation(population, changePop, r, spec, pu, zones, blm);

        logger.ShowWarningMessage("Running population annealing with " + to_string(numSweeps) + " sweeps.");
        for (int i = 0; i < numSweeps; i++)
        {
            // Sweep across all population reserves
            //#pragma omp parallel for schedule(dynamic)
            for (int j = 0; j < populationSize; j++) {
                Sweep(population[i], changePop[i], spec, pu, zones, tpf1, tpf2, costthresh, blm, randomDist, floatRange);
            }
            
            // Resample population.
            ResamplePopulation(population, changePop, spec, zones);

            // Todo - should keep track of the bestR after every sweep.
            // update temperatures. Linearly reduce existing temperature
            settings.temp = settings.temp*settings.Tcool;
        }

        // Replace r with the best reserve solution
        r = FindBestReserve(population);
    }

    // perform a single sweep on the reserve, and only on non-locked pu.
    void Sweep(Reserve& r, schange& change, Species& spec, Pu& pu, Zones& zones, double tpf1, double tpf2, 
        double costthresh, double blm, uniform_int_distribution<int>& randomDist, uniform_real_distribution<double>& floatRange)
    {
        int iZone, iPrevZone, zoneCount = zones.zoneCount;

        for (int ipu: pu.validPuIndices)
        {
            iPrevZone = r.solution[ipu];
            iZone = pu.RtnValidZoneForPu(ipu, iPrevZone, randomDist, rngEngine, zoneCount);
            r.CheckChangeValue(change, ipu, iPrevZone, iZone, pu, zones, spec, costthresh, blm, 
                tpf1, tpf2);
            
            if (GoodChange(change.total, floatRange))
            {
                r.ApplyChange(ipu, iZone, change, pu, zones, spec);
            }
        }
    }

    Reserve& FindBestReserve(vector<Reserve>& population) {
        int bestInd = 0;
        double bestValue = numeric_limits<double>::max();

        for (int i = 0; i < population.size(); i++)
        {
            if (population[i].objective.total < bestValue)
            {
                bestValue = population[i].objective.total;
                bestInd = i;
            }
        }

        return population[bestInd];
    }

    mt19937 &rngEngine;
    int id;
    sanneal settings;
    LoggerBase logger;

    private:
    unsigned populationSize; // starting population size.
    uint64_t numSweeps;

    // based on reserve values, resamples and population and duplicates/removes certain pops.
    void ResamplePopulation(vector<Reserve>& population, vector<schange>& changePop, Species& spec, Zones& zones)
    {
        // get new weights for each replica
        vector<unsigned> newWeights = CalculateResamplingWeights(population);
        vector<Reserve> newPop;
        unsigned newPopSize = 0;

        for (int i =0; i < populationSize; i++)
        {
            if (newWeights[i] > 0) 
            {
                newPopSize += newWeights[i];
                // clone this many replicas.
                for (int j = 0; j < newWeights[i]; j++)
                {
                    Reserve& rTemp(population[i]);
                    newPop.push_back(rTemp);
                }
            }
        }

        // modify changepop if the new pop size is larger than before.
        for (int i = changePop.size(); i < newPopSize; i++) 
        {
            changePop.push_back(population[0].InitializeChange(spec, zones));
        }
        
        populationSize = newPopSize;
        population = newPop;
    }

    // gets the boltzmann resampling weights for each population.
    vector<unsigned> CalculateResamplingWeights(vector<Reserve>& population) {
        double betaDiff = 1/(settings.temp*(1-settings.Tcool));
        vector<unsigned> newWeights(populationSize, 0);

        // compute average ratio
        double averageRatio = 0.0;
        for (int i =0; i < populationSize; i++)
        {
            averageRatio += exp(betaDiff*population[i].objective.total);
        }
        averageRatio /= populationSize;

        // compute new ratio and take the floor to be the new # of copies of this configuration.
        for (int i =0; i < populationSize; i++)
        {
            newWeights[i] = (unsigned) exp(betaDiff*population[i].objective.total)/averageRatio;
        }

        return newWeights;
    }

    // Initializes parameters of a population
    void InitPopulation(vector<Reserve>& population, vector<schange>& changePop, Reserve& r, Species& spec, Pu& pu, Zones& zones, double blm) 
    {
        population.reserve(populationSize);
        changePop.resize(populationSize);

        sfname fnamesDummy = {}; // placeholder
        SimulatedAnnealing base(fnamesDummy, logger, 1, settings, rngEngine, 0, id);
        base.Initialize(spec, pu, zones, 0, blm); // initialize parameters.
        settings = base.settings; // borrow the initialized temperature parameters from first run. 

        // create reserve configurations and add them to population.
        for (int i = 0; i < populationSize; i++) {
            Reserve newR(r);
            population.push_back(newR);
            changePop[i] = r.InitializeChange(spec, zones);
        } 

        // set sweeps consistent with sa
        numSweeps = settings.iterations/settings.Tlen;
    }

    bool GoodChange(double changeTotal, uniform_real_distribution<double>& float_range)
    {
        // avoid math operations if possible
        if (changeTotal < 0) {
            return true;
        }

        return (exp(-changeTotal / settings.temp) > float_range(rngEngine));
    }

};
} // namespace marzone
