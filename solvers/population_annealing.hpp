#pragma once

/*
    Population annealing is an extension of the simulated annealing process.
    The implementation roughly follows the paper described here, with some adjustments: https://arxiv.org/pdf/2102.06611.pdf
    Uses the same parameters as simulated annealing, although with some different interpretations.

    Currently this is experimental code. It is slower than traditional simulated annealing but if left running long enough, can get better results.
*/

#include <limits>
#include <random>

#include "../common.hpp"
#include "../pu.hpp"
#include "../reserve.hpp"
#include "../zones.hpp"
#include "../logger.hpp"

namespace marzone {
using namespace std;

typedef struct scosttrace : scost {
    double temperature;
    unsigned popSize;
} scosttrace;

class PopulationAnnealing {
    public:
    PopulationAnnealing(sanneal& anneal, mt19937& rngEngine, int id, sfname& fnames) 
    : rngEngine(rngEngine), id(id), numSweeps(0), tIterations(0)
    {
        populationSize = max(omp_get_max_threads(), 10);
        settings = anneal;
        settings.type = 2; // only type 2 thermal annealing supported - override others.
        outputPrefix = fnames.savename;
        saveAnnealingTrace = fnames.saveannealingtrace;
    }

    void Run(Reserve& r, Species& spec, Pu& pu, Zones& zones, double tpf1, double tpf2, double costthresh, double blm, Logger& logger) {
        uniform_int_distribution<int> randomDist(0, numeric_limits<int>::max());
        uniform_real_distribution<double> floatRange(0.0, 1.0);
        
        vector<Reserve> population;
        vector<schange> changePop;
        InitPopulation(population, changePop, r, spec, pu, zones, costthresh, blm, tpf1, tpf2);

        logger.ShowWarningMessage("Temp start " + to_string(settings.temp) + " dec: " + to_string(settings.Tcool) + "\n");

        if (saveAnnealingTrace)
        {
            iterationResults.push_back({r.objective, settings.temp, populationSize});
        }
            
        for (unsigned i = 0; i < tIterations; i++)
        {
            // Sweep across all population reserves
            #pragma omp parallel for
            for (unsigned j = 0; j < populationSize; j++) {
                for (unsigned k = 0; k < numSweeps; k++)
                    Sweep(population[j], changePop[j], spec, pu, zones, tpf1, tpf2, costthresh, blm, randomDist, floatRange);
            }

            // store best solution if trace required.
            if (saveAnnealingTrace)
            {
                scost& best = population[FindBestReserve(population)].objective;
                iterationResults.push_back({best, settings.temp, populationSize});
            }
            
            // Resample population.
            ResamplePopulation(population, changePop, spec, zones);

            // Technically we should keep track of the bestR after every sweep.
            // update temperatures. Reduce existing temperature
            settings.temp = settings.temp*settings.Tcool;
        }

        // Write pop trace if needed
        if (saveAnnealingTrace)
            WriteTrace();

        // Replace r with the best reserve solution
        r.Assign(population[FindBestReserve(population)]);
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

    unsigned FindBestReserve(vector<Reserve>& population) {
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

        return bestInd;
    }

    mt19937 &rngEngine;
    int id;
    sanneal settings;

    private:
    unsigned populationSize; // starting population size.
    const unsigned minPopulationSize = omp_get_max_threads(); // minimum population size (such that resampling cant go below this number)
    uint64_t numSweeps;
    uint64_t tIterations; // number of temperature changes.

    // for convergence plotting.
    int saveAnnealingTrace; 
    vector<scosttrace> iterationResults;
    string outputPrefix;

    void WriteTrace() {
        string filename = outputPrefix + "_popanneal_trace" + intToPaddedString(id, 5) + getFileSuffix(saveAnnealingTrace);

        ofstream file;
        file.open(filename);
        file << "t,iteration,total,cost,connection,penalty,shortfall,populationsize\n"; // header

        for (int i = 1; i <= iterationResults.size(); i++) 
        {
            scosttrace& term = iterationResults[i-1];
            file << term.temperature << "," << i << "," << term.total << "," << term.cost << "," 
                << term.connection << "," << term.penalty << "," << term.shortfall << "," << term.popSize << "\n";
        }

        file.close();
    }

    // based on reserve values, resamples and population and duplicates/removes certain pops.
    void ResamplePopulation(vector<Reserve>& population, vector<schange>& changePop, Species& spec, Zones& zones)
    {
        // get new weights for each replica
        vector<unsigned> newWeights = CalculateResamplingWeights(population);
        vector<Reserve> oldPop = population; // copy construct old population.
        unsigned newPopSize = 0;

        population.clear();
        for (int i =0; i < populationSize; i++)
        {
            if (newWeights[i] > 0) 
            {
                newPopSize += newWeights[i];
                // clone this many replicas.
                for (int j = 0; j < newWeights[i]; j++)
                {
                    Reserve& rTemp(oldPop[i]);
                    population.push_back(rTemp);
                }
            }
        }

        // modify changepop if the new pop size is larger than before.
        for (int i = changePop.size(); i < newPopSize; i++) 
        {
            changePop.push_back(population[0].InitializeChange(spec, zones));
        }
        
        populationSize = newPopSize;
    }

    // gets the boltzmann resampling weights for each population.
    vector<unsigned> CalculateResamplingWeights(vector<Reserve>& population) {
        double betaDiff = 1/(settings.temp) - 1/(settings.temp*settings.Tcool);
        double objectiveScale = GetObjectiveScale(population, betaDiff); // normalization factor since we run into values that are too close to 0, or too high.
        vector<double> newWeights(populationSize, 0);
        vector<unsigned> newPopulations(populationSize, 1);

        // compute average ratio
        double averageRatio = 0.0;
        for (int i =0; i < populationSize; i++)
        {
            averageRatio += exp(objectiveScale*betaDiff*population[i].objective.total);
        }
        averageRatio /= populationSize;

        // compute new weights and scale to the min population size the new # of copies of this configuration.
        int newPopSampling = 1;
        for (int i =0; i < populationSize; i++)
        {
            newWeights[i] = exp(objectiveScale*betaDiff*population[i].objective.total)/averageRatio;
            newPopSampling += (unsigned) newWeights[i];
        }

        // if new pop is 0 for whatever reason, we just keep the existing pop
        if (newPopSampling == 0) 
        {
            return newPopulations;
        }

        // The idea behind scaleUpFactor is to keep new pop size above the minimum pop size.
        double scaleUpFactor = newPopSampling > minPopulationSize ? 1 : minPopulationSize/(double) newPopSampling;
        unsigned total = 0;
        for (int i =0; i < populationSize; i++)
        {
            newPopulations[i] = (unsigned) (newWeights[i]*scaleUpFactor);
            total += newPopulations[i];
        }

        return newPopulations;
    }

    // Given a set of population, get the average normalizing factor to deal with values that are too small or big compared to the beta. 
    // betaDiff is a negative value.
    double GetObjectiveScale(vector<Reserve>& population, double betaDiff) {
        double bestObj = numeric_limits<double>::max();;
        double betaMagnitude = -betaDiff;
        for (int i = 0; i < population.size(); i++) 
        {
            if (population[i].objective.total < bestObj)
                bestObj = population[i].objective.total;
        }

        bestObj *= betaMagnitude;
        return 100/bestObj;
    }

    // Initializes parameters of a population
    void InitPopulation(vector<Reserve>& population, vector<schange>& changePop, Reserve& r, Species& spec, Pu& pu, Zones& zones, 
    double costthresh, double blm, double tpf1, double tpf2) 
    {
        population.reserve(populationSize);
        changePop.resize(populationSize);

        sfname fnamesDummy = {}; // placeholder
        InitParameters(r, pu, spec, zones, costthresh, blm, tpf1, tpf2);

        // create reserve configurations and add them to population.
        for (int i = 0; i < populationSize; i++) {
            Reserve newR(r);
            population.push_back(newR);
            changePop[i] = r.InitializeChange(spec, zones);
        } 

        // set sweeps consistent with sa
        tIterations = settings.Titns;
        numSweeps = max((uint64_t) 1, (uint64_t) settings.iterations/tIterations/pu.puno);
    }

    void InitParameters(Reserve& rInit, Pu& pu, Species& spec, Zones& zones, double costthresh, double blm, double tpf1, double tpf2) 
    {
        Reserve r(rInit);
        // Sweep once and set temperatures to averageDelta and averageDelta/tItns
        uniform_int_distribution<int> randomDist(0, numeric_limits<int>::max());
        double averageDelta = 0.0, absChange = 0.0;
        schange change = r.InitializeChange(spec, zones);
        unsigned iZone;

        for (int& ipu: pu.validPuIndices) {
            iZone = pu.RtnValidZoneForPu(ipu, r.solution[ipu], randomDist, rngEngine, zones.zoneCount);
            r.CheckChangeValue(change, ipu, r.solution[ipu], iZone, pu, zones, spec, costthresh, blm, 
                tpf1, tpf2);
            
            if (change.total < 0)
                r.ApplyChange(ipu, iZone, change, pu, zones, spec);

            absChange = abs(change.total);
            averageDelta += absChange/(double) pu.validPuIndices.size();
        }

        settings.Tinit = averageDelta/(double) settings.Titns;
        settings.temp = settings.Tinit;
        settings.Tcool = pow(1.0/averageDelta, 1.0/(double)settings.Titns);
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
