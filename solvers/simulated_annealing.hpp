#include <limits>
#include <random>

#include "../marzone.hpp"
#include "../pu.hpp"
#include "../reserve.hpp"
#include "../zones.hpp"

namespace marzone {
using namespace std;

class SimulatedAnnealing {
    public:
    SimulatedAnnealing(int annealingOn, sanneal& anneal, mt19937& rngEngine) : rngEngine(rngEngine) {
        if (annealingOn)
        {
            settings = anneal; // copy construct
        }
    }

    void Initialize() {
        if (settings.type >= 2)
        {
            if (settings.type == 2)
            {
            }
            else if (settings.type == 3)
            {
            }
        }
    }
    
    void DumpBuffer() {

    }

    // Todo make private
    private:
    void ConnollyInit(Reserve& r, Species& spec, Pu& pu, Zones& zones) {
        double localdelta = numeric_limits<double>::epsilon();
        uniform_int_distribution<int> random_dist(0, numeric_limits<int>::max());
        uniform_int_distribution<int> random_pu_dist(0, pu.puno-1);

        // Set reserve to a random and evaluate 
        r.EvaluateObjectiveValue(pu, spec, zones);

        int ipu, iPreviousR, chosenZoneInd, chosenZone, imode = 1;
        for (int i=1;i<= settings.iterations/100; i++)
        {
            // pick pu at random
            do {
                ipu = random_pu_dist(rngEngine);
            } while ((pu.puList[ipu].status > 1) || (pu.puList[ipu].fPULock == 1)); // ignore locked/ non avail pu

            if (pu.puList[ipu].numZones > 1)
            {
                // pick a random available zone for this pu that is different to the current zone
                chosenZoneInd = random_dist(rngEngine) % pu.puList[ipu].numZones;
                chosenZone = pu.puZone[ipu][chosenZoneInd]-1;

                // if zone is already chosen, just increment it
                if (chosenZone == iPreviousR)
                {
                    chosenZoneInd = (chosenZoneInd + 1)%pu.puList[ipu].numZones;
                    chosenZone = pu.puZone[ipu][chosenZoneInd]-1;
                }
            }
            else if (pu.puList[ipu].numZones == 0) {
                // pu can be in  any zone.
                chosenZone = random_dist(rngEngine) % zones.zoneCount;
                if (chosenZone == iPreviousR)
                {
                    chosenZone = (chosenZone + 1)%zones.zoneCount;
                }
            }

            // Check change in zone

            // acceptance criteria

        }
    }

    void AdaptiveInit() {
        
    }

    sanneal settings;
    string debugbuffer;
    mt19937 &rngEngine;
};

} // namespace marzone