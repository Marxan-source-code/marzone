#pragma once

#include "common.hpp"
#include "util.hpp"
#include "pu.hpp"
#include "zones.hpp"
#include "reserve.hpp"

#include <omp.h>

/*
  Deals with calculation of any output analysis of solutions, reserves and other.
*/

namespace marzone {
using namespace std;

class Analysis {
  public:
  Analysis() {
  }

  void initSumSolution(int puno, int zoneCount) {
    // initialize analysis objects.
    sumSoln.resize(puno, 0);
    zoneSumSoln.resize(puno*zoneCount, 0);

    // Init locks
    omp_init_lock(&sumLock);
  }

  // Writes both sumSoln and zoneSumSoln
  void WriteAllSumSoln(string filename, Pu& pu, Zones& zones, int imode) {
    ofstream myfile;
    myfile.open(filename); /* Imode = 1, REST output, Imode = 2, Arcview output */
    string d = imode > 1 ? "," : " ";

    if (imode > 1)
    {
      myfile << "\"planning unit\",\"number\"";
      for (int j = 0; j < zones.zoneCount; j++)
        myfile << ",\"" << zones.zoneNameIndexed[j].name << "\"";
      myfile << "\n";
    }

    for (int i = 0; i < pu.puno; i++)
    {
      if (imode > 1)
      {
        myfile << pu.puList[i].id << d << sumSoln[i];
        for (int j = 0; j < zones.zoneCount; j++)
          myfile << d << zoneSumSoln[pu.puno * j + i];
        myfile << "\n";
      }
    }

    myfile.close();
  }

  // Should be threadsafe. Applies the solution in Reserve to sums
  void ApplyReserveToSumSoln(Reserve& r) {
    int iZoneSumSolnIndex;
    int rSize = r.solution.size();

    omp_set_lock(&sumLock);
    for (int i = 0; i < rSize; i++) {
      if (r.solution[i] >= 0) {
        sumSoln[i] ++;
        iZoneSumSolnIndex = rSize * (r.solution[i]) + i;
        zoneSumSoln[iZoneSumSolnIndex]++;
      }
    }
    omp_unset_lock(&sumLock);
  }

  private:
  vector<int> sumSoln; // size puno
  vector<int> zoneSumSoln; //size puno*zoneCount

  // thread safety
  omp_lock_t sumLock;
};

} // namespace marzone
