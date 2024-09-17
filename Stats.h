//
// Created by rafael on 09/08/22.
//

#ifndef CAPCLONE_STATS_H
#define CAPCLONE_STATS_H
#include <fstream>
#include <vector>
#include "Instance.h"
#include "Config.h"

using namespace std;
class Stats {
public:
    Stats(int intancesQty, Config *config);

    ~Stats();

    void setStat(int execI, int instanceI, double time, double cost);
    double getTime(int execI, int instanceI);
    double getCost(int execI, int instanceI);
    void printStats(string instanceName, int instanceI);
    void printHeader(Config* config);
    void getLitSolutions();

    int instancesQty;
    int execs;
    double* litSol;
    double* gapsSol;

    double** times;
    double** costs;

    double* bestTimes;
    double* bestCosts;

    double* avgTimes;
    double* avgCosts;

    ofstream statsFile;
};


#endif //CAPCLONE_STATS_H
