//
// Created by rafael on 09/08/22.
//

#ifndef CAPCLONE_STATS_H
#define CAPCLONE_STATS_H
#include <fstream>
#include <vector>
#include "Instance.h"
#include "Config.h"
#include "Antibody.h"

using namespace std;
class Stats {
public:
    Stats(int intancesQty, Config *config, bool isQL);
    Stats(string type, bool isQL);

    ~Stats();

    void setStat(int execI, int instanceI, double time, double cost);
    double getTime(int execI, int instanceI);
    double getCost(int execI, int instanceI);
    void printStats(string instanceName, int instanceI);
    void printHeader(Config* config);
    void printHeaderPrintPop();
    void printHeaderPrintEscolhas();
    void printHeaderPrintEscolha();
    void getLitSolutions();
    void printPop(Antibody** a, int popsize, string instName, int execnum, int gerNum);
    void printEscolhas(int* vetor, int execnum, string instname);
    void printEscolha(int ls, int ls_exec_count, int execnum, string instname);
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

    bool isQL;
};


#endif //CAPCLONE_STATS_H
