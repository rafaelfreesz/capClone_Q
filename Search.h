//
// Created by rafael on 03/08/22.
//

#ifndef CAPCLONE_SEARCH_H
#define CAPCLONE_SEARCH_H


#include "Antibody.h"
#include "Config.h"
#include "Stats.h"
#include "Utils.h"
#define NONE (-1)
#define NBR_SWAP 0
#define NONBR_SWAP 1
#define OPPS_SWAP 2
#define FIRST_ANTIBODY 0

class Search {
public:
    Search(Config *config, Instance *instance, double litSol, bool is_q, bool is_debug, int execNum);
    ~Search();

    void evolve();
    void evolve_debug_mode();
    void evolve_q();
    void evolve_q_debug_mode();
    void buildInitialPopulation();
    void operate();
    void operate_q();
    void maturate(int iArray, int cloneQty);
    void reselect();
    void regenerate();
    void improveMemory();
    void improveMemory_q();

    //Local Searches
    void vns(Antibody* antibody);
    void vns_q(Antibody* antibody);
    bool neighborsSwap(Antibody *antibody);
    bool nonNeiborhsSwap(Antibody *antibody);
    bool opositeSideSwap(Antibody *antibody);

    void print_q();
    void select_local_search();
    void calculate_reward(double delta);

    void printPopulation();
    void printClones();
    void printAll();
    void testAllPopulation();

    //Atributos;
    Config* config;
    Antibody** population;
    Instance* instance;
    double litSol;

    //Atributos QLearning
    int state;
    int action;
    double** q;
    bool* usar_action;
    double espacoAm;
    double last_improvement_delta;
    int execNum;
    bool is_debug;


private:
    void buildAntibody(int index);
    static bool antibodyCriterion(Antibody* a, Antibody* b);
    void sortPopulation();
    void sortClones();
    void swapAntibody(int i, int j);
    int* contaEscolha;
    Stats* myStatsPop;
    Stats* myStatsEscolhas;
    Stats* myStatsEscolha;
    bool is_q;
    int ls_exec_count;

};


#endif //CAPCLONE_SEARCH_H
