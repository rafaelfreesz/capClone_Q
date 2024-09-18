//
// Created by rafael on 03/08/22.
//

#ifndef CAPCLONE_SEARCH_H
#define CAPCLONE_SEARCH_H


#include "Antibody.h"
#include "Config.h"
#include "Utils.h"
#define NBR_SWAP 0
#define NONBR_SWAP 1
#define OPPS_SWAP 2
#define FIRST_ANTIBODY 0

class Search {
public:
    Search(Config *config, Instance *instance, double litSol);
    ~Search();

    void evolve();
    void evolve_q();
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
    void rvnd(Antibody* antibody);
    void rvnd_q(Antibody* antibody);
    bool neighborsSwap(Antibody *antibody);
    bool nonNeiborhsSwap(Antibody *antibody);
    bool opositeSideSwap(Antibody *antibody);


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
    double* Q;


private:
    void buildAntibody(int index);
    static bool antibodyCriterion(Antibody* a, Antibody* b);
    void sortPopulation();
    void sortClones();
    void swapAntibody(int i, int j);

};


#endif //CAPCLONE_SEARCH_H
