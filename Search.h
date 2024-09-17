//
// Created by rafael on 03/08/22.
//

#ifndef CAPCLONE_SEARCH_H
#define CAPCLONE_SEARCH_H


#include "Antibody.h"
#include "Config.h"
#include "Utils.h"

#define FIRST_ANTIBODY 0

class Search {
public:
    Search(Config *config, Instance *instance, double litSol);
    ~Search();

    void evolve();
    void buildInitialPopulation();
    void operate();
    void maturate(int iArray, int cloneQty);
    void reselect();
    void regenerate();
    void improveMemory();

    //Local Searches
    void localSearch(Antibody *antibody);
    void rvnd(Antibody* antibody);
    bool neighborsSwap(Antibody *antibody);
    bool nonNeiborhsSwap(Antibody *antibody);
    bool opositeSideSwap(Antibody *antibody);


    void PathR(Antibody *antibody);


    void printPopulation();
    void printClones();
    void printAll();
    void testAllPopulation();

    //Atributos;
    Config* config;
    Antibody** population;
    Instance* instance;
    double litSol;

private:
    void buildAntibody(int index);
    static bool antibodyCriterion(Antibody* a, Antibody* b);
    void sortPopulation();
    void sortClones();
    void swapAntibody(int i, int j);

};


#endif //CAPCLONE_SEARCH_H
