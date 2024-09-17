//
// Created by rafael on 26/07/22.
//

#ifndef CAPCLONE_ANTIBODY_H
#define CAPCLONE_ANTIBODY_H


#include "Instance.h"

class Antibody {
public:

    Antibody(Instance* instance);
    ~Antibody();

    void swapFacility(int i, int j);
    void shake(int size);
    void adjustP();

    Antibody* clone();
    void clone(Antibody* targetClone);

    void calculateAbcissa();
    void calculatCost();
    void calculateSolution();

    //Swap calculating
    void calculateSwap(int i, int j);

    void neigbohrCalculation(int iMin, int iMax);
    void nonNeigbohrCalculation(int iMin, int iMax);
    void opositeSideCalculation(int iMin, int iMax);

    void testCalculation();

    void print();

    int* layout;
    double* abcissa;
    int p;
    double cost;
    bool improved;
    Instance* instance;

    //calculation testing fucntions
    void sameSideCalc(int iMin, int iMax);
};


#endif //CAPCLONE_ANTIBODY_H
