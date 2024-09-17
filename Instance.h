//
// Created by rafael on 19/07/22.
//

#ifndef CAPCLONE_INSTANCE_H
#define CAPCLONE_INSTANCE_H

#include <iostream>
#include <vector>

using namespace std;
class Instance {

public:

    Instance(int n);
    ~Instance();

    void calculateLayoutLength();

    void print();
    void verify();


    string name;
    int n;
    double* lengths;
    double* halfLengths;
    int** demands;
    double layoutLengh;

};


#endif //CAPCLONE_INSTANCE_H
