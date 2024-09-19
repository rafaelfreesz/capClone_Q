//
// Created by rafael on 09/08/22.
//

#ifndef CAPCLONE_CONFIG_H
#define CAPCLONE_CONFIG_H


#include <cmath>
#include <ctime>

class Config {

public:

    Config(int pSize, int gen, double memorySetRate, double betaCoeff, double regenerationRate, int execs, double alpha, double epsilon) {
        this->pSize=pSize;
        this->gen=gen;
        this->memorySetSize=(int)(memorySetRate * pSize);
        this->betaCoeff=betaCoeff;
        this->regenerationQty= this->pSize * regenerationRate;
        this->executions=execs;
        this->clonesPerI= nullptr;
        this->alpha = alpha;
        this->epsilon = epsilon;

        srand(clock());
        this->seeds=new long [execs];
        for(int i=0;i<execs;i++){
            this->seeds[i]=rand();
        }


        countClonalSelection();
        this->arraySize= this->pSize + this->clonePop;

    };

    ~Config() {
        delete[] this->clonesPerI;
        delete[] this->seeds;
    }

    void countClonalSelection(){

        vector<int> clonesPerIVector;
        this->memorySetSize=0;
        int clones=this->pSize;
        int pos=1;
        double div;

        while(clones>0 && pos<this->pSize){
            clonesPerIVector.push_back(clones);
            div=(double)clones/++pos;
            clones= (int)round(div);
        }

        this->memorySetSize=(int)clonesPerIVector.size();
        this->clonesPerI=new int [this->memorySetSize];

        this->clonePop=0;

        for(int i=0;i<clonesPerIVector.size();i++){
            this->clonePop+=clonesPerIVector.at(i);
            this->clonesPerI[i]=clonesPerIVector.at(i);
        }

        this->betaCoeff=(double)this->clonePop/this->pSize;

    }

    void print(){
        cout<<"pSize: "<<pSize<<endl;
        cout<<"gen: "<<gen<<endl;
        cout << "nClonalSelection: " << memorySetSize << endl;
        cout<<"betaCoeff: "<<betaCoeff<<endl;
        cout<<"clonePop: "<<clonePop<<endl;
        cout<<"arraySize: "<<arraySize<<endl;
        cout << "regQty: " << regenerationQty << endl;
        cout<<"seeds: "<<endl<<"\t";
        for(int i=0; i < executions; i++){
            cout<<this->seeds[i]<<" ";
        }
        cout<<"clonesPerI: ";
        for(int i=0; i < memorySetSize; i++){
            cout<<this->clonesPerI[i]<<" ";
        }
        cout<<endl;

    }

    int pSize;
    int gen;
    int memorySetSize;
    double betaCoeff;
    int clonePop;
    int arraySize;
    int regenerationQty;
    int executions;
    long* seeds;
    int* clonesPerI;
    double alpha;
    double epsilon;


};


#endif //CAPCLONE_CONFIG_H
