//
// Created by rafael on 09/08/22.
//

#include "Stats.h"

Stats::Stats(int intancesQty, Config *config, bool isQL) {
    this->execs=config->executions;
    this->instancesQty=intancesQty;
    this->isQL = isQL;

    this->bestTimes=new double [intancesQty];
    this->bestCosts=new double [intancesQty];
    this->avgTimes=new double [intancesQty];
    this->avgCosts=new double [intancesQty];
    this->times=new double * [intancesQty];
    this->costs=new double * [intancesQty];
    this->costs=new double * [intancesQty];
    this->litSol=new double [intancesQty];
    this->gapsSol=new double [intancesQty];

    for(int i=0;i<intancesQty;i++){
        this->times[i]=new double [execs];
        this->costs[i]=new double [execs];
        this->bestTimes[i]=0.0;
        this->bestCosts[i]=0.0;
        this->avgTimes[i]=0.0;
        this->avgCosts[i]=0.0;
    }

    getLitSolutions();
    printHeader(config);



}

Stats::Stats(string type, bool isQL){
    if(type == "print_pop") {
        printHeaderPrintPop();
    }else if(type == "escolha") {
        printHeaderPrintEscolha();
    }else {
        printHeaderPrintEscolhas();
    }

    this->execs=30;
    this->instancesQty=2;
    this->isQL = isQL;

    this->bestTimes=new double [2];
    this->bestCosts=new double [2];
    this->avgTimes=new double [2];
    this->avgCosts=new double [2];
    this->times=new double * [2];
    this->costs=new double * [2];
    this->costs=new double * [2];
    this->litSol=new double [2];
    this->gapsSol=new double [2];

    for(int i=0;i<2;i++){
        this->times[i]=new double [execs];
        this->costs[i]=new double [execs];
        this->bestTimes[i]=0.0;
        this->bestCosts[i]=0.0;
        this->avgTimes[i]=0.0;
        this->avgCosts[i]=0.0;
    }

}

Stats::~Stats() {

    this->statsFile.close();
    
    for(int i=0;i<this->instancesQty;i++){
        if(this->times[i]!=nullptr) {
            delete []this->times[i];
        }
        if(this->costs[i]!=nullptr) {
            delete []this->costs[i];
        }
    }

    if(this->times!=nullptr) {
        delete []this->times;
    }
    if(this->costs!=nullptr) {
        delete []this->costs;
    }
    if(this->bestTimes!=nullptr) {
        delete []this->bestTimes;
    }
    if(this->bestCosts!=nullptr) {
        delete []this->bestCosts;
    }
    if(this->avgTimes!=nullptr) {
        delete []this->avgTimes;
    }
    if(this->avgCosts!=nullptr) {
        delete []this->avgCosts;
    }
    if(this->litSol!=nullptr) {
        delete []this->litSol;
    }
    if(this->gapsSol!=nullptr) {
        delete []this->gapsSol;
    }

}

void Stats::setStat(int execI, int instanceI, double time, double cost) {
    this->times[instanceI][execI]=time;
    this->costs[instanceI][execI]=cost;

    this->avgTimes[instanceI]+=time;
    this->avgCosts[instanceI]+=cost;

    if(execI==0 || time<this->bestTimes[instanceI]){
        this->bestTimes[instanceI]=time;
    }

    if(execI==0 || cost<this->bestCosts[instanceI]){
        this->bestCosts[instanceI]=cost;
    }

    if(execI==(this->execs-1)){
        this->avgTimes[instanceI]/=this->execs;
        this->avgCosts[instanceI]/=this->execs;
    }

}

double Stats::getTime(int execI, int instanceI) {
    return this->times[instanceI][execI];
}

double Stats::getCost(int execI, int instanceI) {
    return this->costs[instanceI][execI];
}

void Stats::printStats(string instanceName, int instanceI) {

    this->gapsSol[instanceI]=(this->bestCosts[instanceI]-this->litSol[instanceI])*100/this->litSol[instanceI];

    this->statsFile<< instanceName +
        ";" + to_string(this->avgTimes[instanceI]) +
        ";" + to_string(this->avgCosts[instanceI]) +
        ";" + to_string(this->bestTimes[instanceI]) +
        ";"+to_string(this->bestCosts[instanceI]) +
        ";"+to_string(this->gapsSol[instanceI])<<endl;


}

void Stats::printHeader(Config *config) {

    int fileId=0;
    string fileName="Stats//stats_"+ to_string(fileId);

    if(this->isQL){
        fileName = "Stats//stats_Q_"+ to_string(fileId);
    }

    this->statsFile.open(fileName,ofstream::in);
    while(this->statsFile.is_open()){
        this->statsFile.close();
        fileId++;
        fileName="Stats//stats_"+ to_string(fileId);
        if(this->isQL){
            fileName = "Stats//stats_Q_"+ to_string(fileId);
        }
        this->statsFile.open(fileName,ofstream::in);
    }
    this->statsFile.close();
    this->statsFile.open(fileName);

    this->statsFile<<"pSize: "<<config->pSize<<endl;
    this->statsFile<<"gen: "<<config->gen<<endl;
    this->statsFile << "memorySetSizen: " << config->memorySetSize << endl;
    this->statsFile<<"betaCoeff: "<<config->betaCoeff<<endl;
    this->statsFile<<"clonePop: "<<config->clonePop<<endl;
    this->statsFile<<"arraySize: "<<config->arraySize<<endl;
    this->statsFile << "regQty: " << config->regenerationQty << endl;
    this->statsFile<<"seed: \n\t";

    for(int i=0;i<config->executions; i++){
        this->statsFile<<config->seeds[i]<<" ";
    }

    this->statsFile<<endl<<"clonesPerI: \n\t";
    for(int i=0;i<config->memorySetSize; i++){
        this->statsFile<<config->clonesPerI[i]<<" ";
    }
    this->statsFile<<endl<<endl<<"Instancia;avgTime;avgCost;bestTime;bestCost"<<endl;
    
   
}

void Stats::printHeaderPrintPop() {
    int fileId=0;
    string fileName="Stats//stats_det_pop"+ to_string(fileId);

    if(this->isQL){
        fileName = "Stats//stats_Q_det_pop"+ to_string(fileId);
    }

    this->statsFile.open(fileName,ofstream::in);

    if(!this->statsFile.is_open()){

        string header = "Instancia,execucao,geracao";


        for(int i=0;i<10;i++) {
            header+= ",i_"+to_string(i+1);
        }
        this->statsFile<<header<<endl;

    }
    this->statsFile.close();
    this->statsFile.open(fileName,ofstream::app);
}

void Stats::printHeaderPrintEscolhas() {
    int fileId=0;
    string fileName="Stats//stats_det_esc"+ to_string(fileId);

    if(this->isQL){
        fileName = "Stats//stats_Q_det_esc"+ to_string(fileId);
    }


    this->statsFile.open(fileName,ofstream::in);

    if(!this->statsFile.is_open()){

        string header = "Instancia,execucao,ls1,ls2,ls3";
        this->statsFile<<header<<endl;

    }

    this->statsFile.close();
    this->statsFile.open(fileName,ofstream::app);
}

void Stats::printHeaderPrintEscolha() {
    int fileId=0;
    string fileName="Stats//stats_det_esc_det"+ to_string(fileId);

    if(this->isQL){
        fileName = "Stats//stats_Q_det_esc_det"+ to_string(fileId);
    }


    this->statsFile.open(fileName,ofstream::in);

    if(!this->statsFile.is_open()){

        string header = "instancia,exec,ls_exec,ls1,ls2,ls3";
        this->statsFile<<header<<endl;

    }

    this->statsFile.close();
    this->statsFile.open(fileName,ofstream::app);
    if(!this->statsFile.is_open()) {
        exit(100);
    }
}

void Stats::getLitSolutions() {
    ifstream litSolFile;
    litSolFile.open("Stats/litSol");

    string line;

    for(int i=0;i<this->instancesQty;i++){
        getline(litSolFile,line);
        this->litSol[i]=stod(line);
    }
    litSolFile.close();
}

void Stats::printPop(Antibody **a, int popsize, string instName, int execnum, int gerNum) {
    string line = instName+","+to_string(execnum)+","+to_string(gerNum);

    for(int i=0; i<popsize;i++) {
        line+=","+to_string(a[i]->cost);
    }

    this->statsFile<<line<<endl;
}

void Stats::printEscolhas(int *vetor, int execnum, string instname) {
    string line = instname+","+to_string(execnum)+","+to_string(vetor[0])+","+to_string(vetor[1])+","+to_string(vetor[2]);

    this->statsFile<<line<<endl;
}

void Stats::printEscolha(int ls, int ls_exec_count, int execnum, string instname) {

    int lss[3] = {0,0,0};
    lss[ls] = 1;

    string line = instname + "," + to_string(execnum) + "," + to_string(ls_exec_count);
    for(int i=0 ;i<3;i++){ line+= "," + to_string(lss[i]);}
    this->statsFile<<line<<endl;
}

