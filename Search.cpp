//
// Created by rafael on 03/08/22.
//

#include <algorithm>
#include <random>
#include "Search.h"

Search::Search(Config *config, Instance *instance, double litSol) {
    this->config=config;
    this->population= new Antibody* [this->config->arraySize];
    this->instance=instance;
    this->litSol=litSol;

    this->q = new double*[3];
    for(int i=0;i<3;i++) {
        this->q[i] = new double[3];
        for(int j=0;j<3;j++) {
            if(i!=j) {
                this->q[i][j] = 1000.0;
            }else {
                this->q[i][j] = 0.0;
            }
        }
    }

    this->usar_action = new bool[3];
    this->usar_action[0] = true;
    this->usar_action[1] = true;
    this->usar_action[2] = true;

    this->espacoAm = (this->q[0][0]+this->q[0][1]+this->q[0][2]);
    this->action = -1;
    this->state = -1;
    this->last_improvement_delta = -1;
}

Search::~Search() {
    delete [] this->q;
    delete [] this->usar_action;
    if(this->population!= nullptr){
        for(int i=0;i< this->config->pSize;i++){
            if(this->population[i]!= nullptr) {
                delete this->population[i];
            }else{
                break;
            }
        }
        delete[] this->population;

    }

}

void Search::evolve() {
    buildInitialPopulation();

   for(int g=0;g<this->config->gen;g++){
        operate();
        reselect();
        regenerate();
   }
    improveMemory();
}

void Search::evolve_q() {
    buildInitialPopulation();

   for(int g=0;g<this->config->gen;g++){
        operate_q();
        reselect();
        regenerate();
   }
    improveMemory_q();
}

void Search::improveMemory() {
    for(int i=0;i<this->config->memorySetSize;i++){
        if(!this->population[i]->improved) {
            vns(this->population[i]);
            this->population[i]->improved= true;
        }
    }
    sortPopulation();
}

void Search::improveMemory_q() {
    for(int i=0;i<this->config->memorySetSize;i++){
        if(!this->population[i]->improved) {
            vns_q(this->population[i]);
            this->population[i]->improved= true;
        }
    }
    sortPopulation();
}
//Processo de clonagem
void Search::operate() {

    int iArray=this->config->pSize;

    for(int i=0;i<this->config->memorySetSize; i++){

        for(int j=0;j<this->config->clonesPerI[i];j++){
            this->population[i]->clone(this->population[iArray]);
            this->population[iArray]->improved= false;
            maturate(iArray++, this->config->clonesPerI[i]); //Maturação (mutação) do indivíduo

        }

    }

    sortClones();


    vns(this->population[this->config->pSize]);
    this->population[this->config->pSize]->improved= true;

}

void Search::operate_q() {

    int iArray=this->config->pSize;

    for(int i=0;i<this->config->memorySetSize; i++){

        for(int j=0;j<this->config->clonesPerI[i];j++){
            this->population[i]->clone(this->population[iArray]);
            //this->population[iArray]=this->population[i]->clone(); //Clonagem do anticorpo
            this->population[iArray]->improved= false;//Clonagem do anticorpo
            maturate(iArray++, this->config->clonesPerI[i]); //Maturação (mutação) do indivíduo

        }

    }

    sortClones();


    vns_q(this->population[this->config->pSize]);
    this->population[this->config->pSize]->improved= true;

}

//Processo de maturação (mutação)
void Search::maturate(int iArray, int cloneQty) {

    double probSwap = (1.0-(double)cloneQty/this->config->clonePop);

    double randProb;
    do{
        this->population[iArray]->swapFacility(rand() % this->instance->n, rand() % this->instance->n);
        randProb=(double)(rand()%100)/100.0;
    }while(randProb<probSwap);

    this->population[iArray]->adjustP(); //todo rever
    this->population[iArray]->calculateSolution();

}

//Processo de Re-seleção. Substitui na população original os clones maturados de melhor afinidade
void Search::reselect() {

    int iClone=this->config->pSize;
    int iPop=this->config->pSize-1;

    while(iClone<this->config->arraySize && iPop>=0 && this->population[iClone]->cost<this->population[iPop]->cost){
        swapAntibody(iClone++, iPop--);
    }

    sortPopulation();
}

//Processo de troca dos d piores indivíduos por novos individuos gerados
void Search::regenerate() {

    int iReg=this->config->pSize-1;

    for(int i=0;i<this->config->regenerationQty; i++){
        shuffle(this->population[iReg]->layout,this->population[iReg]->layout+this->instance->n,std::default_random_engine(rand()));
        this->population[iReg]->adjustP();
        this->population[iReg]->calculateSolution();

        iReg--;
    }

    sortPopulation();
}

//Constroi população inicial de anticorpos
void Search::buildInitialPopulation() {

    if(this->config->pSize>0){

        this->population = new Antibody*[this->config->arraySize];

        for(int i=0;i<this->config->pSize;i++){
            buildAntibody(i);
        }

        sortPopulation();

        for(int i=this->config->pSize;i<this->config->arraySize;i++){
            this->population[i]=new Antibody(this->instance);
        }

    }else{
        cout<<"ERROR: config->pSize<=0"<<endl;
        exit(1);
    }
}

//Constroi um anticorpo aleatoriamente com base no outro
void Search::buildAntibody(int index) {

    if(index != FIRST_ANTIBODY){
        this->population[index]=this->population[index-1]->clone();

    }else{
        //Instanciação do primeiro anticorpo
        this->population[index] = new Antibody(this->instance);

        for(int i=0;i<this->instance->n;i++){
            this->population[index]->layout[i]=i;
        }
    }

    shuffle(this->population[index]->layout,this->population[index]->layout+this->instance->n,std::default_random_engine(rand()));
    this->population[index]->adjustP();
    this->population[index]->calculateSolution();

}

//Impressao da população
void Search::printPopulation() {

    cout<<"----POPULACAO---- "<<endl<<this->instance->name<<endl;
    for(int i=0;i<this->config->pSize;i++){
        cout<<"-- "<<i<<" ";
        this->population[i]->print();
        if(i>0 && this->population[i]->cost<this->population[i-1]->cost){
            cout<<"Erro de ordenação"<<endl;
            exit(1);
        }
    }

}

void Search::printAll(){
    cout<<"--------Population:"<<endl;
    for(int i=0;i<this->config->pSize;i++){
        if(i>0 && this->population[i]->cost<this->population[i-1]->cost){
            cout<<"Erro de ordenação aqui-> ";
        }
        cout<<i<<" - "<<this->population[i]->cost<<endl;
    }
    cout<<"--------Clones:"<<endl;
    for(int i=this->config->pSize;i<this->config->arraySize;i++){
        if(i>this->config->pSize && this->population[i]->cost<this->population[i-1]->cost){
            cout<<"Erro de ordenação aqui-> ";
        }
        cout<<i<<" - "<<this->population[i]->cost<<endl;
    }
}
//Impressao dos Clones
void Search::printClones() {

    cout<<"----CLONES---- "<<endl<<this->instance->name<<endl;
    for(int i=this->config->pSize;i<this->config->arraySize; i++){
        this->population[i]->print();
        if(i>this->config->pSize && this->population[i]->cost<this->population[i-1]->cost){
            cout<<"Erro de ordenação"<<endl;
            exit(1);
        }
    }

}

//Critério de comparação de anticorpos (decrescente)
bool Search::antibodyCriterion(Antibody *a, Antibody *b) {
    return a->cost<b->cost;
}

//Ordenação da população
void Search::sortPopulation() {
    stable_sort(this->population, this->population + config->pSize, antibodyCriterion);
}

//Ordenação dos clones
void Search::sortClones() {
    stable_sort(this->population + config->pSize, this->population + config->arraySize, antibodyCriterion);
}

//Faz o swapFacility de anticorpos no vetor de população
void Search::swapAntibody(int i, int j) {
    Antibody* antibody=this->population[i];

    this->population[i]=this->population[j];
    this->population[j]=antibody;
}

void Search::vns(Antibody *antibody) {

    int ls=rand()%3;
    switch (ls) {
        case NBR_SWAP:
            neighborsSwap(antibody);
            break;
        case NONBR_SWAP:
            nonNeiborhsSwap(antibody);
            break;
        case OPPS_SWAP:
            opositeSideSwap(antibody);
            break;
        default:
            cout<<"Algo deu errado"<<endl;
    }

}

void Search::select_local_search() {

    if(this->state == NONE) {
        this->state = rand()%3;
        this->action = this->state;
    }else {

        double rand_eps = (double)(rand()%1001)/1000;

        bool all_actions_blocked = !this->usar_action[0] && !this->usar_action[1] && !this->usar_action[2];

        if(rand_eps < this->config->epsilon || all_actions_blocked) {
            for(int i=0;i<3;i++) {
                this->usar_action[i] = true;
            }
            this->action = rand()%3;
        }else {

            this->espacoAm = 0.0;

            for(int i=0;i<3;i++) {
                if(this->usar_action[i] && this->action!=i) {
                    this->espacoAm += this->q[this->state][i];
                }
            }

            double prob = rand()%((int)ceil(this->espacoAm));
            double x = 0.0;


            int i=-1;

            do {
                if(this->usar_action[i+1] && this->action!=(i+1)) {
                    x+=this->q[this->state][i+1];
                }
                i++;
            }while (x < prob && this->action < 2);
            this->action=i;

        }
    }
}

void Search::vns_q(Antibody *antibody) {

    select_local_search();

    double delta = antibody->cost;

    if(this->action == NBR_SWAP) {
        neighborsSwap(antibody);
    }else if(this->action == NONBR_SWAP) {
        nonNeiborhsSwap(antibody);
    }else if(this->action == OPPS_SWAP) {
        opositeSideSwap(antibody);
    }

    delta = antibody->cost - delta;

    if(this->state != this->action) {
        calculate_reward(delta);
    }

    this->state = this->action;

}
void Search::print_q() {
    cout<<endl;
    for(int i=0;i<3;i++) {
        for(int j=0;j<3;j++) {
            cout<<to_string(this->q[i][j])+",";
        }
        cout<<endl;
    }
}

void Search::calculate_reward(double delta) {

    double best_next_q = -1.0;

    // print_q(); //TODO APAGAR

    for(int i=0;i<3;i++) {
        if(this->q[this->action][i] > best_next_q) {
            best_next_q = this->q[this->action][i];
        }
    }

    double reward;

    if(delta<0) { //Recompensa
        reward = -delta;
        this->last_improvement_delta = -delta;
    }else { //Punição
        reward = -this->last_improvement_delta;
    }

    q[this->state][this->action] = q[this->state][this->action] + this->config->alpha*(reward + this->config->epsilon*best_next_q);
    // print_q();
}

//Fase de Trocas lado a lado //TODO Analisar melhor o movimento
bool Search::neighborsSwap(Antibody *antibody) {

    float bestCost=antibody->cost;

    int iE=1;
    int iD=antibody->p+1;
    while(iE<antibody->p && iD<antibody->instance->n){

        antibody->swapFacility(iE-1,iE);
        antibody->calculateSwap(iE-1,iE);


        if(antibody->cost<bestCost){
            return true;
        }else{
            antibody->swapFacility(iE-1,iE);
            antibody->calculateAbcissa();
            antibody->cost=bestCost;
            iE++;
        }

        antibody->swapFacility(iD-1,iD);
        antibody->calculateSwap(iD-1,iD);

        if(antibody->cost<bestCost){
            return true;
        }else{
            antibody->swapFacility(iD-1,iD);
            antibody->calculateAbcissa();
            antibody->cost=bestCost;
            iD++;
        }

    }

    while(iE<antibody->p){
        antibody->swapFacility(iE-1,iE);
        antibody->calculateSwap(iE-1,iE);

        if(antibody->cost<bestCost){
            return true;
        }else{
            antibody->swapFacility(iE-1,iE);
            antibody->calculateAbcissa();
            antibody->cost=bestCost;
            iE++;
        }
    }

    while(iD<antibody->instance->n){
        antibody->swapFacility(iD-1,iD);
        antibody->calculateSwap(iD-1,iD);

        if(antibody->cost<bestCost){
            return true;
        }else{
            antibody->swapFacility(iD-1,iD);
            antibody->calculateAbcissa();
            antibody->cost=bestCost;
            iD++;
        }
    }
    return false;

}

bool Search::nonNeiborhsSwap(Antibody *antibody) {

    float bestCost=antibody->cost;

    for(int i=0;i<(antibody->p-2);i++){

        for(int j=i+2;j<antibody->p;j++) {
            antibody->swapFacility(i,j);
            antibody->calculateSwap(i,j);

            if (antibody->cost < bestCost) {
                return true;
            } else {
                antibody->swapFacility(i,j);
                antibody->calculateAbcissa();
                antibody->cost=bestCost;
            }
        }
    }

    for(int i=antibody->p;i<(antibody->instance->n-2);i++){

        for(int j=i+2;j<antibody->instance->n;j++) {
            antibody->swapFacility(i,j);
            antibody->calculateSwap(i,j);

            if (antibody->cost < bestCost) {
                return true;
            } else {
                antibody->swapFacility(i,j);
                antibody->calculateAbcissa();
                antibody->cost=bestCost;
            }
        }
    }

    return false;
}
bool Search::opositeSideSwap(Antibody *antibody) {

    float bestCost=antibody->cost;
    for(int i=0;i<antibody->p;i++){
        for(int j=antibody->p;j<antibody->instance->n;j++){
            antibody->swapFacility(i,j);
            antibody->calculateSwap(i,j);

            if (antibody->cost < bestCost) {
                return true;
            } else {
                antibody->swapFacility(i,j);
                antibody->calculateAbcissa();
                antibody->cost=bestCost;
            }
        }
    }
    return false;
}

void Search::testAllPopulation() {
    for(int i=0;i<this->config->pSize;i++){
        this->population[i]->testCalculation();
    }

}




