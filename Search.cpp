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

}

Search::~Search() {

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
            rvnd(this->population[i]);
            this->population[i]->improved= true;
        }
    }
    sortPopulation();
}

void Search::improveMemory_q() {
    for(int i=0;i<this->config->memorySetSize;i++){
        if(!this->population[i]->improved) {
            rvnd_q(this->population[i]);
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
            //this->population[iArray]=this->population[i]->clone(); //Clonagem do anticorpo
            this->population[iArray]->improved= false;//Clonagem do anticorpo
            maturate(iArray++, this->config->clonesPerI[i]); //Maturação (mutação) do indivíduo

        }

    }

    sortClones();


    rvnd(this->population[this->config->pSize]);
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


    rvnd_q(this->population[this->config->pSize]);
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

void Search::localSearch(Antibody *antibody) {

    double bestCost=antibody->cost;

    for(int i=0;i<this->instance->n-1;i++){

        for(int j=i+1;j<this->instance->n;j++){
            antibody->swapFacility(i, j);
            antibody->calculateSwap(i,j);

            if(antibody->cost<bestCost){
                bestCost= antibody->cost;
                i=-1;
                break;
            }else{
                antibody->swapFacility(i, j);
                antibody->calculateSwap(i,j);
            }

        }
    }
}

void Search::rvnd(Antibody *antibody) {

    int searchSequence[3]={0, 1, 2};
    int lastImproved;
    shuffle(searchSequence,searchSequence+3,std::default_random_engine(rand()));

   // 1 0 2
    int fase=0;
    while(fase < 3){

        switch (searchSequence[fase]) {
            case 0:
                if(searchSequence[fase]!=lastImproved) {
                    if (neighborsSwap(antibody)) {
                        lastImproved = searchSequence[fase];
                        fase = 0;
                    } else {
                        fase++;
                    }
                }else{
                    fase++;
                }
                break;
            case 1:
                if(searchSequence[fase]!=lastImproved) {
                    if (nonNeiborhsSwap(antibody)) {
                        lastImproved = searchSequence[fase];
                        fase = 0;
                    } else {
                        fase++;
                    }
                }else{
                    fase++;
                }
                break;
            case 2:
                if(searchSequence[fase]!=lastImproved) {
                    if (opositeSideSwap(antibody)) {
                        lastImproved = searchSequence[fase];
                        fase = 0;
                    } else {
                        fase++;
                    }
                }else{
                    fase++;
                }
                break;
        }

    }

}
void Search::rvnd_q(Antibody *antibody) {


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

void Search::PathR(Antibody *antibody) {

   /* Antibody* bestAntibody=antibody->clone();

    float foCorrente=antibody->cost;

    bool* trocada=new bool[antibody->instance->n];
    for(int i=0;i<antibody->instance->n;i++) trocada[i]=false;

    int iE=0;
    int iD= antibody->p - 1;
    int melhortroca[2]={-1,-1};
    float melhorParcial=-1;

    int nTrocas=antibody->instance->n;
    if(antibody->p%2==1){
        nTrocas--;
    }
    if((antibody->instance->n-antibody->p)%2==1){
        nTrocas--;
    }

    int trocas=0;
    while(trocas<nTrocas) {
        while (iE < iD) {
            if(!trocada[iE]) {

                antibody->swapFacility(iE,iD);
                antibody->calculateSwap(iE, iD);

                if (antibody->cost < melhorParcial || melhorParcial == -1) {
                    melhortroca[0] = iE;
                    melhortroca[1] = iD;
                    melhorParcial = antibody->cost;
                }

                antibody->swapFacility(iE,iD);
                antibody->calculateAbcissa();
                antibody->cost=foCorrente;
            }else{
                iE++;
                iD--;
            }
        }
        iE=antibody->p;
        iD=antibody->instance->n-1;

        while (iE < iD) {
            if(!trocada[iE]) {
                antibody->swapFacility(iE,iD);
                antibody->calculateSwap(iE, iD);

                if (antibody->cost < melhorParcial || melhorParcial == -1) {
                    melhortroca[0] = iE;
                    melhortroca[1] = iD;
                    melhorParcial = antibody->cost;
                }

                antibody->swapFacility(iE,iD);
                antibody->calculateAbcissa();
                antibody->cost=foCorrente;
            }else{
                iE++;
                iD--;
            }
        }

        trocas+=2;
        antibody->swapFacility(melhortroca[0], melhortroca[1]);
        antibody->calculateSwap(melhortroca[0], melhortroca[1]);
        trocada[melhortroca[0]]=true;
        trocada[ melhortroca[1]]=true;
        foCorrente=antibody->cost;

        melhortroca[0] = -1;
        melhortroca[1] = -1;
        melhorParcial = -1;

        iE = 0;
        iD = antibody->p - 1;
        if (antibody->cost < bestAntibody->cost) {
            bestAntibody->armazenar(this->corredorSolucao, antibody->p, antibody->cost);
        }


    }
    melhorSolucao->restaurar(this->corredorSolucao, &antibody->p, &antibody->cost);
    montarAbcissas();

    delete melhorSolucao;
    delete[] trocada;*/
}

void Search::testAllPopulation() {
    for(int i=0;i<this->config->pSize;i++){
        this->population[i]->testCalculation();
    }

}




