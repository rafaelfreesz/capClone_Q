//
// Created by rafael on 26/07/22.
//

#include "Antibody.h"

Antibody::Antibody(Instance *instance) {
    this->instance=instance;
    this->layout=new int [instance->n];
    this->abcissa=new double[instance->n];
    this->p=0;
    this->cost=0;
    this->improved= false;
}

Antibody::~Antibody() {

    delete []this->layout;
    delete []this->abcissa;
}

void Antibody::swapFacility(int i, int j) {

    int fI=this->layout[i];

    this->layout[i]=this->layout[j];
    this->layout[j]=fI;
}

Antibody *Antibody::clone() {
    Antibody* agClone = new Antibody(this->instance);

    agClone->p=this->p;
    agClone->cost=this->cost;
    agClone->improved=this->improved;

    for(int i=0;i<this->instance->n;i++){
        agClone->abcissa[i]=this->abcissa[i];
        agClone->layout[i]=this->layout[i];
    }

    return agClone;
}

void Antibody::clone(Antibody *targetClone) {

    targetClone->p=this->p;
    targetClone->cost=this->cost;
    targetClone->improved=this->improved;

    for(int i=0;i<this->instance->n;i++){
        targetClone->abcissa[i]=this->abcissa[i];
        targetClone->layout[i]=this->layout[i];
    }
}

//Ajusta o valor de p balanceado TODO tentar verificar pelos dois lados
void Antibody::adjustP() {
    this->p=0;
    double layoutLength=0.0;
    float halfLayoutLenght=this->instance->layoutLengh/2;

    while(layoutLength<halfLayoutLenght){
        layoutLength+=this->instance->lengths[this->layout[++this->p]];
    }
}

void Antibody::calculateAbcissa() {

    /*
    Representing L={0,1,2,3p,4,5,6}
    E:0,1,2
    D:4,5,6
    */

    int lI=0;
    int rI=this->p;
    float lLength=0.0;
    float rLength=0.0;

    while(lI<this->p && rI<this->instance->n){

        int lFacility=this->layout[lI++];
        this->abcissa[lFacility]=lLength+this->instance->halfLengths[lFacility];
        lLength+=this->instance->lengths[lFacility];

        int rFacility=this->layout[rI++];
        this->abcissa[rFacility]=rLength+this->instance->halfLengths[rFacility];
        rLength+=this->instance->lengths[rFacility];

    }

    while(lI<this->p){

        int lFacility=this->layout[lI++];
        this->abcissa[lFacility]=lLength+this->instance->halfLengths[lFacility];
        lLength+=this->instance->lengths[lFacility];

    }

    while(rI<this->instance->n){

        int rFacility=this->layout[rI++];
        this->abcissa[rFacility]=rLength+this->instance->halfLengths[rFacility];
        rLength+=this->instance->lengths[rFacility];

    }

    /* Representing L={0,1,2,3p,4,5,6}
    E:0,1,2
    D:6,5,4

     float tamAcum=0;

    for(int i=0;i<this->p;i++){
        float length=(float)this->instance->lengths[this->layout[i]];
        tamAcum+=length;
        this->abcissa[this->layout[i]]=tamAcum-length/2;
    }

    tamAcum=0;

    for(int i=this->instance->n-1;i>=this->p;i--){
        float length=(float)this->instance->lengths[this->layout[i]];
        tamAcum+=length;
        this->abcissa[this->layout[i]]=tamAcum-length/2;
    }*/



}

void Antibody::calculatCost() {

    this->cost=0.0;

    int limI=this->instance->n-1;
    int limJ=this->instance->n;
    for(int i=0;i<limI;i++){
        for(int j=i+1;j<limJ;j++){
            int fI=this->layout[i];
            int fJ=this->layout[j];
            this->cost+= abs(this->abcissa[fI]-this->abcissa[fJ])*this->instance->demands[fI][fJ];
        }
    }
}

void Antibody::calculateSolution() {
    calculateAbcissa();
    calculatCost();
}

void Antibody::print() {

    cout<<"Layout:"<<endl;
    cout<<"E: ";
    for(int i=0;i<this->p;i++){
        cout<<to_string(this->layout[i])+" ";
    }
    cout<<endl;

    cout<<"D: ";
    for(int i=this->p;i<this->instance->n;i++){
        cout<<to_string(this->layout[i])+" ";
    }
    cout<<endl;

    cout<<"Inline: ";
    for(int i=0;i<this->instance->n;i++){
        cout<<to_string(this->layout[i])+" ";
    }
    cout<<endl;

    cout<<"Abcissa:"<<endl;
    for(int i=0;i<this->instance->n;i++){
        cout<<to_string(this->abcissa[this->layout[i]])+" ";
    }
    cout<<endl;

    cout<<"Cost: "+ to_string(this->cost)<<endl;
}

void Antibody::calculateSwap(int i, int j) {
    int iMin, iMax;
    
    if(i<j){
        iMin=i;
        iMax=j;
    }else{
        iMax=i;
        iMin=j;
    }
    
    if((iMin<this->p && iMax<this->p) || (iMin>=this->p && iMax>=this->p)){
        
        if(iMax-iMin==1){
            neigbohrCalculation(iMin,iMax);
        }else{
            nonNeigbohrCalculation(iMin,iMax);
        }

    }else{//Caso 1
         opositeSideCalculation(iMin, iMax);
    }
}

//Calculo do caso 3 da troca
void Antibody::neigbohrCalculation(int iMin, int iMax) {

    int iMinFac=this->layout[iMin];
    int iMaxFac=this->layout[iMax];

    float iMinOldAbcissa=this->abcissa[iMinFac];
    float iMaxOldAbcissa=this->abcissa[iMaxFac];

    //Redefinição das abcissas
    this->abcissa[iMaxFac]+=this->instance->lengths[iMinFac];
    this->abcissa[iMinFac]-=this->instance->lengths[iMaxFac];

    float deltaDistance;

    for(int i=0;i<this->instance->n;i++){
        if(i != iMin){
            int iFac=this->layout[i];

            deltaDistance= abs(this->abcissa[iFac] - this->abcissa[iMinFac]) - abs(this->abcissa[iFac] - iMinOldAbcissa);
            this->cost+= deltaDistance * this->instance->demands[iFac][iMinFac];

            deltaDistance= abs(this->abcissa[iFac] - this->abcissa[iMaxFac]) - abs(this->abcissa[iFac] - iMaxOldAbcissa);
            this->cost+= deltaDistance * this->instance->demands[iFac][iMaxFac];

        }else{
            i++;
        }
    }
}
//Calculo do caso 2 da troca
void Antibody::nonNeigbohrCalculation(int iMin, int iMax) {
    int iMinFac=this->layout[iMin];
    int iMaxFac=this->layout[iMax];

    float iMaxOldAbcissa=this->abcissa[iMaxFac];
    float iMinOldAbcissa=this->abcissa[iMinFac];

    float iMinCurrenAbcissa;
    float iMaxCurrentAbcissa;

    //Redefinição das abcissas
    if(iMin != 0 && iMin != this->p) {
        int previousFac=this->layout[iMin - 1];
        iMinCurrenAbcissa= this->abcissa[previousFac] + (float)(this->instance->lengths[iMinFac] + this->instance->lengths[previousFac]) / 2;
        this->abcissa[iMinFac] = iMinCurrenAbcissa;
    }else{
        iMinCurrenAbcissa= this->instance->lengths[iMinFac] / 2;
        this->abcissa[iMinFac]=iMinCurrenAbcissa;
    }

    int previousFac=this->layout[iMax - 1];
    iMaxCurrentAbcissa= iMinOldAbcissa + (float)(this->instance->lengths[iMinFac] - this->instance->lengths[iMaxFac]) / 2;
    this->abcissa[iMaxFac] = iMaxCurrentAbcissa;

    //Recalculo

    //Recalculo para intervalo entre indiceMenor e indiceMaior
    float deltaAbcissa=(float)(-this->instance->lengths[iMaxFac] + this->instance->lengths[iMinFac]);

    for(int i= iMin + 1; i < iMax; i++){
        int iFac=this->layout[i];

        float abcissaBackup=this->abcissa[iFac];

        //Inclusao do custo de i com a indiceMenor + atualizacao da abcissa
        this->abcissa[iFac]+=deltaAbcissa;
        float deltaDistance= abs(this->abcissa[iFac] - iMinCurrenAbcissa) - abs(abcissaBackup - iMinOldAbcissa);
        this->cost+= deltaDistance*this->instance->demands[iMinFac][iFac];

        //Inclusao do custo de i com a indiceMaior
        deltaDistance= abs(this->abcissa[iFac] - iMaxCurrentAbcissa) - abs(abcissaBackup - iMaxOldAbcissa);
        this->cost+= deltaDistance*this->instance->demands[iMaxFac][iFac];

        //Recalculo de i com j<indiceMenor
        for(int j=0; j < iMin; j++){
            int jFac=this->layout[j];
            deltaDistance= abs(this->abcissa[iFac] - this->abcissa[jFac]) - abs(abcissaBackup - this->abcissa[jFac]);

            this->cost+=deltaDistance*this->instance->demands[iFac][jFac];

        }
        //Recalculo de i com j>indiceMaior
        for(int j= iMax + 1; j < this->instance->n; j++){
            int jFac=this->layout[j];
            deltaDistance= abs(this->abcissa[iFac] - this->abcissa[jFac]) - abs(abcissaBackup - this->abcissa[jFac]);

            this->cost+=deltaDistance*this->instance->demands[iFac][jFac];

        }
    }

    //Recalculo das salas que estao fora do centro com iFacndiceMenor e iMaxFac
    //Recalculo de i com j<indiceMenor
    for(int j=0; j < iMin; j++){
        int jFac=this->layout[j];
        //Inclusao do custo de i com a indiceMenor
        float deltaDistance= abs(this->abcissa[jFac] - iMinCurrenAbcissa) - abs(this->abcissa[jFac] - iMinOldAbcissa);
        this->cost+= deltaDistance*this->instance->demands[iMinFac][jFac];

        //Inclusao do custo de i com a indiceMaior
        deltaDistance= abs(this->abcissa[jFac] - iMaxCurrentAbcissa) - abs(this->abcissa[jFac] - iMaxOldAbcissa);
        this->cost+= deltaDistance*this->instance->demands[iMaxFac][jFac];
    }

    //Recalculo de i com j>indiceMaior
    for(int j= iMax + 1; j < this->instance->n; j++){
        int jFac=this->layout[j];
        //Inclusao do custo de i com a indiceMenor
        float deltaDistance= abs(this->abcissa[jFac] - iMinCurrenAbcissa) - abs(this->abcissa[jFac] - iMinOldAbcissa);
        this->cost+= deltaDistance*this->instance->demands[iMinFac][jFac];

        //Inclusao do custo de i com a indiceMaior
        deltaDistance= abs(this->abcissa[jFac] - iMaxCurrentAbcissa) - abs(this->abcissa[jFac] - iMaxOldAbcissa);
        this->cost+= deltaDistance*this->instance->demands[iMaxFac][jFac];
    }

}

//Calculo do caso 1 da troca
void Antibody::opositeSideCalculation(int iMin, int iMax){

    int iMinFac=this->layout[iMin];
    int iMaxFac=this->layout[iMax];

    float oldAbcissaL=this->abcissa[iMinFac];
    float oldAbcissaR=this->abcissa[iMaxFac];

    float newAbcissaL= this->abcissa[iMaxFac] + (float)(this->instance->lengths[iMinFac] - this->instance->lengths[iMaxFac]) / 2;
    float newAbcissaR= oldAbcissaL + (float)(this->instance->lengths[iMaxFac] - this->instance->lengths[iMinFac]) / 2;

    //Redefinição das abcissas para iMin e iMax
    this->abcissa[iMinFac]= newAbcissaL;
    this->abcissa[iMaxFac]= newAbcissaR;

    //Recalculo E e D
    float deltaDistance= abs(newAbcissaL - newAbcissaR) - abs(oldAbcissaL - oldAbcissaR);
    this->cost+=deltaDistance*this->instance->demands[iMinFac][iMaxFac];

    float deltaAbcissaE= this->instance->lengths[iMinFac] - this->instance->lengths[iMaxFac];
    float deltaAbcissaD= this->instance->lengths[iMaxFac] - this->instance->lengths[iMinFac];

    //Calculo de i>iMin && i<this->p
    for(int i= iMin + 1; i < this->p; i++){
        int iFac=this->layout[i];
        float abcissaBackup=this->abcissa[iFac];

        this->abcissa[iFac]+=deltaAbcissaE;
        //recalculo com iMinFac
        deltaDistance= abs(this->abcissa[iFac]-this->abcissa[iMinFac]) - abs(abcissaBackup - oldAbcissaL);
        this->cost+=deltaDistance*this->instance->demands[iFac][iMinFac];

        //recalculo com iMaxFac
        deltaDistance= abs(this->abcissa[iFac]-this->abcissa[iMaxFac]) - abs(abcissaBackup - oldAbcissaR);
        this->cost+=deltaDistance*this->instance->demands[iFac][iMaxFac];

        //Recalculo de i com j<iMin
        for(int j=0; j < iMin; j++){
            int jFac=this->layout[j];
            this->cost+=deltaAbcissaE*this->instance->demands[iFac][jFac];
        }

        //Recalculo de i com j<iMax
        for(int j=this->p; j < iMax; j++){
            int jFac=this->layout[j];
            deltaDistance=abs(this->abcissa[jFac]-this->abcissa[iFac])-abs(this->abcissa[jFac]-abcissaBackup);
            this->cost+=deltaDistance*this->instance->demands[iFac][jFac];
        }

    }

    //Calculo de i>iMax && i<this->qtdSalasAlocadas
    for(int i= iMax + 1; i < this->instance->n; i++){
        int iFac=this->layout[i];
        float abcissaBackupI=this->abcissa[iFac];

        this->abcissa[iFac]+=deltaAbcissaD;
        //recalculo com iMinFac
        deltaDistance= abs(this->abcissa[iFac]-this->abcissa[iMinFac]) - abs(abcissaBackupI - oldAbcissaL);
        this->cost+=deltaDistance*this->instance->demands[iFac][iMinFac];

        //recalculo com iMaxFac
        deltaDistance= abs(this->abcissa[iFac]-this->abcissa[iMaxFac]) - abs(abcissaBackupI - oldAbcissaR);
        this->cost+=deltaDistance*this->instance->demands[iFac][iMaxFac];


        //Recalculo de i com j<iMin
        for(int j=0; j < iMin; j++){
            int jFac=this->layout[j];
            deltaDistance=abs(this->abcissa[jFac]-this->abcissa[iFac])-abs(this->abcissa[jFac] - abcissaBackupI);
            this->cost+=deltaDistance*this->instance->demands[iFac][jFac];
        }

        //Recalculo de i com j<iMax
        for(int j=this->p; j < iMax; j++){
            int jFac=this->layout[j];
            this->cost+=deltaAbcissaD*this->instance->demands[iFac][jFac];
        }

        //Recalculo com j>iMin+1 && j<this->p
        for(int j= iMin + 1; j < this->p; j++){
            int jFac=this->layout[j];
            float abcissaBackupJ=this->abcissa[jFac]-deltaAbcissaE;

            deltaDistance=abs(this->abcissa[jFac]-this->abcissa[iFac])-abs(abcissaBackupJ - abcissaBackupI);
            this->cost+=deltaDistance*this->instance->demands[iFac][jFac];
        }

    }

    //Calculo de i<iMin com iMin e iMax
    for(int i=0; i < iMin; i++){
        int iFac=this->layout[i];

        //recalculo com iMinFac
        deltaDistance= abs(this->abcissa[iFac]-this->abcissa[iMinFac]) - abs(this->abcissa[iFac] - oldAbcissaL);
        this->cost+=deltaDistance*this->instance->demands[iFac][iMinFac];

        //recalculo com iMaxFac
        deltaDistance= abs(this->abcissa[iFac]-this->abcissa[iMaxFac]) - abs(this->abcissa[iFac] - oldAbcissaR);
        this->cost+=deltaDistance*this->instance->demands[iFac][iMaxFac];

    }

    //Calculo de i>=this->particoa && i<iMax com iMin e iMax
    for(int i=this->p; i < iMax; i++){
        int iFac=this->layout[i];

        //recalculo com iMinFac
        deltaDistance= abs(this->abcissa[iFac]-this->abcissa[iMinFac]) - abs(this->abcissa[iFac] - oldAbcissaL);
        this->cost+=deltaDistance*this->instance->demands[iFac][iMinFac];

        //recalculo com iMaxFac
        deltaDistance= abs(this->abcissa[iFac]-this->abcissa[iMaxFac]) - abs(this->abcissa[iFac] - oldAbcissaR);
        this->cost+=deltaDistance*this->instance->demands[iFac][iMaxFac];

    }

}

void Antibody::testCalculation() {

    Antibody* clone = this->clone();

    calculateSolution();

    for(int i=0;i<this->instance->n;i++){
        if(this->abcissa[this->layout[i]]!=clone->abcissa[this->layout[i]]){
            cout<<"Abcissa differs:\nCorrect:"<<endl;
            print();
            cout<<"Given:"<<endl;
            clone->print();
            exit(2);
        }
    }

    if(this->cost!=clone->cost){
        cout<<"Cost differs!"<<endl<<this->cost<<endl<<clone->cost<<endl;
    }

    delete clone;

}

void Antibody::sameSideCalc(int iMin, int iMax) {
    int iMinFacility=this->layout[iMin];
    int iMaxFacility=this->layout[iMax];

    double iMinNewAbcissa=this->abcissa[iMaxFacility]-this->instance->halfLengths[iMinFacility]+this->instance->halfLengths[iMaxFacility];
    double iMaxNewAbcissa=this->abcissa[iMinFacility]-this->instance->halfLengths[iMinFacility]+this->instance->halfLengths[iMaxFacility];

    int zoneLimits[8]; //{zone1,zone2,zone3,zone4} *Zone 2 may not exist
    //Setting Zone Limits
    if(iMinFacility<this->p){
        zoneLimits[0]=0;
        zoneLimits[1]=iMin;           //Zone 1 =[0,..,iMin)
        zoneLimits[2]=iMin+1;
        zoneLimits[3]=iMax;             //Zone 2 =[iMin+1,...,iMax)
        zoneLimits[4]=iMax+1;
        zoneLimits[5]=this->p;          //Zone 3= [iMax+1,...,p)
        zoneLimits[6]=this->p;
        zoneLimits[7]=this->instance->n;   //Zone 4= [p,...,assigned)
    }else{
        zoneLimits[0]=this->p;
        zoneLimits[1]=iMin;           //Zone 1 =[p,..,iMin)
        zoneLimits[2]=iMin+1;
        zoneLimits[3]=iMax;             //Zone 2 =[iMin+1,...,iMax)
        zoneLimits[4]=iMax+1;
        zoneLimits[5]=this->instance->n;   //Zone 3= [iMax+1,...,assigned)
        zoneLimits[6]=0;
        zoneLimits[7]=this->p;          //Zone 3= [0,...,p)
    }

    float deltaDist, deltaCost;
    int iFacility, jFacility;
    //Delta calculating by iMin & iMax with zones 1, 3, and 4
    int z=0;
    while(z<8){
        for(int i=zoneLimits[z];i<zoneLimits[z+1];i++){
            iFacility=this->layout[i];

            //iMin x zone z
            deltaDist=abs(iMinNewAbcissa-this->abcissa[iFacility])-abs(this->abcissa[iMinFacility]-this->abcissa[iFacility]);
            deltaCost=deltaDist*this->instance->demands[iMinFacility][iFacility];
            this->cost+=deltaCost;

            //iMax x zone z
            deltaDist=abs(iMaxNewAbcissa-this->abcissa[iFacility])-abs(this->abcissa[iMaxFacility]-this->abcissa[iFacility]);
            deltaCost=deltaDist*this->instance->demands[iMaxFacility][iFacility];
            this->cost+=deltaCost;
        }

        if(z==0){
            z+=4;
        }else{
            z+=2;
        }
    }

    //Delta calculating by zone 2 with iMin & iMax *Zone 2 may not exist
    double zone2Delta=(float)(this->instance->lengths[iMaxFacility]-this->instance->lengths[iMinFacility]);
    for(int i=zoneLimits[2];i<zoneLimits[3];i++){
        iFacility=this->layout[i];

        //zone 2 -> iMin
        deltaDist=abs(iMinNewAbcissa-(this->abcissa[iFacility]+zone2Delta))-abs(this->abcissa[iMinFacility]-this->abcissa[iFacility]);
        deltaCost=deltaDist*this->instance->demands[iMinFacility][iFacility]; //TODO maybe detacost+= and this->cost+= at the end
        this->cost+=deltaCost;

        //zone 2 -> iMax
        deltaDist=abs(iMaxNewAbcissa-(this->abcissa[iFacility]+zone2Delta))-abs(this->abcissa[iMaxFacility]-this->abcissa[iFacility]);
        deltaCost=deltaDist*this->instance->demands[iMaxFacility][iFacility];
        this->cost+=deltaCost;

        //Delta calculating by zone 2 with zones 1, 3 and 4 TODO needed only if l(iMinFacility)!= l(iMaxFacility)
        if(zone2Delta!=0.0) {
            z = 0;
            while (z < 8) {

                for (int j = zoneLimits[z]; j < zoneLimits[z + 1]; j++) {
                    jFacility = this->layout[j];
                    deltaDist = abs(this->abcissa[jFacility] - (this->abcissa[iFacility] + zone2Delta)) -
                                abs(this->abcissa[jFacility] - this->abcissa[iFacility]);
                    deltaCost = deltaDist * this->instance->demands[iFacility][jFacility];
                    this->cost += deltaCost;
                }

                if (z == 0) {
                    z += 4;
                } else {
                    z += 2;
                }
            }
        }
    }

    //UpdatingAbcissa //TODO PAREI AQUI
    this->abcissa[iMaxFacility]=iMaxNewAbcissa;
    this->abcissa[iMinFacility]=iMinNewAbcissa;

    float deltaAbcissa = this->instance->lengths[iMaxFacility]-this->instance->lengths[iMinFacility];
    for(int i=iMin+1;i<iMax;i++){
        this->abcissa[this->layout[i]]+=deltaAbcissa;
    }
}













