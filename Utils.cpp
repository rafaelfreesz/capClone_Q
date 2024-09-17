//
// Created by rafael on 22/07/22.
//

#include <fstream>
#include "Utils.h"
#include "Antibody.h"

vector<Instance *> Utils::loadInstances(string strFileName) {
    vector<Instance*> instances;
    string line;
    ifstream fileInstName;

    fileInstName.open(strFileName.c_str());

    //Tratando cada instância
    while(getline(fileInstName,strFileName)){
        line="Instances/"+strFileName;
        ifstream fileInst;
        fileInst.open(line.c_str());

        //Nome e n
        getline(fileInst,line);
        Instance* instance = new Instance(stoi(line));
        instance->name=strFileName;

        //Comprimentos
        getline(fileInst,line);
        vector<string> tokenizedLengths;
        tokenize(line,tokenizedLengths,",");
        if(tokenizedLengths.size()!=instance->n){
            cout<<"Diferença de tamanho, instancia " + instance->name<<endl;
            exit(2);
        }

        for(int i=0;i<instance->n;i++){
            if(getline(fileInst,line)){
                vector<string> tokenizedDemands;
                tokenize(line,tokenizedDemands,",");

                if(tokenizedDemands.size()==instance->n){
                    for(int j=0;j<instance->n;j++){
                        instance->demands[i][j]= stoi(tokenizedDemands.at(j));
                    }
                }else{
                    cout<<"Diferença de tamanho, instancia " + instance->name<<endl;
                    exit(2);
                }
                }else{
                    cout<<"Linhas não batem"<<endl;
                }
            }

        for(int i=0;i<instance->n;i++){
            instance->lengths[i]= stoi(tokenizedLengths.at(i));
            instance->halfLengths[i]=(double)instance->lengths[i]/2.0;
        }

        instance->calculateLayoutLength();
        instance->verify();
        fileInst.close();
        instances.push_back(instance);
    }

    fileInstName.close();

    return instances;
}

void Utils::tokenize(string str, vector<string> &token_v, string DELIMITER) {
    size_t start = str.find_first_not_of(DELIMITER), end=start;

    while (start != std::string::npos) {
        end = str.find(DELIMITER, start);
        token_v.push_back(str.substr(start, end-start));
        start = str.find_first_not_of(DELIMITER, end);
    }
}

void Utils::testInstance(vector<Instance *> instances) {
    int i=0;
    while(instances.at(i)->name!="S10"){i++;}

    Instance* instance=instances.at(i);
    instance->print();

    Antibody* antibody = new Antibody(instance);

    string layoutConf="5,1,3,4,6,8,2,0,9,7,6";

    vector<string> tokens;
    Utils::tokenize(layoutConf,tokens,",");

    for(i=0;i<instance->n;i++){
        antibody->layout[i]= stoi(tokens.at(i));
    }
    antibody->p=stoi(tokens.at(instance->n));

    antibody->calculateSolution();
    antibody->print();
    delete antibody;
}

void Utils::swapInt(int i, int j, int *array) {
    int aux=array[i];

    array[i]=array[j];
    array[j]=aux;
}
