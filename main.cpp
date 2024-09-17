#include <iostream>
#include "Instance.h"
#include "Utils.h"
#include "Search.h"
#include "Config.h"
#include "Stats.h"
#include <vector>

using namespace std;

void testSolution();
int main() {
    //testSolution();exit(0);

    vector<Instance*> instances = Utils::loadInstances("Instances/Instances");
    //500, 5000, 0.2, 0.50, 0.1,30
    Config* config = new Config(10, 500, 0.9, 0.2, 0.1,30);
    Stats* stats= new Stats((int) instances.size(), config, false);
    Stats* stats_q= new Stats((int) instances.size(), config, true);


    for(int i=0;i<instances.size();i++){
        cout<<"Instância "+instances.at(i)->name<<endl;
        cout<<"Execução Normal"<<endl;
        for(int j=0; j < config->executions; j++) {
            srand(config->seeds[j]);
            cout<<"\t"<<j<<" - ";

            Search* search = new Search(config, instances.at(i), stats->litSol[i]);

            clock_t time=clock();
            search->evolve();
            time=clock()-time;

            search->testAllPopulation();

            stats->setStat(j, i, ((double) time / CLOCKS_PER_SEC), search->population[0]->cost);

            cout<<"time: "<<to_string(stats->getTime(j,i))<<"s | cost: "<< to_string(stats->getCost(j,i))<<endl;


            delete search;
        }
        stats->printStats(instances.at(i)->name,i);
        cout<<endl<<"RESUME: AVG TIME: "<<to_string(stats->avgTimes[i])<<"s | BEST COST: "<<to_string(stats->bestCosts[i])<<" | LitSol: "<<to_string(stats->litSol[i])<<" | GAP: "<<to_string(stats->gapsSol[i])<<endl<<endl;

        cout<<"Execução Q-learning"<<endl;
        for(int j=0; j < config->executions; j++) {
            srand(config->seeds[j]);
            cout<<"\t"<<j<<" - ";

            Search* search = new Search(config, instances.at(i), stats_q->litSol[i]);

            clock_t time=clock();
            search->evolve();
            time=clock()-time;

            search->testAllPopulation();

            stats_q->setStat(j, i, ((double) time / CLOCKS_PER_SEC), search->population[0]->cost);

            cout<<"time: "<<to_string(stats->getTime(j,i))<<"s | cost: "<< to_string(stats->getCost(j,i))<<endl;


            delete search;
        }

        stats_q->printStats(instances.at(i)->name,i);
        cout<<endl<<"RESUME: AVG TIME: "<<to_string(stats_q->avgTimes[i])<<"s | BEST COST: "<<to_string(stats_q->bestCosts[i])<<" | LitSol: "<<to_string(stats_q->litSol[i])<<" | GAP: "<<to_string(stats_q->gapsSol[i])<<endl<<endl;

    }


    delete config;
    delete stats;
    delete stats_q;

    return 0;
}


void testSolution(){

    vector<Instance*> instances= Utils::loadInstances("Instances/Instances");

    string instanceName="S9";
    Instance* instance;
    for(Instance* i: instances){
        if(i->name==instanceName){
            instance=i;
            break;
        }
    }

    Antibody* layout=new Antibody(instance);
    int corridor[9]={2,1,4,0,8,7,5,3,6};
    layout->layout=corridor;
    layout->p=5;
    layout->calculateSolution();

    int k=1;
    for(int i=0;i<layout->p-1;i++){
        while(i+k<layout->p-1) {
            for (int j = i + k; j < layout->p; j++) {
                //layout->sameSideCalc(i, j);
                layout->swapFacility(i, j);
                layout->testCalculation();
            }
            k++;
        }
    }
    k=1;
    for(int i=layout->p;i<layout->instance->n-1;i++){
        while(i+k<layout->instance->n-1) {
            for (int j = i + k; j < layout->p; j++) {
                //layout->sameSideCalc(i, j);
                layout->swapFacility(i, j);
                layout->testCalculation();
            }
            k++;
        }
    }

    layout->print();

    exit(0);
}