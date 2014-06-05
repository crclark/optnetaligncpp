#include "localSearch.h"

#include "Alignment.h"
#include "Network.h"
#include "nsga-ii.h"
#include <random>
#include <iostream>
#include <vector>
#include <string>
using namespace std;

//local search improving one specified objective
Alignment* hillClimb(mt19937& prng, Alignment* orig, bool total,
	                 int maxIters, const vector<string>& fitnessNames,
	                 int objectiveToImprove){
	int size = orig->actualSize;
	double mutpb = 1.0/(double(size));
	//for stats tracking, get start objective value
	double startObj = orig->fitness.at(objectiveToImprove);
	int numMuts = 0;
	int numImprovingMuts = 0;
	Alignment* curr = orig;
	for(int i = 0; i < maxIters; i++){
		//generate up to size alignments each iteration.
		//change to first that is better
		for(int j = 0; j < 100; j++){
			Alignment* temp = new Alignment(*curr);
			temp->mutate(prng, mutpb, total);
			temp->computeFitness(fitnessNames);
			numMuts++;
			if(temp->fitness.at(objectiveToImprove) > 
			   curr->fitness.at(objectiveToImprove)){
				if(curr != orig){
					delete curr;
				}
				numImprovingMuts++;
				curr = temp;
				continue;
			}
			else{
				delete temp;
			}
			//cout<<"evaluated "<<j<<" mutations"<<endl;
		}
	}

	double endObj = curr->fitness.at(objectiveToImprove);
	cout<<"hillClimb improvement to "<<fitnessNames.at(objectiveToImprove)
	    <<" is "<<(endObj - startObj)<<endl;
	cout<<"Mutation improved "<<fitnessNames.at(objectiveToImprove)
	    <<" in "<<numImprovingMuts<<"/"<<numMuts<<" = "
	    <<(double(numImprovingMuts)/double(numMuts))<<" cases."<<endl;
	return curr;
}