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

void fastHillClimb(mt19937& prng, Alignment* aln, bool total,
	               int maxIters, const vector<string>& fitnessNames,
	               int objectiveToImprove){
	auto randIndex = uniform_int_distribution<int>(0,aln->aln.size()-1);

	if(fitnessNames.at(objectiveToImprove) == "Size"){

		//make sure not already total
		if(v1Unaligned.size() == 0){
			return;
		}

		for(int i = 0; i < maxIters; i++){
			node x = randIndex(prng);

			//todo: this loop is O(n). Fix!
			while(aln->alnMask[x]){
				x = randIndex(prng);
			}

			vector<double> delta = aln->onBitHypothetical(x);

			bool noWorseAtOthers = true;
			for(int j = 0; j < delta.size(); j++){
				noWorseAtOthers = noWorseAtOthers && delta.at(j) >= 0;
			}

			if(noWorseAtOthers){
				aln->onBit(x);
				aln->computeFitness(fitnessNames);
			}
			
		}
	}
	else{
		for(int i = 0; i < maxIters; i++){
			node x = randIndex(prng);
			node y = x;
			while(y == x){
				y = randIndex(prng);
			}

			
			vector<double> delta = aln->doSwapHypothetical(x,y);
			bool betterAtTgt = false;
			bool noWorseAtOthers = true;
			for(int j = 0; j < delta.size(); j++){
				if(j == objectiveToImprove){
					betterAtTgt = delta.at(j) > 0;
				}
				else{
					noWorseAtOthers = noWorseAtOthers && delta.at(j) >= 0;
				}
			}
			if(betterAtTgt && noWorseAtOthers){
				aln->doSwap(x,y);
				aln->computeFitness(fitnessNames);
			}
			
		}
	}
}