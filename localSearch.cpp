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
Alignment* hillClimb(RandGenT& prng, Alignment* orig, bool total,
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

void fastHillClimb(RandGenT& prng, Alignment* aln, bool total,
	               int maxIters, const vector<string>& fitnessNames,
	               int objectiveToImprove, bool worsenOthers){
	auto randIndex = uniform_int_distribution<int>(0,aln->aln.size()-1);

	if(fitnessNames.at(objectiveToImprove) == "Size"){

		for(int i = 0; i < maxIters; i++){
			//check alignment not total
			if(aln->v1Unaligned.size() == 0){
			return;
			}

			auto randCountdown = uniform_int_distribution<int>(0,aln->v1Unaligned.size());

			int counter = randCountdown(prng);

			auto iter = aln->v1Unaligned.begin();
			
			node x;

			while(counter>0){
				x = *iter;
				counter--;
			}

			vector<double> delta = aln->onBitHypothetical(x);

			bool noWorseAtOthers = true;

			if(!worsenOthers){
				for(int j = 0; j < delta.size(); j++){
					noWorseAtOthers = noWorseAtOthers && delta.at(j) >= 0;
				}
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
				else if(!worsenOthers){
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


//instead of using the buggy hypotheticals, do an actual swap and
//check for improvement, undoing if worse.
void correctHillClimb(RandGenT& prng, Alignment* aln, bool total,
					  int maxIters, const vector<string>& fitnessNames){
	
	auto randIndex = uniform_int_distribution<int>(0,aln->aln.size()-1);

	for(int i = 0; i < maxIters; i++){
		node x = randIndex(prng);
		node y = x;
		while(y == x){
			y = randIndex(prng);
		}

		vector<double> currFit = aln->fitness;

		aln->doSwap(x,y);
		aln->computeFitness(fitnessNames);

		vector<double> newFit = aln->fitness;

		bool improved = true;

		for(int j = 0; j < currFit.size(); j++){
			if(newFit[j] < currFit[j]){
				improved = false;
			}
		}

		//if not an improvement, undo the swap
		if(!improved){
			aln->doSwap(x,y);
			aln->computeFitness(fitnessNames);
		} 
	}
}

//optimizes objectives 0 and 1, with proportion to 0 determined by last arg
void proportionalSearch(RandGenT& prng, Alignment* aln, bool total,
	                    int iters, const vector<string>& fitnessNames,
	                    double proportion){

	auto prob = uniform_real_distribution<double>(0,1);

	for(int i = 0; i < iters; i++){
		double res = prob(prng);
		if(res < proportion){
			correctHillClimb(prng, aln, total,
	               500, fitnessNames, 0, true);
		}
		else{
			correctHillClimb(prng, aln, total,
	               500, fitnessNames, 1, true);	
		}
	}

}
