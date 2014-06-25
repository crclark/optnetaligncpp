#include "localSearch.h"

#include "Alignment.h"
#include "Network.h"
#include "nsga-ii.h"
#include <random>
#include <iostream>
#include <vector>
#include <string>
using namespace std;

//instead of using the buggy hypotheticals, do an actual swap and
//check for improvement, undoing if worse.
//when obj is specified, only checks for improvement in objective obj
void correctHillClimb(RandGenT& prng, Alignment* aln, bool total,
					  int maxIters, const vector<string>& fitnessNames,
					  int obj){
	
	auto randIndex = uniform_int_distribution<int>(0,aln->aln.size()-2);

	for(int i = 0; i < maxIters; i++){
		node x = randIndex(prng);
		node y = randIndex(prng);
		if(y >= x){
			y++;
		}

		vector<double> currFit = aln->fitness;

		aln->doSwap(x,y);
		aln->computeFitness(fitnessNames);

		vector<double> newFit = aln->fitness;

		bool improved = true;

		if(obj == -1){
			for(int j = 0; j < currFit.size(); j++){
				if(newFit[j] < currFit[j]){
					improved = false;
				}
			}
		}
		else{
			if(newFit[obj] <= currFit[obj]){
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
	               500, fitnessNames,0);
		}
		else{
			correctHillClimb(prng, aln, total,
	               500, fitnessNames,1);	
		}
	}

}
