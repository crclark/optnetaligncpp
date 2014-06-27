#include "localSearch.h"

#include "Alignment.h"
#include "Network.h"
#include "nsga-ii.h"
#include <random>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <thread>
using namespace std;

//instead of using the buggy hypotheticals, do an actual swap and
//check for improvement, undoing if worse.
//when obj is specified, only checks for improvement in objective obj
void correctHillClimb(RandGenT& prng, Alignment* aln, bool total,
					  int maxIters, const vector<string>& fitnessNames,
					  int obj){
	auto randNonDummyIndex = uniform_int_distribution<int>(0, aln->actualSize-1);
	auto randIndex = uniform_int_distribution<int>(0,aln->aln.size()-2);
	for(int i = 0; i < maxIters; i++){
		node x = randNonDummyIndex(prng);
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

//optimizes given objective with time proportional to given proportion,
//and evenly distributes the rest of hillclimbing time to the other
//objectives.
void proportionalSearch(RandGenT& prng, Alignment* aln, bool total,
	                    int iters, const vector<string>& fitnessNames,
	                    int obj, double proportion){

	auto prob = uniform_real_distribution<double>(0,1);

	auto randObj = uniform_int_distribution<int>(0,fitnessNames.size()-1);

	for(int i = 0; i < iters; i++){
		double res = prob(prng);
		if(res < proportion){
			correctHillClimb(prng, aln, total,
	               500, fitnessNames,obj);
		}
		else{
			int robj = obj;
			while(robj == obj){
				robj = randObj(prng);
			}
			correctHillClimb(prng, aln, total,
	               500, fitnessNames,robj);	
		}
	}

}

//returns -1.0 if any objectives are worsened. 
//Otherwise, returns magnitude of pct improvement vector
double swapNormalizedDelta(Alignment& aln, const vector<string>& 
						   fitnessNames, node x, node y){
	vector<double> currFit = aln.fitness;

	aln.doSwap(x,y);
	aln.computeFitness(fitnessNames);
	vector<double> newFit = aln.fitness;
	aln.doSwap(x,y);
	aln.computeFitness(fitnessNames);

	vector<double> pctDeltas;

	bool oneNeg = false;
	for(int i = 0; i < newFit.size(); i++){
		double pctDelt = (newFit[i] - currFit[i])/currFit[i];
		pctDeltas.push_back(pctDelt);
		if(pctDelt < 0){
			oneNeg = true;
		}
	}

	if(oneNeg){
		return -1.0;
	}
	else{
		double sumSq = 0.0;
		for(int i = 0; i < pctDeltas.size(); i++){
			sumSq += pctDeltas[i]*pctDeltas[i];
		}

		return sqrt(sumSq);
	}
}

void steepestAscentHillClimb(Alignment* aln, 
							 vector<string>& fitnessNames,
							 int nthreads, bool verbose){

	vector<node> bestXs(nthreads,-1);
	vector<node> bestYs(nthreads,-1);
	vector<double> bestDeltas(nthreads,-1.0);

	auto worker = [&](int tid, int xmin, int xmax, int ymin, int ymax){
		Alignment localAln(*aln);
		node bestX = -1;
		node bestY = -1;
		double bestDelta = -1.0;

		for(node x = xmin; x < xmax; x++){
			for(node y = ymin; y < ymax; y++){
				double newDelta = 
					swapNormalizedDelta(localAln, fitnessNames,
										x,y);
				if(newDelta > bestDelta){
					bestDelta = newDelta;
					bestX = x;
					bestY = y;
				}
			}
		}

		bestXs[tid] = bestX;
		bestYs[tid] = bestY;
		bestDeltas[tid] = bestDelta;
	};

	bool done = false;
	int numiters = 0;
	while(!done){
		//launch worker threads
		int grainsize = aln->actualSize / nthreads;
		vector<thread> ts;
		for(int tid = 0; tid < nthreads; tid++){
			int tminx = tid*grainsize;
			int tmaxx;
			if(tid == nthreads-1){
				tmaxx = aln->actualSize;
			}
			else{
				tmaxx = tminx + grainsize;
			}
			ts.push_back(thread(worker,tid,tminx,tmaxx,0,aln->aln.size()));
		}

		//join threads
		for(int i = 0; i < nthreads; i++){
			ts.at(i).join();
		}

		//find absolute best swap and commit to it
		node bestX = -1;
		node bestY = -1;
		double bestDelt = -1.0;

		for(int i = 0; i < bestDeltas.size(); i++){
			if(bestDeltas.at(i) > bestDelt){
				bestDelt = bestDeltas.at(i);
				bestX = bestXs.at(i);
				bestY = bestYs.at(i);
			}
		}

		if(bestDelt <= 0.0){
			done = true;
		}
		else{
			aln->doSwap(bestX, bestY);
			aln->computeFitness(fitnessNames);
			if(verbose){
				reportStats({aln}, fitnessNames, true);
				cout<<(++numiters)<<" swaps performed."<<endl;
			}
		}

	}
}

VelocityTracker::VelocityTracker(){
	nextSpot = 0;
	size = 0;

	recentDeltas = vector<vector<double> >(500);

}

void VelocityTracker::reportDelta(const vector<double>& in){
	recentDeltas[nextSpot] = vector<double>(in);
	nextSpot = (nextSpot + 1) % 500;

	if(size < 500){
		size++;
	}
}

vector<double> VelocityTracker::getRecentVel() const{
	vector<double> toReturn(recentDeltas[0].size(), 0.0);
	
	for(int i = 0; i < size; i++){
		for(int j = 0; j < recentDeltas[i].size(); j++){
			toReturn[j] += recentDeltas[i][j];
		}
	}
	
	double n = double(size);
	for(int j = 0; j < toReturn.size(); j++){
		toReturn[j] /= n;
	}

	return toReturn;
}