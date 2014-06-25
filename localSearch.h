#pragma once

#include "Alignment.h"
#include "Network.h"

#include <random>
#include <vector>

void correctHillClimb(RandGenT& prng, Alignment* aln, bool total,
					  int maxIters, const vector<string>& fitnessNames,
					  int obj = -1);

void proportionalSearch(RandGenT& prng, Alignment* aln, bool total,
	                    int iters, const vector<string>& fitnessNames,
	                    int obj, double proportion);

//todo: we could do this in place with one int.
//Might be worth it if this is slow.
class VelocityTracker{
public:
	VelocityTracker();
	void reportDelta(vector<double>& in);
	vector<double> getRecentVel() const;

private:
	int nextSpot;
	int size;
	vector<vector<double> > recentDeltas;
};