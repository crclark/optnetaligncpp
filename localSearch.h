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


class VelocityTracker{
public:
	VelocityTracker();
	void reportDelta(const vector<double>& in);
	vector<double> getRecentVel() const;

private:
	int nextSpot;
	int size;
	vector<vector<double> > recentDeltas;
};