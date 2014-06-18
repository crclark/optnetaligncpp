#pragma once

#include "Alignment.h"
#include "Network.h"

#include <random>
#include <vector>

Alignment* hillClimb(RandGenT& prng, Alignment* orig, bool total,
	                 int maxIters, const vector<string>& fitnessNames,
	                 int objectiveToImprove);

void fastHillClimb(RandGenT& prng, Alignment* aln, bool total,
	               int maxIters, const vector<string>& fitnessNames,
	               int objectiveToImprove, bool worsenOthers);

void proportionalSearch(RandGenT& prng, Alignment* aln, bool total,
	                    int iters, const vector<string>& fitnessNames,
	                    double proportion);