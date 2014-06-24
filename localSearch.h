#pragma once

#include "Alignment.h"
#include "Network.h"

#include <random>
#include <vector>

void correctHillClimb(RandGenT& prng, Alignment* aln, bool total,
					  int maxIters, const vector<string>& fitnessNames);

void proportionalSearch(RandGenT& prng, Alignment* aln, bool total,
	                    int iters, const vector<string>& fitnessNames,
	                    double proportion);