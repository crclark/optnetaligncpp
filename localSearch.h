#pragma once

#include "Alignment.h"
#include "Network.h"

#include <random>
#include <vector>

Alignment* hillClimb(mt19937& prng, Alignment* orig, bool total,
	                 int maxIters, const vector<string>& fitnessNames,
	                 int objectiveToImprove);

void fastHillClimb(mt19937& prng, Alignment* aln, bool total,
	               int maxIters, const vector<string>& fitnessNames,
	               int objectiveToImprove);