#pragma once

#include "Alignment.h"
#include "Network.h"
#include <vector>
#include <string>
#include <random>

/*
The functions below, taken together, can be used to implement
all of NSGA-II. This entails:
-doing a non-dominated sort of the population
-for each front, assigning crowding distances
-sorting by the crowded comparison operator
*/

//todo: consider passing vector by reference
//todo: add more pass-by-reference in general
vector<vector<Alignment*> > nonDominatedSort(const vector<Alignment*>& in);

//in is assumed to be a front produced by nonDominatedSort
void setCrowdingDists(vector<Alignment*>& in);

//returns true iff aln1 dominates aln2
bool dominates(Alignment* aln1, Alignment* aln2);

//implements crowded-comparison operator from Deb et al. 2002
bool crowdedComp(Alignment* aln1, Alignment* aln2);

//precondition: in alns have crowdDist, domCount set.
vector<Alignment*> binSel(mt19937& prng,
	                      const vector<Alignment*>& in, unsigned int tournSize);

void reportStats(const vector<Alignment*>& in,
                 const vector<string> fitnessNames, bool verbose);