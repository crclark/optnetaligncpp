#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <string>

#include "blastinfo.h"
#include "Network.h"

class Alignment{
public:	
	Alignment(const Network& net1, const Network& net2);
	Alignment(const Network& net1, const Network& net2, string filename);
	void shuf(mt19937& prng); //shuffles the alignment to make it completely random
	void mutate(mt19937& prng, float mutswappb);
	void becomeChild(mt19937& prng, float cxswappb, 
		             const Alignment& p1,
		             const Alignment& p2);
	void computeFitness(const Network& net1,
		                const Network& net2,
		                const BLASTDict& bitscores,
		                const BLASTDict& evalues,
		                const vector<string>& fitnessNames);
	void save(const Network& net1,
		      const Network& net2,
		      string filename) const;
	double ics(const Network& net1, const Network& net2) const;
	double sumBLAST(const Network& net1,
		            const Network& net2,
		            const BLASTDict& d) const;	
	vector<node> aln;
	vector<bool> alnMask; //indicates whether the corresponding node
	                      //of V2 is allowed to be aligned to.
	                      //The idea here is that a node in the larger
	                      //net may have no ortholog in the smaller,
	                      //so we might want to give up on aligning it.
						  //This presents some difficulties in that the
		                  //same alignment has many representations.
						  //A node in V2 can be unaligned because it is
	                      //aligned to a dummy node or because its mask
	                      //bit is switched off.
		                  // but if we separate searching through permutations
	                      //from searching through masks, it should be 
	                      //workable.
	bool fitnessValid;
	vector<double> fitness; //all fitnesses stored s.t. larger is better.
	int domCount;
	double crowdDist;
};

/*
The functions below, taken together, can be used to implement
all of NSGA-II. This entails:
-doing a non-dominated sort of the population
-for each front, assigning crowding distances
-sorting by the crowded comparison operator
-selecting the individuals to populate the next generation (todo)
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

void reportStats(const vector<Alignment*>& in);