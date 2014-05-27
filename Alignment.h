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
	void mutate(mt19937& prng, float mutswappb);
	void becomeChild(mt19937& prng, float cxswappb, 
		             const Alignment& p1,
		             const Alignment& p2);
	void computeFitness(const Network& net1,
		                const Network& net2,
		                const BLASTDict d,
		                const vector<string> fitnessNames, 
		                const vector<double> fitnessWeights);
	void save(const Network& net1,
		      const Network& net2,
		      string filename) const;
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

	double ics(const Network& net1, const Network& net2) const;
	double sumBLAST(const Network& net1,
		            const Network& net2,
		            const BLASTDict d) const;
};

vector<vector<Alignment*> > nonDominatedSort(vector<Alignment*> in);

//returns true iff aln1 dominates aln2
bool dominates(Alignment* aln1, Alignment* aln2);