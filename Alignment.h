#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include <unordered_set>

#include "blastinfo.h"
#include "Network.h"

class Alignment{
public:	
	Alignment(const Network* n1, const Network* n2, 
		      const BLASTDict* bit);
	Alignment(const Network* n1, const Network* n2, string filename,
		      const BLASTDict* bit);
	Alignment(mt19937& prng, float cxswappb, 
		             const Alignment& p1,
		             const Alignment& p2,
		             bool total = true);
	void shuf(mt19937& prng, bool total =true); //shuffles the alignment to make it completely random
	void mutate(mt19937& prng, float mutswappb, bool total = true);
	void doSwap(node x, node y);
	void computeFitness(const BLASTDict& bitscores,
		                const BLASTDict& evalues,
		                const vector<string>& fitnessNames);
	void save(string filename) const;
	double ics() const;
	
	double sumBLAST() const;
	double alnSize() const;	
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
	int actualSize; //number of nodes in net1
	int domRank; //which front this aln is in. 0 is best, 1 is 2nd best, etc.
	int numThatDominate; //how many others in the population dominate this one.
	unordered_set<Alignment*> dominated;//set of alns this aln dominates.
	double crowdDist;

	//new fitness system
	double currBitscore;
	const BLASTDict* bitscores;
	void updateBitscore(node n1, node n2old, node n2new, bool oldMask,
		                bool newMask);

	//stored info version of ICS for fast computation
	//(incrementally update as the alignment changes)
	double fastICSDenominator() const;
	double fastICSNumerator() const;
	double fastICS() const;
	const Network* net1;
	const Network* net2;
	//conservedCounts counts the number of conserved edges that node i
	//is incident to. Summing it and dividing by 2 (to fix double-counting)
	//gives the number of conserved edges in the alignment.
	vector<int> conservedCounts;
	void updateConservedCount(node n1, node n2old, node n2new, bool oldMask,
		                     bool newMask, node ignore);
	void initConservedCount(node n1, node n2, bool mask);
};

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

void reportStats(const vector<Alignment*>& in, bool verbose);