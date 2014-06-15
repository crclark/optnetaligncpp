#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include <unordered_set>

#include "blastinfo.h"
#include "Network.h"
#include "goc.h"

class Alignment{
public:	
    //simplest constructor creates an arbitrary non-random alignment
	Alignment(const Network* n1, const Network* n2, 
		      const BLASTDict* bit, const GOCDict* goc);
    //load alignment from file
	Alignment(const Network* n1, const Network* n2, string filename,
		      const BLASTDict* bit, const GOCDict* goc);
	//crossover constructor (UPMX)
    Alignment(mt19937& prng, float cxswappb, 
		             const Alignment& p1,
		             const Alignment& p2,
		             bool total = true);
	void greedyBitscoreMatch();
	void shuf(mt19937& prng, bool uniformsize,
	          bool smallstart, bool total); //shuffles the alignment to make it completely random
	void mutate(mt19937& prng, float mutswappb, bool total = true);
	void doSwap(node x, node y);
	vector<double> doSwapHypothetical(node x, node y) const;
	void onBit(node x);
	vector<double> onBitHypothetical(node x) const;
	void computeFitness(const vector<string>& fitnessNames);
	void save(string filename) const;
	double ics() const;
	
	double sumBLAST() const;
    double sumGOC() const;
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
	unordered_set<node> v1Unaligned; //tracks which nodes are unaligned
	                                 //DANGER: tracks dummy nodes, too!
	bool fitnessValid;
	vector<double> fitness; //all fitnesses stored s.t. larger is better.
	vector<double> fitnessNormalized; //stores fitnesses normalized so
	                                  //that crowded comparison works right.
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

	//returns by how much updateBitscore would change currBitscore
	double hypotheticalBitscoreDelta(node n1, node n2old, node n2new,
		                             bool oldMask, bool newMask) const;
    
    //same stuff for GOC
    double currGOC;
    const GOCDict* gocs;
    void updateGOC(node n1, node n2old, node n2new, bool oldMask,
                   bool newMask);
    double hypotheticalGOCDelta(node n1, node n2old, node n2new,
                                bool oldMask, bool newMask) const;
    
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
	int hypotheticalConservedCountDelta(node n1, node n2old, node n2new, 
		            bool oldMask, bool newMask, node ignore) const;
	void initConservedCount(node n1, node n2, bool mask);
};

