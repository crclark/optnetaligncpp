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
	vector<node> aln;
	bool fitnessValid;
	vector<double> fitness; //all fitnesses stored s.t. larger is better.

	double ics(const Network& net1, const Network& net2) const;
	double sumBLAST(const Network& net1,
		            const Network& net2,
		            const BLASTDict d) const;
};