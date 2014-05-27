#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include "Network.h"

class Alignment{
public:	
	Alignment(const Network& net1, const Network& net2, bool random = true);
	Alignment(const Network& net1, const Network& net2, string filename);
	void mutate(mt19937& prng, float mutswappb);
	void becomeChild(mt19937& prng, float cxswappb, 
		             const Alignment& p1,
		             const Alignment& p2);
	vector<node> aln;
	bool fitnessValid;
	vector<double> fitness;
};