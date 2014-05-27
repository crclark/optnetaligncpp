#pragma once

#include <vector>
#include <algorithm>
#include <random>
#include "Network.h"

class Alignment{
public:	
	Alignment(Network* net1, Network* net2, bool random = true);
	void mutate(mt19937& prng, float mutswappb);
	vector<node> aln;
};