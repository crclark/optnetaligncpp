#pragma once

#include <vector>
#include <algorithm>

#include "Network.h"

class Alignment{
public:	
	Alignment(Network* net1, Network* net2, bool random = true);

	vector<node> aln;
};