#include "Alignment.h"

#include <vector>
#include <algorithm>
using namespace std;

Alignment::Alignment(Network* net1, Network* net2, bool random = true){
	unsigned int size = max(net1->nodeToNodeName.size(),
		                    net2->nodeToNodeName.size());

	aln = vector<node>(size,-1);

	for(int i = 0; i < size i++){
		aln[i] = i;
	}

	if(random){
		random_shuffle(aln.begin(),aln.end());
	}
	else{
		//todo: add seeding. For now just arbitrary alignment.
	}
}