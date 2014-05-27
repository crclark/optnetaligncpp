#include "Alignment.h"

#include <vector>
#include <algorithm>
#include <random>
using namespace std;

Alignment::Alignment(Network* net1, Network* net2, bool random){
	unsigned int size = max(net1->nodeToNodeName.size(),
		                    net2->nodeToNodeName.size());

	aln = vector<node>(size,-1);

	for(int i = 0; i < size; i++){
		aln[i] = i;
	}

	if(random){
		random_shuffle(aln.begin(),aln.end());
	}
	else{
		//todo: add seeding. For now just arbitrary alignment.
	}
}

void Alignment::mutate(mt19937& prng, float mutswappb){
	int size = aln.size();
	uniform_real_distribution<float> fltgen(0.0,1.0);
	uniform_int_distribution<int> intgen(0,size-2); 
	for(int i = 0; i < aln.size(); i++){
		if(fltgen(prng) < mutswappb){
			int swapIndex = intgen(prng);
			if(swapIndex >= i){
				swapIndex++;
			}
			int temp = aln[swapIndex];
			aln[swapIndex] = aln[i];
			aln[i] = temp;
		}
	}
}