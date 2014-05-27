#include "Alignment.h"

#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;

Alignment::Alignment(const Network& net1, const Network& net2, bool random){
	unsigned int size = max(net1.nodeToNodeName.size(),
		                    net2.nodeToNodeName.size());

	aln = vector<node>(size,-1);
	fitnessValid = false;
	
	for(int i = 0; i < size; i++){
		aln[i] = i;
	}

	if(random){
		random_shuffle(aln.begin(),aln.end());
	}
	else{
		//todo: add seeding. For now just arbitrary alignment.
		random_shuffle(aln.begin(),aln.end());
	}
}

Alignment::Alignment(const Network& net1, const Network& net2, 
	                 string filename){
	int size = net1.nodeNameToNode.size();
	aln = vector<node>(size,-1);
	fitnessValid = false;

	ifstream infile(filename);

	string line;
	while(getline(infile,line)){
		istringstream iss(line);
		string a, b;
		node u, v;

		if (!(iss >> a >> b)){
			u = net1.nodeNameToNode.at(a);
			v = net2.nodeNameToNode.at(b);
			aln[u] = v;
		}
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

/*
This does UPMX crossover in a (hopefully) more efficient way. An existing
alignment becomes the child of two alignments. Standard UPMX modifies the
two existing alignments in-place. Since our algorithm is elitist, however,
we need those old alignments. Instead, the given alignment becomes what one
of those two resulting alignment would be without modifying the parents.
Which modified parent to become is chosen at random. 
*/
void Alignment::becomeChild(mt19937& prng, float cxswappb, 
	                        const Alignment& p1,
		                    const Alignment& p2){
	const Alignment* par1;
	const Alignment* par2;
	int size = aln.size();
	uniform_int_distribution<int> coinFlip(0,1);
	uniform_real_distribution<float> fltgen(0.0,1.0);
	if(coinFlip(prng)){
		par1 = &p1;
		par2 = &p2;
	}
	else{
		par1 = &p2;
		par2 = &p1;
	}

	vector<node> par1Indices(size,-1);
	for(int i = 0; i < size; i++){
		par1Indices[par1->aln[i]] = i;
	}

	aln = par1->aln;

	for(int i = 0; i<size; i++){
		if(fltgen(prng) < cxswappb){
			//swap nodes
			int temp1 = aln[i];
			int temp2 = par2->aln[i];
			aln[i] = temp2;
			aln[par1Indices[temp2]] = temp1;

			//swap index records 
			int itemp = par1Indices[temp1];
			par1Indices[temp1] = par1Indices[temp2];
			par1Indices[temp2] = itemp; 
		}

	}
}