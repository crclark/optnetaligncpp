#include "Alignment.h"
#include "blastinfo.h"
#include "Network.h"

#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <set>
using namespace std;

//todo: note these assume net2 is larger. Ensure that when loading nets.
Alignment::Alignment(const Network& net1, const Network& net2){
	unsigned int size = net2.nodeToNodeName.size();
	aln = vector<node>(size,-1);
	alnMask = vector<bool>(size,true);
	fitnessValid = false;
	
	for(int i = 0; i < size; i++){
		aln[i] = i;
	}

	random_shuffle(aln.begin(),aln.end());
}

Alignment::Alignment(const Network& net1, const Network& net2, 
	                 string filename){
	int size = net2.nodeNameToNode.size();
	aln = vector<node>(size,-1);
	alnMask = vector<bool>(size,true);
	fitnessValid = false;
	ifstream infile(filename);
	cout<<"loading from "<<filename<<endl;
	if(!infile){
		throw LineReadException(string("Alignment file ")+filename+
			                    string(" failed to open!"));
	}
	string line;
	//keep track of which nodes in V2 aren't aligned so we can
	//add them to the permutation somewhere after loading the file.
	set<node> v2Unaligned;
	for(int i = 0; i < net2.nodeToNodeName.size(); i++){
		v2Unaligned.insert(i);
	}

	//process each line
	while(getline(infile,line)){
		//cout<<"parsing line: "<<line<<endl;
		istringstream iss(line);
		string a, b;
		node u, v;

		

		if (!(iss >> a >> b)){
			throw LineReadException(string("Parse error in network: ") + filename +
				                   string("on line: ") + line + "\n");
		
		}
		
		u = net1.nodeNameToNode.at(a);
		if(!net2.nodeNameToNode.count(b)){
			cout<<"node "<<b<<" not found in net2!"<<endl;
		}
		v = net2.nodeNameToNode.at(b);
		
		aln[u] = v;
		v2Unaligned.erase(v);
		//cout<<"node "<<a<<" (int: "<<u<<") aligned to node "
		//    <<b<<" (int: "<<v<<")"<<endl;
	
	}

	//look for unaligned nodes in V1 and align them to unaligned nodes
	//in V2. 
	for(int i = 0; i < aln.size(); i++){
		if(aln[i] == -1){
			node arbNode = *(v2Unaligned.begin());
			aln[i] = arbNode;
			alnMask[i] = false;
			v2Unaligned.erase(arbNode);
		}
	}
}

//todo: add secondary mutate op for changing mask

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

			bool tempb = alnMask[swapIndex];
			alnMask[swapIndex] = alnMask[i];
			alnMask[i] = tempb;
		}
	}
	fitnessValid = false;
}

/*
This does UPMX crossover in a (hopefully) more efficient way. An existing
alignment becomes the child of two alignments. Standard UPMX modifies the
two existing alignments in-place. Since our algorithm is elitist, however,
we need those old alignments. Instead, the given alignment becomes what one
of those two resulting alignment would be without modifying the parents.
Which modified parent to become is chosen at random.

todo: add secondary crossover op for changing mask. Or something.
todo: this currently ignores the mask of the parents and maintains
      the mask of the current alignment.
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

			//swap mask bools
			//todo: add tests that this is right
			bool tempb = alnMask[i];
			alnMask[i] = par2->alnMask[i];
			alnMask[par1Indices[temp2]] = tempb;

			//swap index records 
			int itemp = par1Indices[temp1];
			par1Indices[temp1] = par1Indices[temp2];
			par1Indices[temp2] = itemp; 
		}

	}
	fitnessValid = false;
}

void Alignment::computeFitness(const Network& net1,
	                const Network& net2,
	                const BLASTDict d,
	                const vector<string> fitnessNames, 
		            const vector<double> fitnessWeights){
	fitness = vector<double>(fitnessNames.size(),0.0);
	for(int i = 0; i < fitnessNames.size(); i++){
		if(fitnessNames.at(i) == "ICS"){
			fitness.at(i) = ics(net1,net2);
		}
		if(fitnessNames.at(i) == "BitscoreSum"){
			fitness.at(i) = sumBLAST(net1,net2,d);
		}
		if(fitnessNames.at(i) == "EvalsSum"){
			fitness.at(i) = -1.0*sumBLAST(net1,net2,d);
		}
	}
	fitnessValid = true;
}

double Alignment::ics(const Network& net1, const Network& net2) const{


	set<node> v1Unaligned;
	set<node> v2Unaligned; //set of unaligned nodes in V2
	for(int i = 0; i < aln.size(); i++){
		if(!alnMask[i]){
			v1Unaligned.insert(i);
			v2Unaligned.insert(aln[i]);
		}
		//second way for a node in V2 to be unaligned: by
		//being aligned to a dummy node.
		if(i > net1.nodeToNodeName.size()){
			v2Unaligned.insert(aln[i]);
		}
	}

	//induced subgraph of net2 that has been mapped to:
	set<Edge> inducedES;
	for(Edge e : net2.edges){
		if(!(v2Unaligned.count(e.u()) || v2Unaligned.count(e.v()))){
			inducedES.insert(e);
		}
	}

	set<Edge> mappedNet1ES;
	for(Edge e : net1.edges){
		if(!(v1Unaligned.count(e.u()) || v1Unaligned.count(e.v()))){
			Edge mapped(aln[e.u()],aln[e.v()]);
			mappedNet1ES.insert(mapped);
		}
	}

	double denominator = double(inducedES.size());
	if(denominator == 0.0){
		return 0.0;
	}
	else{
		set<Edge> intersect;
		set_intersection(inducedES.begin(), inducedES.end(),
			             mappedNet1ES.begin(), mappedNet1ES.end(),
			             inserter(intersect,intersect.begin()));
		double numerator = double(intersect.size());
		return numerator/denominator;
	}
}

double Alignment::sumBLAST(const Network& net1,
		            const Network& net2,
		            const BLASTDict d) const{
	double toReturn = 0.0;
	for(node i = 0; i < net1.nodeNameToNode.size(); i++){
		if(d.count(i) && d.at(i).count(aln[i]) &&
		   alnMask[i]){
			toReturn += d.at(i).at(aln[i]);
		}
	}

	return toReturn;
}

void Alignment::save(const Network& net1,
	                 const Network& net2,
	                 string filename) const{
	ofstream ofile(filename);
	for(int i = 0; i < net1.nodeToNodeName.size(); i++){
		if(alnMask[i]){
			string u = net1.nodeToNodeName.at(i);
			string v = net2.nodeToNodeName.at(aln[i]);
			ofile<<u<<' '<<v<<endl;
		}
	}
}

vector<vector<Alignment*> > nonDominatedSort(vector<Alignment*> in){
	vector<Alignment*> temp = in;
	sort(temp.begin(), temp.end(), dominates);

	//do domination counts (naively-- O(n^2))
	vector<int> domCount(temp.size(),0);
	for(int i = 0; i < temp.size(); i++){
		for(int j = 0; j < i; j++){
			if(dominates(temp[j], temp[i])){
				domCount[i]++;
			}
		}
	}

	//use domination counts to gather results
	vector<vector<Alignment*> > toReturn;
	int lastDomCount = -1;
	for(int i = 0; i < temp.size(); i++){
		if(domCount[i] == lastDomCount){
			toReturn.back().push_back(temp[i]);
		}
		else{
			toReturn.push_back(vector<Alignment*>());
			lastDomCount = domCount[i];
			toReturn.back().push_back(temp[i]);
		}
	}

	return toReturn;
}

bool dominates(Alignment* aln1, Alignment* aln2){

	bool oneBigger = false;

	for(int i = 0; i < aln1->fitness.size(); i++){
		if(aln1->fitness[i] < aln2->fitness[i]){
			return false;
		}
		if(aln1->fitness[i] > aln2->fitness[i]){
			oneBigger = true;
		}
	}

	return oneBigger;
}