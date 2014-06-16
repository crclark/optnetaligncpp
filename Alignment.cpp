#include "Alignment.h"
#include "blastinfo.h"
#include "Network.h"

#include <assert.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unordered_set>
#include <limits>
using namespace std;

//todo: note these assume net2 is larger. Ensure that when loading nets.
//this constructor just creates an arbitrary non-random alignment.
//Make it a random one by calling shuf()
Alignment::Alignment(const Network* n1, const Network* n2,
	                 const BLASTDict* bit, const GOCDict* goc){
	unsigned int size = n2->nodeToNodeName.size();
	aln = vector<node>(size,-1);
	alnMask = vector<bool>(size,true);
	fitnessValid = false;
	v1Unaligned = unordered_set<node>();
	v1Unaligned.reserve(n2->nodeToNodeName.size());
	for(int i = 0; i < size; i++){
		aln[i] = i;
	}

	net1 = n1;
	net2 = n2;
	domRank = -1;
	crowdDist = -1.0;
	bitscores = bit;
    gocs = goc;
	actualSize = net1->nodeToNodeName.size();
	//initialize conserved counts for fast ICS computation
	conservedCounts = vector<int>(actualSize);
	for(int i = 0; i < actualSize; i++){
		initConservedCount(i,aln[i], alnMask[i]);
	}	
}

Alignment::Alignment(const Network* n1, const Network* n2, 
	                 string filename,
	                 const BLASTDict* bit, const GOCDict* goc){
	int size = n2->nodeToNodeName.size();
    //todo: nodeNameToNode is sometimes a different size than 
    //nodeToNodeName and I can't figure out why.
	aln = vector<node>(size,-1);
	alnMask = vector<bool>(size,true);
	fitnessValid = false;
	domRank = -1;
	crowdDist = -1.0;
	bitscores = bit;
    gocs = goc;
	net1 = n1;
	net2 = n2;
	actualSize = net1->nodeToNodeName.size();
	ifstream infile(filename);
	if(!infile){
		throw LineReadException(string("Alignment file ")+filename+
			                    string(" failed to open!"));
	}
	string line;
	//keep track of which nodes in V2 aren't aligned so we can
	//add them to the permutation somewhere after loading the file.
	unordered_set<node> v2Unaligned;
	for(int i = 0; i < net2->nodeToNodeName.size(); i++){
		v2Unaligned.insert(i);
	}
    assert(v2Unaligned.size() == aln.size());
    
    int count = 0;
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
		
        if(!net1->nodeNameToNode.count(a)){
            continue;
        }
		u = net1->nodeNameToNode.at(a);
        
		if(!net2->nodeNameToNode.count(b)){
            continue;
		}
		v = net2->nodeNameToNode.at(b);
		
        count++;
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
			v1Unaligned.insert(i);
			v2Unaligned.erase(arbNode);
		}
	}

	//initialize conserved counts for fast ICS computation
	conservedCounts = vector<int>(actualSize);
	for(int i = 0; i < actualSize; i++){
		initConservedCount(i,aln[i], alnMask[i]);
	}
    
	if(bitscores)
		currBitscore = sumBLAST();
        
    if(gocs)
        currGOC = sumGOC();
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
Alignment::Alignment(mt19937& prng, float cxswappb, 
	                        const Alignment& p1,
		                    const Alignment& p2,
		                    bool total){
	const Alignment* par1;
	const Alignment* par2;
	int size = p1.aln.size();
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

	//do standard initialization
	aln = par1->aln;
	alnMask = par1->alnMask;
	fitnessValid = false;
	domRank = -1;
	numThatDominate = -1;
	crowdDist = -1.0;
	bitscores = par1->bitscores;
    gocs = par1->gocs;
	actualSize = par1->actualSize;
	currBitscore = par1->currBitscore;
    currGOC = par1->currGOC;
	conservedCounts = par1->conservedCounts;
	net1 = par1->net1;
	net2 = par1->net2;
	v1Unaligned = par1->v1Unaligned;


	//do crossover
	for(int i = 0; i<size; i++){
		if(fltgen(prng) < cxswappb){
			//swap nodes
			int temp1 = aln[i];
			int temp2 = par2->aln[i];
			aln[i] = temp2;
			aln[par1Indices[temp2]] = temp1;

			bool temp1Mask = alnMask[i];
			bool par1bool = alnMask[par1Indices[temp2]];
			bool par2bool = par2->alnMask[i];
			//swap mask bools
			if(!total){ //if partial alns allowed, either | or & the bit here.
				alnMask[par1Indices[temp2]] = temp1Mask;
				if(coinFlip(prng)){
					alnMask[i] = par1bool && par2bool;
					if(!alnMask[i]){
						v1Unaligned.insert(i);
					}
				}
				else{
					alnMask[i] = par1bool || par2bool;
					if(!alnMask[i]){
						v1Unaligned.insert(i);
					}
				}
			}

			//update currBitscore and conservedCounts
			updateBitscore(i,temp1,temp2,temp1Mask,alnMask[i]);
			updateConservedCount(i, temp1, temp2, temp1Mask, alnMask[i],
				                 -1);
			updateGOC(i,temp1,temp2,temp1Mask,alnMask[i]);
			updateBitscore(par1Indices[temp2],temp2,temp1,
				           par1bool, alnMask[par1Indices[temp2]]);
			updateConservedCount(par1Indices[temp2],temp2,temp1,
				                 par1bool, alnMask[par1Indices[temp2]],i);
			updateGOC(par1Indices[temp2],temp2,temp1,par1bool,
				      alnMask[par1Indices[temp2]]);
			//swap index records 
			int itemp = par1Indices[temp1];
			par1Indices[temp1] = par1Indices[temp2];
			par1Indices[temp2] = itemp; 
		}

	}
	fitnessValid = false;
}

//performs a fast, greedy matching based on bitscore
void Alignment::greedyBitscoreMatch(){
	if(bitscores == nullptr){
		assert(1==2);
	}

	aln = vector<node>(net2->nodeToNodeName.size(), -1);
	v1Unaligned = unordered_set<node>();

	for(int i = 0; i < net1->nodeToNodeName.size(); i++){
		v1Unaligned.insert(i);
	}

	unordered_set<node> v2Unaligned;
	for(int i = 0; i < net2->nodeToNodeName.size(); i++){
		v2Unaligned.insert(i);
	}

	unordered_set<node> v1Aligned, v2Aligned;

	for(auto n1 : v1Unaligned){
		node bestN2 = -1;
		double bestScore = 0.0;
		for(auto n2 : v2Unaligned){
			if(bitscores->count(n1) && bitscores->at(n1).count(n2)
			   && bitscores->at(n1).at(n2) > bestScore ){
				bestN2 = n2;
				bestScore = bitscores->at(n1).at(n2);
			}
			else if(bestN2 == -1){
				bestN2 = n2;
			}
		}

		aln[n1] = bestN2;
		v1Aligned.insert(n1);
		v2Aligned.insert(bestN2);
		v2Unaligned.erase(bestN2);
		//NOTE: can't change v1Unaligned here because we're iterating over it.
		//instead changed in the for loop updating the mask below.
	}

	//update mask
	for(int i = 0; i < alnMask.size(); i++){
		if(v1Aligned.count(i)){
			alnMask[i] = true;
			v1Unaligned.erase(i);
		}
		else{
			alnMask[i] = false;
		}
	}

	//look for unaligned nodes in V1 and align them to unaligned nodes
	//in V2. 
	for(int i = 0; i < aln.size(); i++){
		if(aln[i] == -1){
			node arbNode = *(v2Unaligned.begin());
			aln[i] = arbNode;
			alnMask[i] = false;
			v1Unaligned.insert(i);
			v2Unaligned.erase(arbNode);
		}
	}

	//initialize conserved counts for fast ICS computation
	conservedCounts = vector<int>(actualSize);
	for(int i = 0; i < actualSize; i++){
		initConservedCount(i,aln[i], alnMask[i]);
	}

	currBitscore = sumBLAST();
}

void Alignment::shuf(mt19937& prng, bool uniformsize, 
	                 bool smallstart, bool total){
	shuffle(aln.begin(),aln.end(), prng);
	if(!total){
		if(!uniformsize && !smallstart){
			uniform_int_distribution<int> coinFlip(0,1);
			for(int i =0; i < alnMask.size(); i++){
				if(coinFlip(prng)){
					alnMask[i] = false;
					v1Unaligned.insert(i);
				}
			}
		}
		else if(uniformsize){
			//actually randomly generating number of mask cells to
			//turn OFF because they are already all on.
			uniform_int_distribution<int> sizeDist(1,actualSize-1);
			int numToDeactivate = sizeDist(prng);
			//generate indices to flip on
			uniform_int_distribution<int> indexDist(0,actualSize-1);
			for(int i = 0; i < numToDeactivate; i++){
				int index = indexDist(prng);
				//if index already deactivated, find another
				while(!alnMask[index]){
					index = indexDist(prng);
				}
				alnMask[index] = false;
				v1Unaligned.insert(index);
			}
			//deactivate all unaligned nodes aligned to dummy nodes
			for(int i = actualSize; i < aln.size(); i++){
				alnMask[i] = false;
				v1Unaligned.insert(i);
			}
		}
		else{ //smallsize
			int maxsize = net1->nodeToNodeName.size() / 100;
			if(maxsize == 0){
				maxsize = 1;
			}
			uniform_int_distribution<int> sizeDist(1,maxsize);
			int numToDeactivate = actualSize - sizeDist(prng);
			uniform_int_distribution<int> indexDist(0,actualSize-1);
			for(int i = 0; i < numToDeactivate; i++){
				int index = indexDist(prng);
				while(!alnMask[index]){
					index = indexDist(prng);
				}
				alnMask[index] = false;
				v1Unaligned.insert(index);
			}
			//deactivate all unaligned nodes aligned to dummy nodes
			for(int i = actualSize; i < aln.size(); i++){
				alnMask[i] = false;
				v1Unaligned.insert(i);
			}
		}
	}
	if(bitscores)
		currBitscore = sumBLAST();
        
    if(gocs)
        currGOC = sumGOC();

	//initialize conserved counts for fast ICS computation
	conservedCounts = vector<int>(actualSize);
	for(int i = 0; i < actualSize; i++){
		initConservedCount(i,aln[i], alnMask[i]);
	}	
}

//todo: add secondary mutate op for changing mask

void Alignment::mutate(mt19937& prng, float mutswappb, bool total){
	int size = aln.size();
	uniform_real_distribution<float> fltgen(0.0,1.0);
	//size-2 because of increment below
	uniform_int_distribution<int> intgen(0,size-2); 
	for(int i = 0; i < aln.size(); i++){
		//swap probabilistically
		
		if(fltgen(prng) < mutswappb){
			int swapIndex = intgen(prng);
			if(swapIndex >= i){
				swapIndex++;
			}
			doSwap(i,swapIndex);
		}
		
		//switch mask bit on/off probabilistically
		
		if(!total && fltgen(prng) < mutswappb){
			alnMask[i] = !alnMask[i];
			if(alnMask[i]){
				v1Unaligned.erase(i);
			}
			else{
				v1Unaligned.insert(i);
			}
			updateBitscore(i, aln[i], aln[i], !alnMask[i],
		                   alnMask[i]);
            updateGOC(i,aln[i],aln[i],!alnMask[i],alnMask[i]);
			updateConservedCount(i, aln[i], aln[i], !alnMask[i],
		                         alnMask[i], -1);
		}
		
		
	}
	fitnessValid = false;
}

//takes 2 nodes from V1 and swaps the nodes they are aligned to in V2.
//updates currBitscore accordingly.
//updates conservedCounts accordingly as well.
void Alignment::doSwap(node x, node y){
	node temp = aln[x];
	aln[x] = aln[y];
	aln[y] = temp;

	bool tempb = alnMask[x];
	alnMask[x] = alnMask[y];
	alnMask[y] = tempb;

	if(alnMask[x]){
		v1Unaligned.erase(x);
	}
	else{
		v1Unaligned.insert(x);
	}

	if(alnMask[y]){
		v1Unaligned.erase(x);
	}
	else{
		v1Unaligned.insert(x);
	}

	updateBitscore(x, aln[y], aln[x], alnMask[y], alnMask[x]);
    updateGOC(x, aln[y], aln[x], alnMask[y], alnMask[x]);
	updateConservedCount(x, aln[y], aln[x], alnMask[y], alnMask[x],
		                 -1);
	updateConservedCount(y, aln[x], aln[y], alnMask[x], alnMask[y],
		                 x);
	updateBitscore(y, aln[x], aln[y], alnMask[x], alnMask[y]);
    updateGOC(y, aln[x], aln[y], alnMask[x], alnMask[y]);
}

//todo: test this
vector<double> Alignment::doSwapHypothetical(node x, node y) const{
	vector<double> toReturn;
    
	node hypAlnx = aln[y];
	node hypAlny = aln[x];
    
	bool hypAlnMskx = alnMask[y];
	bool hypAlnMsky = alnMask[x];
    
    double ccDelt = 0.0;
	ccDelt += hypotheticalConservedCountDelta(x, hypAlny, hypAlnx,
	               hypAlnMsky, hypAlnMskx, -1);
	ccDelt += hypotheticalConservedCountDelta(y,hypAlnx,hypAlny,
		           hypAlnMskx, hypAlnMsky, x);
    toReturn.push_back(ccDelt);
	if(bitscores != nullptr){
        double bitDelt = 0.0;
		bitDelt += hypotheticalBitscoreDelta(x, hypAlny, hypAlnx,
			           hypAlnMsky, hypAlnMskx);
		bitDelt += hypotheticalBitscoreDelta(y, hypAlnx, hypAlny,
			           hypAlnMskx, hypAlnMsky);
        toReturn.push_back(bitDelt);
	}
    
    if(gocs != nullptr){
        double gocDelt = 0.0;
        gocDelt += hypotheticalGOCDelta(x, hypAlny, hypAlnx,
                    hypAlnMsky, hypAlnMskx);
        gocDelt += hypotheticalGOCDelta(y, hypAlnx, hypAlny,
                       hypAlnMskx, hypAlnMsky);
        toReturn.push_back(gocDelt);
    }

	return toReturn;
}

void Alignment::onBit(node x){
	bool old = alnMask[x];
	alnMask[x] = true;
	v1Unaligned.erase(x);
	updateBitscore(x,aln[x],aln[x],old, alnMask[x]);
    updateGOC(x, aln[x], aln[x], old, alnMask[x]);
	updateConservedCount(x,aln[x],aln[x],old, alnMask[x],-1);
}

//todo: test this
vector<double> Alignment::onBitHypothetical(node x) const{
	vector<double> toReturn;

	toReturn.push_back(hypotheticalConservedCountDelta(x, aln[x], aln[x],
		           false, true, -1));

	if(bitscores != nullptr){
		toReturn.push_back(hypotheticalBitscoreDelta(x, aln[x], aln[x],
			           false, true));
	}
    
    if(gocs != nullptr){
        toReturn.push_back(hypotheticalGOCDelta(x,aln[x],aln[x],
                           false, true));
    }

	return toReturn;
}

//todo: add support for GO annotations
//todo: maybe something more principled than fitnessNames (so ad hoc!)
void Alignment::computeFitness(const vector<string>& fitnessNames){
	fitness = vector<double>(fitnessNames.size(),0.0);
	for(int i = 0; i < fitnessNames.size(); i++){
		if(fitnessNames.at(i) == "ICS"){
			fitness.at(i) = fastICS(); //ics();
		}
		if(fitnessNames.at(i) == "EC"){
			fitness.at(i) = fastICSNumerator() / 
			                double(net1->edges.size());
		}
		if(fitnessNames.at(i) == "BitscoreSum"){
			fitness.at(i) = currBitscore; //sumBLAST();
		}
		if(fitnessNames.at(i) == "EvalsSum"){
			fitness.at(i) = -1.0*sumBLAST();
		}
		if(fitnessNames.at(i) == "Size"){
			fitness.at(i) = alnSize();
		}
        if(fitnessNames.at(i) == "GOC"){
            fitness.at(i) = currGOC;
        }
	}
	fitnessValid = true;
}

double Alignment::ics() const{


	unordered_set<node> v1Unaligned;
	unordered_set<node> v2Unaligned; //set of unaligned nodes in V2
	for(int i = 0; i < aln.size(); i++){
		if(!alnMask[i]){
			v1Unaligned.insert(i);
			v2Unaligned.insert(aln[i]);
		}
		//second way for a node in V2 to be unaligned: by
		//being aligned to a dummy node.
		if(i >= net1->nodeToNodeName.size()){
			v2Unaligned.insert(aln[i]);
		}
	}

	//induced subgraph of net2 that has been mapped to:
	unordered_set<Edge, EdgeHash> inducedES;
	for(Edge e : net2->edges){
		if(!(v2Unaligned.count(e.u()) || v2Unaligned.count(e.v()))){
			inducedES.insert(e);
		}
	}

	unordered_set<Edge, EdgeHash> mappedNet1ES;
	for(Edge e : net1->edges){
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
		unordered_set<Edge,EdgeHash> intersect;
		if(mappedNet1ES.size() < inducedES.size()){
			for(auto it = mappedNet1ES.begin(); it != mappedNet1ES.end(); it++){
				if(inducedES.count(*it)){
					intersect.insert(*it);
				}
			}
		}
		else{
			for(auto it = inducedES.begin(); it != inducedES.end(); it++){
				if(mappedNet1ES.count(*it)){
					intersect.insert(*it);
				}
			}
		}
		double numerator = double(intersect.size());
		return numerator/denominator;
	}
}

//note: this function is probably as fast as it is going to get. 
//my attempts at optimizing it further have only slowed it down.
//todo: if the speed of this function is still unbearable,
//consider tracking it as we go like numerator.
double Alignment::fastICSDenominator() const{
	//construct set of nodes actually mapped to

	vector<bool> mapped(aln.size(),false);
	for(int i = 0; i < actualSize; i++){
		mapped[aln[i]] = alnMask[i]; 
	}

	int count = 0;
	for(int i = 0; i < mapped.size(); i++){
		if(mapped[i]){
			for(auto y : net2->adjList.at(i)){
				if(mapped[y]){
					count++;
				}
			}
		}
	}

	//self loops are only counted once instead of twice like
	//all other edges, so we add the number of self loops to
	//the count so we count them correctly.
	count += net2->numSelfLoops;

	//cout<<"fast denominator count is "<<count<<endl;

	return double(count / 2);
}

double Alignment::fastICSNumerator() const{
	int sum = 0;
	for(int i = 0; i < conservedCounts.size(); i++){
		if(alnMask[i]){
			sum += conservedCounts[i];
		}
	}
	return double(sum/2);
}

double Alignment::fastICS() const{
	double denom = fastICSDenominator();
	if(denom == 0.0){
		return 0.0;
	}
	return fastICSNumerator() / denom;
}

double Alignment::sumBLAST() const{
	double toReturn = 0.0;

	for(node i = 0; i < actualSize; i++){
		if(bitscores->count(i) && bitscores->at(i).count(aln[i]) &&
		   alnMask[i]){
			toReturn += bitscores->at(i).at(aln[i]);
		}
	}
	return toReturn;
}

double Alignment::sumGOC() const{
    double toReturn = 0.0;
    
    for(node i = 0; i < actualSize; i++){
        if(gocs->count(i) && gocs->at(i).count(aln[i]) && alnMask[i]){
            toReturn += gocs->at(i).at(aln.at(i));
        }
    }
    
    return toReturn;
}

double Alignment::alnSize() const{
	double toReturn = 0.0;
	for(int i = 0; i < actualSize; i++){
		if(alnMask[i]){
			toReturn += 1.0;
		}
	}
	return toReturn;
}

void Alignment::save(string filename) const{
	ofstream ofile(filename);
	for(int i = 0; i < net1->nodeToNodeName.size(); i++){
		if(alnMask[i]){
			string u = net1->nodeToNodeName.at(i);
			string v = net2->nodeToNodeName.at(aln[i]);
			ofile<<u<<' '<<v<<endl;
		}
	}
}

inline void Alignment::updateBitscore(node n1, node n2old, node n2new, 
									  bool oldMask, bool newMask){
	if(!bitscores)
		return;
	if(n1 >= actualSize)
		return;

	double oldScore = 0.0;

	//set oldScore
	if(bitscores->count(n1) && bitscores->at(n1).count(n2old) && oldMask)
		oldScore = bitscores->at(n1).at(n2old);

	double newScore = 0.0;

	//set newScore
	if(bitscores->count(n1) && bitscores->at(n1).count(n2new) && newMask)
		newScore = bitscores->at(n1).at(n2new);

	double delta = newScore - oldScore;
	currBitscore += delta;
}

inline void Alignment::updateGOC(node n1, node n2old, node n2new,
                                 bool oldMask, bool newMask){
    
    currGOC += hypotheticalGOCDelta(n1, n2old, n2new, oldMask, newMask);
}

double Alignment::hypotheticalBitscoreDelta(node n1, node n2old, node n2new,
		                             bool oldMask, bool newMask) const{
	if(!bitscores)
		return 0.0;
	if(n1 >= actualSize)
		return 0.0;

	double oldScore = 0.0;

	//set oldScore
	if(bitscores->count(n1) && bitscores->at(n1).count(n2old) && oldMask)
		oldScore = bitscores->at(n1).at(n2old);

	double newScore = 0.0;

	//set newScore
	if(bitscores->count(n1) && bitscores->at(n1).count(n2new) && newMask)
		newScore = bitscores->at(n1).at(n2new);

	return newScore - oldScore;	
}

double Alignment::hypotheticalGOCDelta(node n1, node n2old, node n2new,
                                       bool oldMask, bool newMask) const
{
    if(!gocs)
        return 0.0;
    if(n1 >= actualSize)
        return 0.0;
        
    double oldScore = 0.0;
    
    if(gocs->count(n1) && gocs->at(n1).count(n2old) && oldMask)
        oldScore = gocs->at(n1).at(n2old);
        
    double newScore = 0.0;
    
    if(gocs->count(n1) && gocs->at(n1).count(n2new) && newMask)
        newScore = gocs->at(n1).at(n2new);
        
    return newScore - oldScore;
}

//ignore parameter is for the edge case where two neighbors are 
//swapped with each other. In that case, one of them already has
//an up-to-date conservedCount and so it becomes incorrect when it
//is touched again while its neighbor is being updated. Thus, we
//need to specify that it should be ignored.
void Alignment::updateConservedCount(node n1, node n2old, node n2new, 
	                                 bool oldMask, bool newMask,
	                                 node ignore){
	
	//easy cases first:
	//check that n1 is not a dummy node:
	if(n1 >= actualSize){
		return;
	}

	if(n2old == n2new && oldMask == newMask){
		return;
	}


	conservedCounts[n1] = 0;
	
	//todo: degree 0 currently impossible but may be needed later
	/*
	if(net1->degree(n1) == 0){
		return;
	}
	*/

	if(!newMask && !oldMask){
		//nothing has changed. We were unaligned, we are still unaligned.
		return;
	}
	else if(!newMask && oldMask){
		//our bit has been turned off.
		//if our neighbors had a conserved edge thanks to us,
		//we need to decrement their count.
		for(auto i : net1->adjList.at(n1)){
			if(net2->adjMatrix[aln[i]][n2old]){
				if(conservedCounts[i] > 0){
					conservedCounts[i]--;
				}
			}
		}
		return; //don't fall through to counting conservations that dont exist
	}

	for(auto i : net1->adjList.at(n1)){
		//need to count self loops one extra time so as to
		//ensure they are counted correctly when conserved.
		//otherwise, they will only be counted once
		//when all other conserved edges get counted twice
		//(which is why we divide by two in fastICSNumerator())
		if(n1 == i){
			//NOTE: at this point conservedCounts[n1] is zero!
			//cout<<"in self-loop case!"<<endl;
			if(newMask && net2->adjMatrix[n2new][n2new]){
				conservedCounts[n1]++;
			}
			else {
				continue;
			}
		}

		//todo: we still get small errors in ICS due to the presence of self-
		//loops (ignoring self-loops when loading network makes them go away).

		if(i == ignore){
			//cout<<"i is ignore; continuing"<<endl;
			continue;
		}
		if(!alnMask[i]){
			//cout<<"mask off; continuing"<<endl;
			continue;
		}
		if(net2->adjMatrix[n2new][aln[i]]){
			//cout<<"aln[i] is neighbor to n2new; edge conserved. n1 conserved "
			//      "count++"<<endl;
			conservedCounts[n1]++;
		}
		//conservedCount of i either increases by 1,
		//decreases by 1, or remains unchanged.

		//unchanged case: both n2old and n2new are 
		//neighbors of aln[i], or neither are.
		auto alnINbrs = &(net2->adjMatrix[aln[i]]); 
		if(((*alnINbrs)[n2old] && 
		    (*alnINbrs)[n2new] && oldMask == newMask)
		   || (!(*alnINbrs)[n2old] &&
		   	   !(*alnINbrs)[n2new])){
			//cout<<"both n2old and n2new are neighbors of aln[i], or neither is."
		    //   <<endl;
			//do nothing
		}
		else if((*alnINbrs)[n2old] && (*alnINbrs)[n2new]
			    && !oldMask && newMask){
			//cout<<"both n2old and n2new are neighbors and mask went from 0 to 1"
		     //   <<endl;
			conservedCounts[i]++;
		}
		else if((*alnINbrs)[n2old] &&
			    !(*alnINbrs)[n2new]){
			if(conservedCounts[i]>0){
			//	cout<<"n1 to i used to be conserved and is no longer"<<endl;
				conservedCounts[i]--;
			}
		}
		else if(!(*alnINbrs)[n2old] &&
			    (*alnINbrs)[n2new]){
			//cout<<"n1 to i wasn't conserved but is now"<<endl;
			conservedCounts[i]++;
		}
	}
}

int Alignment::hypotheticalConservedCountDelta(node n1, node n2old, node n2new, 
	                    bool oldMask, bool newMask, node ignore) const{
	//easy cases first:
	//check that n1 is not a dummy node:
    //cout<<"first cond"<<endl;
	if(n1 >= actualSize){
		return 0;
	}

    //cout<<"second cond"<<endl;
	if(n2old == n2new && oldMask == newMask){
		return 0;
	}


	int toReturn = 0;
	
	//todo: degree 0 currently impossible but may be needed later
	/*
	if(net1->degree(n1) == 0){
		return;
	}
	*/
    
    //cout<<"nochange cond"<<endl;
	if(!newMask && !oldMask){
		//nothing has changed. We were unaligned, we are still unaligned.
		return 0;
	}
	else if(!newMask && oldMask){
		//our bit has been turned off.
		//if our neighbors had a conserved edge thanks to us,
		//we need to decrement their count.
		for(auto i : net1->adjList.at(n1)){
			if(net2->adjMatrix[aln[i]][n2old]){
				if(conservedCounts[i] > 0){
					toReturn--;
				}
			}
		}
		return toReturn; //don't fall through to counting conservations that dont exist
	}
    
    //cout<<"for neighbors"<<endl;
	for(auto i : net1->adjList.at(n1)){
		//need to count self loops one extra time so as to
		//ensure they are counted correctly when conserved.
		//otherwise, they will only be counted once
		//when all other conserved edges get counted twice
		//(which is why we divide by two in fastICSNumerator())
        //cout<<"trying n1==i"<<endl;
		if(n1 == i){
            //cout<<"n1 == i"<<endl;
			//NOTE: at this point conservedCounts[n1] is zero!
			//cout<<"in self-loop case!"<<endl;
			if(newMask && net2->adjMatrix[n2new][n2new]){
				toReturn++;
			}
			else {
				continue;
			}
		}

		//todo: we still get small errors in ICS due to the presence of self-
		//loops (ignoring self-loops when loading network makes them go away).
        //cout<<"trying i == ignore"<<endl;
		if(i == ignore){
		//	cout<<"i is ignore; continuing"<<endl;
			continue;
		}
        //cout<<"alnMask[i]"<<endl;
		if(!alnMask[i]){
			//cout<<"mask off; continuing"<<endl;
			continue;
		}
        
		if(net2->adjMatrix[n2new][aln[i]]){
			//cout<<"aln[i] is neighbor to n2new; edge conserved. n1 conserved "
			 //     "count++"<<endl;
			toReturn++;
		}
		//conservedCount of i either increases by 1,
		//decreases by 1, or remains unchanged.

		//unchanged case: both n2old and n2new are 
		//neighbors of aln[i], or neither are.
		auto alnINbrs = &(net2->adjMatrix[aln[i]]); 
		if(((*alnINbrs)[n2old] && 
		    (*alnINbrs)[n2new] && oldMask == newMask)
		   || (!(*alnINbrs)[n2old] &&
		   	   !(*alnINbrs)[n2new])){
			//cout<<"both n2old and n2new are neighbors of aln[i], or neither is."
		    //   <<endl;
			//do nothing
		}
		else if((*alnINbrs)[n2old] && (*alnINbrs)[n2new]
			    && !oldMask && newMask){
			//cout<<"both n2old and n2new are neighbors and mask went from 0 to 1"
		     //   <<endl;
			toReturn++;
		}
		else if((*alnINbrs)[n2old] &&
			    !(*alnINbrs)[n2new]){
			if(conservedCounts[i]>0){
				//cout<<"n1 to i used to be conserved and is no longer"<<endl;
				toReturn--;
			}
		}
		else if(!(*alnINbrs)[n2old] &&
			    (*alnINbrs)[n2new]){
			//cout<<"n1 to i wasn't conserved but is now"<<endl;
			toReturn++;
		}
	}

	return toReturn;
}

void Alignment::initConservedCount(node n1, node n2, bool mask){
	conservedCounts[n1] = 0;
	if(!mask){
		return;
	}
	if(net1->degree(n1) == 0){
		return;
	}
	for(auto i : net1->adjList.at(n1)){
		//deal with self loops
		if(n1 == i){
			if(net2->adjList.at(n2).count(n2)){
				conservedCounts[n1]++;
			}
		}
		if(!alnMask[i]){
			continue;
		}
		if(net2->adjList.at(n2).count(aln[i])){
			conservedCounts[n1]++;
		}
	}
}
