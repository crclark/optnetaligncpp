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
Alignment::Alignment(const Network& net1, const Network& net2){
	unsigned int size = net2.nodeToNodeName.size();
	aln = vector<node>(size,-1);
	alnMask = vector<bool>(size,true);
	fitnessValid = false;
	
	for(int i = 0; i < size; i++){
		aln[i] = i;
	}

	domRank = -1;
	crowdDist = -1.0;
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
	unordered_set<node> v2Unaligned;
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

void Alignment::shuf(mt19937& prng){
	shuffle(aln.begin(),aln.end(), prng);
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


//todo: add support for GO annotations
//todo: maybe something more principled than fitnessNames (so ad hoc!)
void Alignment::computeFitness(const Network& net1,
	                const Network& net2,
	                const BLASTDict& bitscores,
	                const BLASTDict& evalues,
	                const vector<string>& fitnessNames){
	fitness = vector<double>(fitnessNames.size(),0.0);
	for(int i = 0; i < fitnessNames.size(); i++){
		if(fitnessNames.at(i) == "ICS"){
			fitness.at(i) = ics(net1,net2);
		}
		if(fitnessNames.at(i) == "BitscoreSum"){
			fitness.at(i) = sumBLAST(net1,net2,bitscores);
		}
		if(fitnessNames.at(i) == "EvalsSum"){
			fitness.at(i) = -1.0*sumBLAST(net1,net2,evalues);
		}
	}
	fitnessValid = true;
}

double Alignment::ics(const Network& net1, const Network& net2) const{


	unordered_set<node> v1Unaligned;
	unordered_set<node> v2Unaligned; //set of unaligned nodes in V2
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
	unordered_set<Edge, EdgeHash> inducedES;
	for(Edge e : net2.edges){
		if(!(v2Unaligned.count(e.u()) || v2Unaligned.count(e.v()))){
			inducedES.insert(e);
		}
	}

	unordered_set<Edge, EdgeHash> mappedNet1ES;
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

double Alignment::sumBLAST(const Network& net1,
		            const Network& net2,
		            const BLASTDict& d) const{
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


//this function assumes fitnesses have already been assigned
//todo: add check that that's the case.
vector<vector<Alignment*> > nonDominatedSort(const vector<Alignment*>& in){
	vector<vector<Alignment*> > fronts(1); //know there is at least one front.
	for(int i = 0; i < in.size(); i++){
		in[i]->numThatDominate = 0;
		in[i]->dominated.clear();
		for(int j = 0; j < in.size(); j++){
			if(i == j){
				continue;
			}
			if(dominates(in[i],in[j])){
				in[i]->dominated.insert(in[j]);
			}
			else if(dominates(in[j],in[i])){
				in[i]->numThatDominate++;
			}
		}

		if(in[i]->numThatDominate == 0){
			in[i]->domRank = 0;
			fronts[0].push_back(in[i]);
		}
	}

	int i = 0;
	while(!(fronts.size() == i || fronts[i].empty())){
		vector<Alignment*> nextFront;

		for(int j = 0; j < fronts[i].size(); j++){
			for(auto q : fronts[i][j]->dominated){
				q->numThatDominate--;
				if(q->numThatDominate == 0){
					q->domRank = i+1;
					nextFront.push_back(q);
				}
			}
		}
		i++;
		if(!nextFront.empty()){
			fronts.push_back(nextFront);
		}
	}

	int frontsTotal = 0;
	for(int i = 0; i<fronts.size(); i++){
		frontsTotal += fronts[i].size();
	}

	assert(frontsTotal == in.size());

	return fronts;

}


//takes a front as input and assigns crowdDist to each element
//note: results meaningless if input is not non-dominated set
void setCrowdingDists(vector<Alignment*>& in){
	/*
	//debug check: no one in input should dominate anyone else
	for(auto i : in){
		for(auto j : in){
			if(dominates(i,j)){
				cout<<"one aln in a front dominates another!"<<endl;
				cout<<"All fitnesses in this front:"<<endl;
				for(auto x : in){
					for(int y = 0; y<x->fitness.size(); y++){
						cout<<x->fitness[y]<<" ";
					}
					cout<<endl;
				}
				assert(0==1);
			}
		}
	}*/

	//init all to zero
	for(auto i : in){
		i->crowdDist = 0.0;
	}

	int numObjs = in[0]->fitness.size();
	int numAlns = in.size();

	//for each objective m
	for(int m = 0; m < numObjs; m++){
		//sort by objective m
		sort(in.begin(), in.end(),
			[m](const Alignment* a, const Alignment* b){
				return a->fitness[m] < b->fitness[m];
			});

		//set boundary points to max dist
		in[0]->crowdDist = numeric_limits<double>::max();
		in[numAlns-1]->crowdDist = numeric_limits<double>::max();

		//increment crowding dist for the current objective
		for(int i = 1; i< (numAlns-1); i++){
			double numerator = in[i+1]->fitness[m] - in[i-1]->fitness[m];
			double denom = in[numAlns-1]->fitness[m] - in[0]->fitness[m];
			in[i]-> crowdDist += (numerator/denom);
		}

	}
}

//returns true if aln1 Pareto dominates aln2
bool dominates(Alignment* aln1, Alignment* aln2){
	/*
	if(!(aln1->fitnessValid && aln2->fitnessValid)){
		cout<<"danger: dominates called on alns with invalid fitness"<<endl;
	}
	if(aln1->fitness.size() != aln2->fitness.size()){
		cout<<"danger: fitness sizes not equal!"<<endl;
	}
	*/
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

bool crowdedComp(Alignment* aln1, Alignment* aln2){
	if(aln1->domRank == -1 || aln2->domRank == -1){
		cout<<"Danger: domRank uninitialized in crowdedComp!"<<endl;
	}
	if(aln1->crowdDist < 0.0 || aln2->crowdDist < 0.0){
		cout<<"Danger: crowdDist uninitialized in crowdedComp!"<<endl;
	}
	return (aln1->domRank < aln2->domRank)
	        || (aln1->domRank == aln2->domRank &&
	        	aln1->crowdDist > aln2->crowdDist);
}


//preconditions: tournSize smaller than in
//all in elems have crowdDist and domCount calculated
//returns two alignment pointers
vector<Alignment*> binSel(mt19937& prng,
	                      const vector<Alignment*>& in, 
	                      unsigned int tournSize){
	vector<unsigned int> indices(in.size());
	
	for(int i = 0; i < in.size(); i++){
		indices[i] = i;
	}

	shuffle(indices.begin(),indices.end(), prng);

	//grab the best of a random subset of in
	sort(indices.begin(),indices.begin()+tournSize,
		[&in](unsigned int a, unsigned int b){
			return crowdedComp(in.at(a),in.at(b));
		});

	vector<Alignment*> toReturn;
	toReturn.push_back(in.at(indices[0]));
	toReturn.push_back(in.at(indices[1]));

	return toReturn;
}

void reportStats(const vector<Alignment*>& in){

	for(int i =0; i < in[0]->fitness.size(); i++){
		double sum = 0.0;
		double max = 0.0;
		double mean;

		for(auto p : in){
			double temp = p->fitness[i];
			if(temp > max)
				max = temp;
			sum += p->fitness[i];
		}
		cout<<"Max of objective "<<i<<" is "<<max<<endl;
		mean = sum/double(in.size());
		cout<<"Mean of objective "<<i<<" is "<<mean<<endl;

		double std_dev = 0.0;

		for(auto p : in){
			double temp = p->fitness[i] - mean;
			std_dev += temp*temp;
		}

		std_dev /= double(in.size());
		std_dev = sqrt(std_dev);
		cout<<"Std. Dev. of objective "<<i<<" is "<<std_dev<<endl;
	}

}