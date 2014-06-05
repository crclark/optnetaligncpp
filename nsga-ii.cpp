#include "nsga-ii.h"

#include "Alignment.h"
#include "Network.h"
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <iostream>
using namespace std;

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

	return fronts;

}


//takes a front as input and assigns crowdDist to each element
//note: results meaningless if input is not non-dominated set
void setCrowdingDists(vector<Alignment*>& in){

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

void reportStats(const vector<Alignment*>& in, 
	             const vector<string> fitnessNames, 
	             bool verbose){

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

		mean = sum/double(in.size());
		double std_dev = 0.0;

		for(auto p : in){
			double temp = p->fitness[i] - mean;
			std_dev += temp*temp;
		}

		std_dev /= double(in.size());
		std_dev = sqrt(std_dev);
		if(verbose){
			cout<<"Max of objective "<<fitnessNames[i]<<" is "<<max<<endl;
			cout<<"Mean of objective "<<fitnessNames[i]<<" is "<<mean<<endl;
			cout<<"Std. Dev. of objective "<<fitnessNames[i]
			    <<" is "<<std_dev<<endl;
		}
		else{
			cout<<'\t'<<max<<'\t'<<mean<<'\t'<<std_dev;
		}
	}

}