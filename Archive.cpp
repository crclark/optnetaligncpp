#include "Archive.h"
#include "Network.h"
#include "nsga-ii.h"
#include "Alignment.h"

#include <vector>
#include <set>
#include <algorithm>
#include <random>

//any aln given to Archive can be deleted at any time!
void Archive::insert(Alignment* aln){

	vector<Alignment*> dominatedByAln;

	bool alnIsDominated = false;

	for(auto p : nonDominatedSet){

		if(dominates(aln->fitness,p->fitness)){
			dominatedByAln.push_back(p);
		}
		else if(dominates(p->fitness,aln->fitness)){
			alnIsDominated = true;
			break;
		}
	}

	if(alnIsDominated){
		delete aln;
		return;
	}
	else{
		for(auto p : dominatedByAln){
			nonDominatedSet.erase(p);
			delete p;
		}
		nonDominatedSet.insert(aln);
		return;
	}
}

void Archive::shrinkToSize(int size){

	if(size >= nonDominatedSet.size()){
		return;
	}

	vector<Alignment*> nonDominatedVector;
	copy(nonDominatedSet.begin(), nonDominatedSet.end(), back_inserter(nonDominatedVector));

	normalizeFitnesses(nonDominatedVector);
	setCrowdingDists(nonDominatedVector);

	sort(nonDominatedVector.begin(), nonDominatedVector.end(), 
		[](Alignment* p1, Alignment* p2){
			return p1->crowdDist > p2->crowdDist;
		}
	);

	set<Alignment*> newSet;

	for(int i = 0; i < nonDominatedVector.size(); i++){
		if(i < size){
			newSet.insert(nonDominatedVector[i]);
		}
		else{
			delete nonDominatedVector[i];
		}
	}

	nonDominatedSet = newSet;

}