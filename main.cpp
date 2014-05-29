#include <iostream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <unordered_set>
#include <assert.h>
#include <boost/program_options.hpp>

#include "Alignment.h"
#include "argParsing.h"
using namespace std;
namespace po = boost::program_options;



int main(int ac, char* av[])
{
	try{
		argRetVals vals = handleArgs(ac,av);
		po::variables_map vm = get<0>(vals);
		Network* net1 = get<1>(vals);
		Network* net2 = get<2>(vals);
		BLASTDict bitscores = get<3>(vals);
		BLASTDict evalues = get<4>(vals);
		vector<string> fitnessNames = get<5>(vals);

		const int nthreads = vm["nthreads"].as<int>();
		const float mutswappb = vm["mutswappb"].as<float>();
		const float cxswappb = vm["cxswappb"].as<float>();
		//initialize population
		cout<<"creating initial population"<<endl;
		const unsigned int popsize = vm["popsize"].as<int>();
		vector<Alignment*> pop;
		for(int i = 0; i < popsize; i++){	
			Alignment* aln = new Alignment(*net1,*net2);
			aln->shuf();
			aln->computeFitness(*net1,*net2,bitscores,evalues,fitnessNames);
			pop.push_back(aln);
		}

		cout<<"creating initial children"<<endl;
		vector<Alignment*> kids;
		for(int i = 0; i < popsize; i++){
			Alignment* aln = new Alignment(*net1,*net2);
			aln->shuf();
			aln->computeFitness(*net1,*net2,bitscores,evalues,fitnessNames);
			kids.push_back(aln);
		}

		//main loop
		cout<<"starting main loop"<<endl;
		int generations = vm["generations"].as<int>();
		for(int gen = 0; gen < generations; gen++){

			//combinedPtrs is R_t from Deb et al. 2002

			vector<Alignment*> combinedPtrs(2*popsize);
			for(int i = 0; i < popsize; i++){
				combinedPtrs[i] = pop[i];
			}
			for(int i = popsize; i < 2*popsize; i++){
				Alignment* ptr = kids.at(i-popsize);
				combinedPtrs.at(i) = ptr;
			}
			
			vector<vector<Alignment*> > fronts = nonDominatedSort(combinedPtrs);


			//started with best front, add to new population front-by-front
			unordered_set<Alignment*> popNew;
			popNew.reserve(popsize);
			int i = 0;
			while(popNew.size() + fronts[i].size() < popsize){
				setCrowdingDists(fronts[i]);
				popNew.insert(fronts[i].begin(), fronts[i].end());
				i++;
			}

			//add the least-crowded members of the front that doesn't
			//completely fit given our popsize.
			int numLeftToInsert = popsize - popNew.size();
			if(numLeftToInsert > 0){
				assert(fronts[i].size() >= numLeftToInsert);
				sort(fronts[i].begin(),fronts[i].end(),crowdedComp);
				popNew.insert(fronts[i].begin(), 
				              fronts[i].begin() + numLeftToInsert);
			}

			//go through combinedPtrs, deleting those not in popNew
			for(int i = 0; i < combinedPtrs.size(); i++){
				if(!popNew.count(combinedPtrs[i])){
					delete combinedPtrs[i];
				}
			}
			//set pop = popNew
			assert(pop.size() == popNew.size());
			copy(popNew.begin(), popNew.end(), pop.begin());

			//create new kids.
			kids.clear();

			if(nthreads == 1){
				mt19937 g(14);

				for(int i = 0; i < popsize; i++){
					Alignment * aln = new Alignment(*net1,*net2);
					vector<Alignment*> parents = binSel(g, pop, (popsize/10));
					aln->becomeChild(g, cxswappb, *parents[0], *parents[1]);
					aln->mutate(g, mutswappb);
					aln->computeFitness(*net1,*net2,bitscores,evalues,fitnessNames);
					kids.push_back(aln);
				}

			}
			else{ //do multithreaded version of kids creation

			}

			//idea: create a constructor for a non-random arbitrary aln
			//allocate kids using that constructor. Then,
			//in parallel for: binary select tournament, crossover, mutate
			//and evaluate fitness. Maybe also have a chance to just copy
			//and mutate without crossover, too.
			cout<<"finished generation "<<gen<<endl;
		}




	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		cerr <<"Run with --help for help."<<endl;
	}

	return 0;
}