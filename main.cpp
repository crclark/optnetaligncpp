#include <iostream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
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

		//initialize population
		//todo: strongly consider dynamic allocation and avoiding copies
		cout<<"creating initial population"<<endl;
		unsigned int popsize = vm["popsize"].as<int>();
		vector<Alignment> pop;
		for(int i = 0; i < popsize; i++){	
			Alignment aln(*net1,*net2);
			aln.computeFitness(*net1,*net2,bitscores,evalues,fitnessNames);
			pop.push_back(aln);
		}

		cout<<"creating initial children"<<endl;
		vector<Alignment> kids;
		for(int i = 0; i < popsize; i++){
			Alignment aln(*net1,*net2);
			aln.computeFitness(*net1,*net2,bitscores,evalues,fitnessNames);
			kids.push_back(aln);
		}

		//main loop
		cout<<"starting main loop"<<endl;
		int generations = vm["generations"].as<int>();
		for(int gen = 0; gen < generations; gen++){

			vector<Alignment*> combinedPtrs(2*popsize);
			for(int i = 0; i < popsize; i++){
				combinedPtrs[i] = &pop[i];
			}
			for(int i = popsize; i < 2*popsize; i++){
				combinedPtrs.at(i) = &kids.at(i-popsize);
			}
			
			vector<vector<Alignment*> > fronts = nonDominatedSort(combinedPtrs);

			vector<Alignment*> popNew;
			popNew.reserve(popsize);
			int i = 0;
			while(popNew.size() + fronts[i].size() < popsize){
				setCrowdingDists(fronts[i]);
				popNew.insert(popNew.end(), fronts[i].begin(), fronts[i].end());
				i++;
			}

			
			int numLeftToInsert = popsize - popNew.size();
			if(!numLeftToInsert){
				sort(fronts[i].begin(),fronts[i].end(),crowdedComp);
				popNew.insert(popNew.end(), fronts[i].begin(), 
				          fronts[i].begin() + numLeftToInsert);
			}

		}




	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		cerr <<"Run with --help for help."<<endl;
	}

	return 0;
}