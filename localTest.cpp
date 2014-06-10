#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <unordered_set>
#include <ctime>
#include <boost/thread/thread.hpp>
#include <assert.h>
#include <boost/program_options.hpp>

#include "Alignment.h"
#include "argParsing.h"
#include "nsga-ii.h"
#include "localSearch.h"
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

		const int nthreads = vm.count("nthreads") ? vm["nthreads"].as<int>()
		                                          : 1;
		const float mutswappb = vm.count("mutswappb")  
		                             ? vm["mutswappb"].as<float>()
		                             : 0.005;
		const float cxswappb = vm.count("cxswappb") ? vm["cxswappb"].as<float>()
		                                            : 0.1;
		const bool verbose = vm.count("verbose");
		const bool tournsel = vm.count("tournsel");
		const bool total = vm.count("total");
		const bool uniformsize = vm.count("uniformsize");
		const bool smallstart = vm.count("smallstart");
		const bool finalstats = vm.count("finalstats");
		const string outprefix = vm["outprefix"].as<string>();

		const BLASTDict* bitPtr = vm.count("bitscores") ? &bitscores : nullptr;

		const int generations = vm["generations"].as<int>();

		mt19937 g(14);
		//initialize population
		
		/*
		Alignment* aln = new Alignment(net1,net2, &bitscores);
		//aln->greedyBitscoreMatch();
		aln->shuf(g,false,false,total);
		aln->computeFitness(fitnessNames);
		
		//todo: instead of just flipping obj, switch according to some
		//input time proportion.
		//fast hill climb version
		cout<<"starting main loop"<<endl;
		int obj = 1;
		for(int i = 0; i < generations; i++){
			fastHillClimb(g, aln, total,
	               500, fitnessNames, obj,true);
			for(int j = 0; j <fitnessNames.size();j++){
				cout<<"current "<<fitnessNames.at(j)<<" is "
				    <<aln->fitness.at(j)<<endl;
			}
			//obj = (obj + 1) % 3;
			cout<<"Generation "<<i<<" complete."<<endl;
		}
		aln->save("localTest.aln");
		*/

		vector<Alignment*> pop(10, nullptr);
		double prop = 0.0;
		for(int i = 0; i < 10; i++){
			prop = 1.0 / (double(i));
			pop[i] = new Alignment(net1,net2,&bitscores);
			pop[i]->shuf(g,false,false,total);
			pop[i]->computeFitness(fitnessNames);
			cout<<"calling proportionalSearch"<<endl;
			proportionalSearch(g, pop[i], total,
	                    10000, fitnessNames,
	                    prop);
			cout<<"First alignment: "<<endl;
			for(int j = 0; j <fitnessNames.size();j++){
				cout<<fitnessNames.at(j)<<" is "
				    <<pop[i]->fitness.at(j)<<endl;
			}
			cout<<endl;
		}
	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		cerr <<"Run with --help for help."<<endl;
	}



	return 0;
}