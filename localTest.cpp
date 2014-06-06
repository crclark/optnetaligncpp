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

		mt19937 g(14);
		//initialize population
		
		Alignment* aln = new Alignment(net1,net2, &bitscores);
		aln->shuf(g,false,false,total);
		aln->computeFitness(fitnessNames);

		//slow hill climb version
		/*
		cout<<"starting main loop"<<endl;
		for(int i = 0; i < 10000; i++){
			cout<<"calling hillClimb"<<endl;
			Alignment* temp = hillClimb(g, aln, total,
	                  5, fitnessNames,0);
			if(temp != aln)
				delete aln;
			aln = temp;
			cout<<"curr ICS: "<<aln->fitness.at(0)<<endl;
		}
		*/
		
		//fast hill climb version
		cout<<"starting main loop"<<endl;
		int obj = 0;
		for(int i = 0; i < 100000; i++){
			fastHillClimb(g, aln, total,
	               500, fitnessNames, obj);
			for(int j = 0; j <fitnessNames.size();j++){
				cout<<"current "<<fitnessNames.at(j)<<" is "
				    <<aln->fitness.at(j)<<endl;
			}
			obj = obj ? 0 : 1;
		}
		aln->save("localTest.aln");

	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		cerr <<"Run with --help for help."<<endl;
	}

	return 0;
}