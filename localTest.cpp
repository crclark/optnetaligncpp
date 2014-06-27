#include <iostream>
#include <fstream>
#include <string>
#include <tuple>
#include <array>
#include <algorithm>
#include <unordered_set>
#include <ctime>
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
		GOCDict gocs = get<5>(vals);
		vector<string> fitnessNames = get<6>(vals);

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
		const GOCDict* gocPtr = vm.count("annotations1") ? &gocs : nullptr;
		const int generations = vm["generations"].as<int>();

		mt19937 g(14);
		//initialize population
		
		
		Alignment* aln = new Alignment(net1,net2, &bitscores, &gocs);
		//aln->greedyMatch(false);
		aln->shuf(g,false,false,total);
		aln->computeFitness(fitnessNames);
		
		//todo: instead of just flipping obj, switch according to some
		//input time proportion.
		//fast hill climb version
		cout<<"starting main loop"<<endl;
		VelocityTracker velTracker;
		int gensWithoutImprovement = 0;
		int mutsDone = 0;
		vector<double> initSpeed;
		vector<double> bestFits(fitnessNames.size(),0.0);
		for(int i = 0; i < generations; i++){
			auto oldFit = aln->fitness;
			correctHillClimb(g, aln, total,
               500, fitnessNames);	
			auto newFit = aln->fitness;
			for(int q = 0; q < newFit.size(); q++){
				if(newFit[q] > bestFits[q]){
					bestFits[q] = newFit[q];
				}
			}
			
			vector<double> deltaFit;
			for(int j = 0; j < newFit.size(); j++){
				deltaFit.push_back(newFit[j] - oldFit[j]);
			}
			velTracker.reportDelta(deltaFit);

			if(!dominates(newFit,oldFit)){
				gensWithoutImprovement++;
			}
			else{
				gensWithoutImprovement = 0;
			}
			if(i > 10000 && gensWithoutImprovement == (i/100)){
				aln->mutate(g, 0.001, total);
				aln->computeFitness(fitnessNames);
				gensWithoutImprovement = 0;
				mutsDone++;
			}
			for(int j = 0; j <fitnessNames.size();j++){
				cout<<"current "<<fitnessNames.at(j)<<" is "
				    <<aln->fitness.at(j)<<endl;
				cout<<"best "<<fitnessNames.at(j)<<" is "
				    <<bestFits.at(j)<<endl;
			}
            cout<<"EC is "<<((double)(aln->currConservedCount))/((double)(net1->edges.size()))<<endl;
            cout<<"Velocity is ";
            vector<double> currVel = velTracker.getRecentVel();
            for(int j = 0; j < currVel.size(); j++){
            	cout<<currVel[j]<<' ';
            }
            cout<<endl;
            cout<<"Velocity-based optimum detected: ";
            bool belowThresh = true;
            if(i > 50){			
				for(int j = 0; j < currVel.size(); j++){
					belowThresh &= currVel[j] 
	                                < 0.001*initSpeed[j];
				}
			}
			else{
				belowThresh = false;
			}

			if(i == 50){
				initSpeed = velTracker.getRecentVel();
			}
			if(belowThresh){
				cout<<"yes."<<endl;
			}
			else{
				cout<<"no."<<endl;
			}
            cout<<"Mutations performed: "<<mutsDone<<endl;
			cout<<"Generation "<<i<<" complete."<<endl;
		}
		aln->save(outprefix + "_localTest.aln");
		
	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		cerr <<"Run with --help for help."<<endl;
	}



	return 0;
}
