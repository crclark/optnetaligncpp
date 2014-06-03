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
		const bool finalstats = vm.count("finalstats");
		const string outprefix = vm["outprefix"].as<string>();

		const BLASTDict* bitPtr = vm.count("bitscores") ? &bitscores : nullptr;

		mt19937 g(14);
		//initialize population
		if(verbose){
			cout<<"creating initial population"<<endl;
		}
		const unsigned int popsize = vm.count("popsize") 
		                             ? vm["popsize"].as<int>()
		                             : 100;
		vector<Alignment*> pop;
		for(int i = 0; i < popsize; i++){	
			Alignment* aln = new Alignment(net1,net2, bitPtr);
			aln->shuf(g, uniformsize, total);
			aln->computeFitness(bitscores,evalues,fitnessNames);
			pop.push_back(aln);
		}
		if(verbose){
			cout<<"creating initial children"<<endl;
		}
		vector<Alignment*> kids;
		for(int i = 0; i < popsize; i++){
			Alignment* aln = new Alignment(net1,net2, bitPtr);
			aln->shuf(g, uniformsize, total);
			aln->computeFitness(bitscores,evalues,fitnessNames);
			kids.push_back(aln);
		}

		//main loop
		if(verbose){
			cout<<"starting main loop"<<endl;
		}
		const int generations = vm["generations"].as<int>();
		for(int gen = 0; gen < generations; gen++){

			//combinedPtrs is R_t from Deb et al. 2002

			//todo: do this as a more idiomatic concatenation
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
				setCrowdingDists(fronts[i]);
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

			//check that pop has sane contents
			assert(pop.size() == popsize);
			for(int i = 0; i < pop.size(); i++){
				assert(popNew.count(pop[i]));
			}

			//do multithreaded version of kids creation
			//first resize kids to the proper size
			kids = vector<Alignment*>(popsize);

			//this will be executed by each thread.
			auto worker = [&](vector<Alignment*>::iterator begin,
				              vector<Alignment*>::iterator end,
				              int seed){
				mt19937 tg(seed + clock());
				for(auto it = begin; it != end; ++it){
					uniform_real_distribution<double> dist(0.0,1.0);
					double prob = dist(tg);
					if(prob <= 0.7){
						vector<Alignment*> parents;
						if(tournsel){
							parents = binSel(tg,pop,(popsize/10));
						}
						else{
							uniform_int_distribution<int> rint(0,popsize-1);
							int par1 = rint(tg);
							int par2 = par1;
							while(par2 == par1){
								par2 = rint(tg);
							}
							parents.push_back(pop[par1]);
							parents.push_back(pop[par2]);
						}
						*it = new Alignment(tg,cxswappb,*parents[0],
							               *parents[1], total);
						if(prob > 0.2){
							(*it)->mutate(tg,mutswappb,total);
						}
					}
					else{
						vector<Alignment*> parents = binSel(tg,pop,(popsize/10));
						(*it) = new Alignment(*parents[0]);
						(*it)->mutate(tg,mutswappb,total);
					}
					(*it)->computeFitness(bitscores,evalues,fitnessNames);
				}

			};

			const int grainsize = popsize / nthreads;
			if(nthreads == 1){
				worker(kids.begin(), kids.end(),g());
			}
			else{
				vector<boost::thread> threads(nthreads);
				auto work_iter = kids.begin();
				for(int i = 0; i < threads.size() - 1; i++){
					threads[i] = boost::thread(worker, work_iter, work_iter + grainsize,
						                i);
				    work_iter += grainsize;									
				}
				threads.back() = boost::thread(worker, work_iter, kids.end(), nthreads);

				for(auto&& i : threads){
					i.join();
				}
			}

			//idea: create a constructor for a non-random arbitrary aln
			//allocate kids using that constructor. Then,
			//in parallel for: binary select tournament, crossover, mutate
			//and evaluate fitness. Maybe also have a chance to just copy
			//and mutate without crossover, too.

			
			if(verbose){
				cout<<"Finished generation "<<gen<<endl;
				reportStats(pop,true);
			}
		}

		if(verbose){
			cout<<"Finished!"<<endl;
			cout<<"Writing alignments in Pareto front"<<endl;
		}
		if(finalstats){
			cout<<popsize;
			cout<<'\t'<<generations;
			cout<<'\t'<<mutswappb;
			cout<<'\t'<<cxswappb;
			cout<<'\t'<<tournsel;
			cout<<'\t'<<uniformsize;
			reportStats(pop,false);
		}
		vector <Alignment*> allAlns;
		allAlns.reserve(popsize*2);
		allAlns.insert(allAlns.end(), pop.begin(), pop.end());
		allAlns.insert(allAlns.end(), kids.begin(), kids.end());
		vector<vector<Alignment* > > fronts = nonDominatedSort(allAlns);

		string infoFilename = outprefix + ".info";
		ofstream infoFile(infoFilename);

		//make infoFile column labels
		infoFile << "filename";

		for(auto str : fitnessNames){
			infoFile << '\t' << str;
		}

		infoFile << endl;

		//output all alignments in the first front
		for(int i = 0; i < fronts[0].size(); i++){
			string filename = outprefix + "_" + to_string(i) + ".aln";
			fronts[0][i]->save(filename);

			infoFile << filename;

			//write summary info to infoFile
			for(int j = 0; j < fronts[0][i]->fitness.size(); j++){
				infoFile << '\t' << fronts[0][i]->fitness[j];
			}

			infoFile << endl;
		}


	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		cerr <<"Run with --help for help."<<endl;
	}

	return 0;
}