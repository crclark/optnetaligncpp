
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
#include <thread>

#include "Alignment.h"
#include "argParsing.h"
#include "nsga-ii.h"
#include "localSearch.h"
#include "Archive.h"

#include "tbb/spin_mutex.h"
#include "tbb/atomic.h"
using namespace std;
using namespace tbb;
namespace po = boost::program_options;

typedef tbb::spin_mutex ArchiveMutexType;

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
		const float cxrate = vm.count("cxrate") ? vm["cxrate"].as<float>() : 0.7;
		const bool verbose = vm.count("verbose");
		const bool tournsel = vm.count("tournsel");
		const bool total = vm.count("total");
		const bool uniformsize = vm.count("uniformsize");
		const bool smallstart = vm.count("smallstart");
		const bool finalstats = vm.count("finalstats");
        const bool nooutput = vm.count("nooutput");
		const string outprefix = vm["outprefix"].as<string>();

		const BLASTDict* bitPtr = vm.count("bitscores") ? &bitscores : nullptr;
		const GOCDict* gocPtr = vm.count("annotations1") ? &gocs : nullptr;
		const int generations = vm["generations"].as<int>();
		const int hillclimbiters = vm.count("hillclimbiters") 
		                           ? vm["hillclimbiters"].as<int>() 
		                           : 0; 
		const unsigned int popsize = vm.count("popsize") 
		                             ? vm["popsize"].as<int>()
		                             : 100;
		const int randseed = vm.count("randseed") ? vm["randseed"].as<int>() : clock();
			
		RandGenT g(randseed);

		tbb::atomic<int> numAlnsGenerated; //this will be our stopping condition
		                                  //also use this as condition to shrink archive
										  //todo: perhaps make this count only
										  //non-dominated alns generated, instead of all.
		
		numAlnsGenerated.store(0);
		
		tbb::atomic<int> numPropSearch;
		numPropSearch.store(0);

		tbb::atomic<int> numCx;
		numCx.store(0);

		tbb::atomic<int> numNonDominatedGenerated;
		numNonDominatedGenerated.store(0);

		tbb::atomic<int> nonDomPropSearch;
		nonDomPropSearch.store(0);

		tbb::atomic<int> nonDomCx;
		nonDomCx.store(0);

		ArchiveMutexType archiveMutex;
		Archive archive;


		int numThreads = vm.count("nthreads")
                         ? vm["nthreads"].as<int>()
                         : thread::hardware_concurrency();

		if(!numThreads){
			if(verbose){
				cout<<"Warning: failed to detect number of cores. "
					  "Falling back to single-threaded mode."<<endl;
			}
			numThreads++;
		}

		if(verbose){
			cout<<"Launching "<<numThreads<<" thread(s)."<<endl;
		}
		vector<int> randSeeds;
		for(int i = 0; i < numThreads; i++){
			randSeeds.push_back(g());
		}

		auto worker = [&](int threadNum){
			RandGenT tg(randSeeds[threadNum] + clock());
			uniform_real_distribution<double> prob(0.0,1.0);
			uniform_int_distribution<int> randObj(0, fitnessNames.size()-1);
			while(numAlnsGenerated < popsize*generations){

				//grab 2 existing alns from archive.
				//if less than 2 in archive, just create a new one.
				//todo: these constructors waste time except for the first iter
				Alignment par1(net1, net2, bitPtr, gocPtr);
				Alignment par2(net1, net2, bitPtr, gocPtr);

				bool foundAln = true;
				bool didCx = false;
				bool didPropSearch = false;
				//lock and get from archive, aborting if archive too small
				{
					ArchiveMutexType::scoped_lock lock(archiveMutex);
					int archsize = archive.nonDominated.size();
					if(archsize == 0){
						foundAln = false;
					}
					else{
						uniform_int_distribution<int> r(0, archsize-1);
						int i1 = r(tg);
						int i2 = r(tg);

						par1 = Alignment(*archive.nonDominated.at(i1));
						par2 = Alignment(*archive.nonDominated.at(i2));
					}
				}

				//if we didn't get 2 alns from archive, init with random ones
				if(!foundAln){
					par1.shuf(tg, uniformsize, smallstart, total);
					par2.shuf(tg, uniformsize, smallstart, total);
				}

				Alignment* child;

				//do crossover or mutation
				if(prob(tg) < cxrate){
					child = new Alignment(tg, cxswappb, par1, par2, total);
					didCx = true;
				}
				else{
					child = new Alignment(par1);
					child->mutate(tg, mutswappb, total);
				}

				//now, the local search...
				//Should we do proportional or standard?
				//And what should our stopping criterion be?
				VelocityTracker veltracker;
				child->computeFitness(fitnessNames);
				vector<double> initSpeed;
				//for proportional search, decide which way to search
				int rObj = randObj(tg);

				//now prefering propsearch when seeding. See if that helps.
				if(!foundAln || prob(tg) < 0.01){ //todo: make this probability a param
					didPropSearch = true;
				}

				for(int i = 0; i < hillclimbiters; i++){
					vector<double> currFit = child->fitness;
					
					if(didPropSearch){
						//todo: used to decide whether to do proportional
						//or standard hillclimb on each hillclimbiter.
						//Decide whether to keep that old behavior or
						//switch to this new behavior.
						proportionalSearch(tg, child, total,
	                    	1, fitnessNames,
	                    	rObj, 0.95);
					}
					else{
					//note: 1 iter of proportionalSearch = 500 iters of hillclimb
					correctHillClimb(tg, child, total, 500, fitnessNames);
					}
					vector<double> newFit = child->fitness;

					vector<double> delta(newFit.size(), 0.0);


					for(int j = 0; j < delta.size(); j++){
						delta[j] = newFit[j] - currFit[j];
					}
					veltracker.reportDelta(delta);

					if(i == 500){
						initSpeed = veltracker.getRecentVel();
					}

					if(i > 5000 && i % 500 == 0){
						bool belowThresh = true;
						vector<double> currAvgVel = veltracker.getRecentVel();
						for(int j = 0; j < currAvgVel.size(); j++){
							belowThresh &= currAvgVel[j] 
                                            < 0.01*initSpeed[j];
						}
                        
                        if(belowThresh){
                            break;
                        }
					}
				}
                
				//since proportional search can leave us somewhere suboptimal,
				//do correctHillClimb thereafter
				correctHillClimb(tg, child, total, 500*hillclimbiters,
								 fitnessNames);

                //insert in archive
				{
                    ArchiveMutexType::scoped_lock lock(archiveMutex);
                    bool wasNonDominated = archive.insert(child);
                    if(wasNonDominated){
                    	numNonDominatedGenerated++;
                    	if(didPropSearch){
                    		nonDomPropSearch++;
                    	}
                    	if(didCx){
                    		nonDomCx++;
                    	}
                    }
                    if(archive.nonDominated.size() > popsize){
                        archive.shrinkToSize(popsize);
                    }
                }
                

                if(numAlnsGenerated % 100 == 0 && verbose){
                    ArchiveMutexType::scoped_lock lock(archiveMutex);
                    reportStats(archive.nonDominated, fitnessNames,
                                verbose);
                    cout<<archive.nonDominated.size()<<" non-dominated."<<endl;
                    cout<<numAlnsGenerated<<" created total."<<endl;
                    cout<<numPropSearch<<" created with proportionalSearch."<<endl;
                    cout<<numCx<<" created with crossover."<<endl;
                    cout<<numNonDominatedGenerated<<" of created were non-dominated."<<endl;
                    cout<<nonDomCx<<" of non-dominated were made with crossover."<<endl;
                    cout<<nonDomPropSearch<<" of non-dominated were made with proportionalSearch."<<endl;
                }
                
				numAlnsGenerated++;
				if(didPropSearch){
					numPropSearch++;
				}
				if(didCx){
					numCx++;
				}
			}
		};

		vector<thread> ts;

		for(int i = 0; i < numThreads; i++){
			ts.push_back(thread(worker,i));
		}

		for(int i = 0; i < numThreads; i++){
			ts[i].join();
		}
        
        if(verbose){
            cout<<"Found "<<archive.nonDominated.size()
                <<" non-dominated alignments."<<endl;
            cout<<"Writing to disk!"<<endl;
        }
        
        if(finalstats){
			cout<<popsize;
			cout<<'\t'<<generations;
			cout<<'\t'<<mutswappb;
			cout<<'\t'<<cxswappb;
			cout<<'\t'<<tournsel;
			cout<<'\t'<<uniformsize;
			cout<<'\t'<<cxrate;
			cout<<'\t'<<hillclimbiters;
			cout<<'\t'<<(double(numNonDominatedGenerated)/double(numAlnsGenerated));
			cout<<'\t'<<double(nonDomPropSearch)/double(numNonDominatedGenerated);
			cout<<'\t'<<double(nonDomCx)/double(numNonDominatedGenerated);
			reportStats(archive.nonDominated,fitnessNames,false);
			cout<<endl;
		}
		
        if(!nooutput){
            string infoFilename = outprefix + ".info";
            ofstream infoFile(infoFilename);

            //make infoFile column labels
            infoFile << "filename";

            for(auto str : fitnessNames){
                infoFile << '\t' << str;
            }

            infoFile << endl;

            //output all alignments in the first front
            for(int i = 0; i < archive.nonDominated.size(); i++){
                string filename = outprefix + "_" + to_string(i) + ".aln";
                archive.nonDominated[i]->save(filename);

                infoFile << filename;

                //write summary info to infoFile
                for(int j = 0; j < archive.nonDominated[i]->fitness.size(); j++){
                    infoFile << '\t' << archive.nonDominated[i]->fitness[j];
                }

                infoFile << endl;
            }
        }
        
	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		cerr <<"Run with --help for help."<<endl;
	}



	return 0;
}
