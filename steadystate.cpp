
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
		const string outprefix = vm["outprefix"].as<string>();

		const BLASTDict* bitPtr = vm.count("bitscores") ? &bitscores : nullptr;
		const GOCDict* gocPtr = vm.count("annotations1") ? &gocs : nullptr;
		const int generations = vm["generations"].as<int>();
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
		
		ArchiveMutexType archiveMutex;
		
		Archive archive;


		int numThreads = thread::hardware_concurrency();

		if(!numThreads){
			if(verbose){
				cout<<"Warning: failed to detect number of cores. "
					  "Falling back to single-threaded mode."<<endl;
			}
			numThreads++;
		}

		vector<int> randSeeds;

		for(int i = 0; i < numThreads; i++){
			randSeeds.push_back(g());
		}

		auto worker = [&](int threadNum){
			RandGenT tg(randSeeds[threadNum] + clock());
			while(numAlnsGenerated < popsize*generations){


				numAlnsGenerated++;
			}
		};

	
	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		cerr <<"Run with --help for help."<<endl;
	}



	return 0;
}
