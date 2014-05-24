#pragma once

#include <exception>
#include <string>
#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

class ArgError : public exception{
public:

	ArgError(string msg){
		errMsg = msg;
	}

	virtual const char* what() const throw(){
		return errMsg.c_str();
	}

private:
	string errMsg;
};

//todo: this is a terrible crime, but I am going to use globals for 
//the networks, bit scores, E-values, and annotations. Put them here.

po::variables_map handleArgs(int ac, char* av[]){
	po::options_description desc("Options");
	//todo: add to help descriptions to indicate which args optional.
	desc.add_options()
		("help", "Displays this help message.")
		("int", po::value<int>(), "this is an int.")
		("net1", po::value<string>(), "Path to first network file. "
			                          "One interaction per line,"
			                          "represented as two protein names separated "
			                          "by whitespace. See documentation for an "
			                          "example. The network may at most have as "
			                          "many nodes as net2.")
		("net1", po::value<string>(), "Path to the second network file. "
			                          "Same format expected as first network "
			                          "file.")
		("outprefix", po::value<string>(), "Prefix for all output files.")
		("nthreads", po::value<int>(), "Number of threads to use. More will "
			                           "generally allow the program to run "
			                           "faster. If unspecified, defaults to "
			                           "the number of cores in your machine.")
		("total", "When set, restricts alignments to total alignment only.")
		("popsize", po::value<int>(), "The number of alignments to maintain "
			                          "in the genetic algorithm's population.")
		("generations", po::value<int>(), "The max number of generations to "
			                              "run the algorithm.")
		("initlist", po::value<string>(), "Path to a list of paths to existing "
			                              "alignments to use to seed the search. "
			                              "Optional.")
		("bitscores", po::value<string>(), "Path to a list of BLAST bit scores "
			                               "for the given networks. Each line "
			                               "should be the names of two nodes and "
			                               "their bitscore, whitespace-separated."
			                               "When this argument is used, the bit "
			                               "score information will be used in "
			                               "creating alignments.")
		("evalues", po::value<string>(), "Path to E-values for the given "
			                             "networks. Uses the same format as "
			                             "bit scores data.")
		("annotations1", po::value<string>(), "Path to a list of GO annotations "
			                                  "for network 1. Each line should "
			                                  "be a node name followed by its "
			                                  "GO annotations, separated by "
			                                  "whitespace. When a path to "
			                                  "annotations is specified for both "
			                                  "networks, this info is used as "
			                                  "an objective to guide the "
			                                  "alignment.")
		("annotations2", po::value<string>(), "Path to GO annotations for net 2.")
		("ics", "When set, integrated conserved structure score will be used as "
			    "an alignment objective.")
		("verbose", "When set, extra information about the progress of the "
			        "alignment is printed to stdout.")
		("mutswappb", po::value<float>(), "Sets the probability of swap "
			                              "mutation. See documentation.")
		("cxswappb", po::value<float>(), "Sets the probability of swapping "
			                             "during crossover. See documentation.")
		("finalstats", "When set, prints comma-separated stats to stdout "
			           "when execution completes. Used to find good parameter "
			           "values through experimentation.")
		("seeding", "When set, attempts to start with better-than-random "
			        "alignments by creating randomized seed-and-extend "
			        "alignments to serve as the initial population. "
			        "Must be used in conjunction with bit scores or E-"
			        "values.")
		("tournsel", "When set, use domination-based tournament selection "
			         "to choose parents for crossover. Otherwise, uses "
			         "random selection.")
	;
	
	po::variables_map vm;
	po::store(po::parse_command_line(ac,av,desc), vm);
	po::notify(vm);

	if(vm.count("help")){
		cout<<desc<<endl;
		throw ArgError("Please run with valid arguments.");
	}


	if(!(vm.count("ics") || vm.count("bitscores") || vm.count("evalues")
		|| (vm.count("annotations1") && vm.count("annotations2")) )){
		throw ArgError("At least one objective must be specified!")
	}

	if((vm.count("annotations1") && !vm.count("annotations2")) ||
	   (vm.count("annotations2") && !vm.count("annotations1"))){
		throw ArgError("Please specify annotations for both networks "
			           "or neither.")
	}

	if(vm.count("seeding") && !(vm.count("bitscores") || vm.count("evalues"))){
		throw ArgError("Seeding requires either bitscore or E-value data.")
	}

	
	if(!vm.count("int")){
		throw ArgError("testtesttest");
	}

	if(vm.count("net1")){
		throw ArgError(vm["net1"].as<string>());
	}

	return vm;
}