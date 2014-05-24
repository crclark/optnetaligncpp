#include <iostream>
#include <string>
#include <tuple>
#include <boost/program_options.hpp>

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
		net1->edges.size();

	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		cerr <<"Run with --help for help."<<endl;
	}

	return 0;
}