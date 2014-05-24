#include <iostream>
#include <string>
#include <boost/program_options.hpp>

#include "argParsing.h"
using namespace std;
namespace po = boost::program_options;



int main(int ac, char* av[])
{
	try{
		
		po::variables_map vm = handleArgs(ac, av);

	}
	catch(exception& e){
		cerr << "error: " << e.what() << endl;
		cerr <<"Run with --help for help."<<endl;
	}

	return 0;
}