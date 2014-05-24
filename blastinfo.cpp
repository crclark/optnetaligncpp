#include <unordered_map>
#include <tuple>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "Network.h"
#include "blastinfo.h"

BLASTDict loadBLASTInfo(Network* net1, Network* net2, string filename){

	ifstream infile(filename);

	BLASTDict toReturn;

	string line;
	unsigned int count = 0;

	while(getline(infile,line)){
		istringstream iss(line);

		double blastScore;
		string a, b;
		node u, v;

		if (!(iss >> a >> b >> blastScore)){
			throw LineReadException(string("Parse error in network: ") + filename +
				                   string("on line: ") + line + "\n");
		}
		u = net1->nodeNameToNode[a];
		v = net2->nodeNameToNode[b];
		toReturn[u][v] = blastScore;
		count++;
	}
	cout<<count<<" lines of blast info parsed from "<<filename<<endl;
	return toReturn;
}