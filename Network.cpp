#include "Network.h"

#include <set>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;


Network::Network(string filename){
	unsigned int count = 0;
	set<string> alreadySeen;

	ifstream infile(filename);

	string line;
	while(getline(infile,line)){
		istringstream iss(line);
		string a, b;
		node u, v;

		if (!(iss >> a >> b)){
			throw LineReadException(string("Parse error in network: ") + filename +
				                   string("on line: ") + line + "\n");
		}

		if(!alreadySeen.count(a)){
			alreadySeen.insert(a);
			count++;
			nodeToNodeName[count] = a;
			nodeNameToNode[a] = count;
			u = count;
		}
		else{
			u = nodeNameToNode[a];
		}

		if(!alreadySeen.count(b)){
			alreadySeen.insert(b);
			count++;
			nodeToNodeName[count] = b;
			nodeNameToNode[b] = count;
			v = count;
		}
		else{
			v = nodeNameToNode[b];
		}

		Edge e = Edge(u,v);

		edges.insert(e);

	}

	cout<<"Network "<<filename<<" has "<<count<<" nodes and "
	    <<edges.size()<<" edges"<<endl;
}