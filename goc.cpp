#include <unordered_map>
#include <set>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <assert.h>
using namespace std;

#include "goc.h"
#include "Network.h"
GOCDict loadGOC(Network* net1, Network* net2,
                        string file1, string file2){

    auto anns1 = loadAnnotations(net1,file1);
    auto anns2 = loadAnnotations(net2,file2);
    
    GOCDict toReturn;
    for(auto pair1 : anns1){
        for(auto pair2 : anns2){
            set<int> un;
            set_union(pair1.second.begin(),pair1.second.end(),
                      pair2.second.begin(),pair2.second.end(),
                      inserter(un, un.begin()));
                      
            set<int> in;
            set_intersection(pair1.second.begin(),pair1.second.end(),
                             pair2.second.begin(),pair2.second.end(),
                             inserter(in,in.begin()));
                             
            double toStore = double(in.size())/double(un.size());
            toReturn[pair1.first][pair2.first] = toStore;
        }
    }
    return toReturn;
}


unordered_map<node, set<int> > loadAnnotations(Network* net, string fp){
    
    ifstream infile(fp);
	if(!infile.good()){
		throw LineReadException("Failed to load annotations! Check filename.");
	}
	unordered_map<node,set<int> > toReturn;

	string line;
	unsigned int count = 0;

	while(getline(infile,line)){
		istringstream iss(line);

        string a;
        
		node u;
        
        iss >> a;
        
        set<int> anns;
        
        int x; 
        
        while(iss >> x){
            anns.insert(x);
        }
        
        
		u = net->nodeNameToNode[a];
		toReturn[u] = anns;
		count++;
	}
	return toReturn;
    
}
