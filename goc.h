#pragma once

#include <unordered_map>
#include <set>
#include <tuple>
#include <string>
using namespace std;

#include "Network.h"

typedef unordered_map<node, unordered_map<node,double> > GOCDict;

GOCDict loadGOC(Network* net1, Network* net2,
                        string file1, string file2);

unordered_map<node, set<int> > loadAnnotations(Network* net, string fp);
