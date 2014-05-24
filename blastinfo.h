#pragma once

#include<unordered_map>
#include<tuple>
#include<string>
using namespace std;

#include "Network.h"

typedef unordered_map<node, unordered_map<node,double> > BLASTDict;

BLASTDict loadBLASTInfo(Network* net1, Network* net2, string filename);