#pragma once

#include <unordered_set>
#include <string>
#include <unordered_map>
using namespace std;

typedef int node;

//for giving parse errors
class LineReadException: public exception
{

public:
	LineReadException(string msg){
		wt = msg;
	}

	virtual const char* what() const throw()
	{
	return wt.c_str();
	}

	string wt;
};

//Edges are undirected, so they are always normalized such that
//the first node is a smaller int than the second.
class Edge{
public:
	Edge(node u, node v){
		if(u <= v){
			n1 = u;
			n2 = v;
		}
		else{
			n1 = v;
			n2 = u;
		}
	}

	node u(){
		return n1;
	}

	node v(){
		return n2;
	}

	bool operator==(const Edge& other) const{
		return n1 == other.n1 && n2 == other.n2;
	}

	bool operator!=(const Edge& other) const{
		return !(*this == other);
	}

	bool operator<(const Edge& other) const{
		/*if(n1 < other.n1){
			return true;
		}
		else if(n1 == other.n1){
			return n2 < other.n2;
		}
		else{
			return false;
		}*/
		return (n1 < other.n1) || (n1 == other.n1 && n2 < other.n2);
	}

	friend class EdgeHash;

private:
	node n1;
	node n2;
};

class EdgeHash{
public:
	size_t operator()(const Edge& e) const {
		return e.n1 ^ e.n2;
	}
};

//A network is, for now, just a set of edges. If we
//end up needing real network functionality, I'll add a
//real network library.

class Network{
public:
	Network(string filename);
	//set<Edge> edges;
	unordered_set<Edge, EdgeHash> edges;
	unordered_map<node,string> nodeToNodeName;
	unordered_map<string,node> nodeNameToNode;
};