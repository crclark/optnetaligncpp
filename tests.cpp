//  (C) Copyright Gennadiy Rozental 2005.
//  Distributed under the Boost Software License, Version 1.0.
//  (See accompanying file LICENSE_1_0.txt or copy at 
//  http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/test for the library home page.

// Boost.Test

// each test module could contain no more then one 'main' file with init function defined
// alternatively you could define init function yourself
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <utility>
#include <cmath>
#include <iostream>
#include <set>
#include <random>
#include <limits>
#include "Alignment.h"
#include "Network.h"
#include "blastinfo.h"
#include "nsga-ii.h"
using namespace std;

bool approxEqual(double x, double y){
	return abs(x-y) < 0.000000001;
}

//____________________________________________________________________________//

// most frequently you implement test cases as a free functions with automatic registration
BOOST_AUTO_TEST_CASE( bitscores_correct_length )
{
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	BLASTDict d = loadBLASTInfo(&net1,&net2,"../optnetalign/tests/cg1.sim");
	int size = 0;
	for(auto subd : d){
		size += std::get<1>(subd).size();
	}
    BOOST_CHECK( size == 18181 );
}

//____________________________________________________________________________//

// each test file may contain any number of test cases; each test case has to have unique name

BOOST_AUTO_TEST_CASE( bitscoreSum_match_1 )
{
    Network net1("../optnetalign/tests/small.net");
    Network net2("../optnetalign/tests/small.net");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/small.bitScores");
    Alignment a= Alignment(&net1,&net2,"../optnetalign/tests/small.aln",&b);
    double sumbit = a.sumBLAST();
    //std::cout<<"sumbit is "<<sumbit<<std::endl;
    BOOST_CHECK(approxEqual(2.0,a.sumBLAST()));
    // reports 'error in "test2": check i == 2 failed [0 != 2]'
    //BOOST_CHECK_EQUAL( i, 2 );

    //BOOST_CHECK_EQUAL( i, 0 );
}

BOOST_AUTO_TEST_CASE( bitscoreSum_match_2 )
{
	Network net1("../optnetalign/tests/small.net");
    Network net2("../optnetalign/tests/small.net");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/small.bitScores");
    Alignment a2 = Alignment(&net1,&net2,"../optnetalign/tests/small2.aln",&b);
    double sumbit2 = a2.sumBLAST();
    std::cout<<"sumbit2 is "<<sumbit2<<std::endl;
    BOOST_CHECK(approxEqual(sumbit2,0.0));
}

BOOST_AUTO_TEST_CASE( currBitscore_consistent ){
	Network net1("../optnetalign/tests/small.net");
    Network net2("../optnetalign/tests/small.net");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/small.bitScores");
    Alignment a2 = Alignment(&net1,&net2,&b);
    mt19937 g1(12);
    a2.shuf(g1,false,false,true);
    double sumbit2 = a2.sumBLAST();
    BOOST_CHECK(approxEqual(sumbit2,a2.currBitscore));
}


BOOST_AUTO_TEST_CASE( bitscoreSum_match_3 )
{
	Network net1("../optnetalign/tests/bittest1.net");
    Network net2("../optnetalign/tests/bittest2.net");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/bittest.sim");
    Alignment a2 = Alignment(&net1,&net2,&b);
    double sumbit2 = a2.sumBLAST();
    std::cout<<"sumbit2 is "<<sumbit2<<std::endl;
    BOOST_CHECK(approxEqual(sumbit2,11.0+22.0+33.0+44.0));
}

BOOST_AUTO_TEST_CASE( currBitscore_consistent_2 ){
	Network net1("../optnetalign/tests/bittest1.net");
    Network net2("../optnetalign/tests/bittest2.net");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/bittest.sim");
    Alignment a2 = Alignment(&net1,&net2,&b);
    
    Alignment a3 = Alignment(&net1,&net2,&b);
    mt19937 g1(12);
    a2.shuf(g1,false,false,false);
    a3.shuf(g1,false,false,false);
    Alignment a4 = Alignment(g1,0.5,a2,a3,false);
    double sumbit2 = a4.sumBLAST();
    cout<<"sumbit2: "<<sumbit2<<endl;
    cout<<"currBitscore: "<<a4.currBitscore<<endl;
    BOOST_CHECK(approxEqual(sumbit2,a4.currBitscore));
}

BOOST_AUTO_TEST_CASE( currBitscore_consistent_3 ){
	Network net1("../optnetalign/tests/bittest1.net");
    Network net2("../optnetalign/tests/bittest2.net");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/bittest.sim");
    Alignment a2 = Alignment(&net1,&net2,&b);
    a2.currBitscore = a2.sumBLAST();
    mt19937 g1(12);
    a2.mutate(g1,0.25);
    double sumbit2 = a2.sumBLAST();
 	cout<<"sumbit2: "<<sumbit2<<endl;
    cout<<"currBitscore: "<<a2.currBitscore<<endl;   
    BOOST_CHECK(approxEqual(sumbit2,a2.currBitscore));
}

BOOST_AUTO_TEST_CASE( currBitscore_consistent_4 ){
	Network net1("../optnetalign/tests/bittest1.net");
    Network net2("../optnetalign/tests/bittest2.net");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/bittest.sim");
    Alignment a2 = Alignment(&net1,&net2,&b);
    a2.currBitscore = a2.sumBLAST();
    node temp = a2.aln[1];
    a2.aln[1] = a2.aln[2];
    a2.aln[2] = temp;

    a2.updateBitscore(1, temp, a2.aln[1], true, true);
    a2.updateBitscore(2, a2.aln[1], a2.aln[2], true, true);

    double sumbit2 = a2.sumBLAST();
    BOOST_CHECK(approxEqual(sumbit2,a2.currBitscore));
}

BOOST_AUTO_TEST_CASE( currBitscore_consistent_5 ){
	Network net1("../optnetalign/tests/cg1a.net");
    Network net2("../optnetalign/tests/cg1b.net");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/cg1.sim");
    Alignment a2 = Alignment(&net1,&net2,&b);
    
    Alignment a3 = Alignment(&net1,&net2,&b);
    mt19937 g1(12);
    a2.shuf(g1,false,false,false);
    a3.shuf(g1,false,false,false);
    Alignment a4 = Alignment(g1,0.5,a2,a3,false);
    double sumbit2 = a4.sumBLAST();
    cout<<"sumbit2: "<<sumbit2<<endl;
    cout<<"currBitscore: "<<a4.currBitscore<<endl;
    BOOST_CHECK(approxEqual(sumbit2,a4.currBitscore));
}

BOOST_AUTO_TEST_CASE( currBitscore_consistent_6 ){
	Network net1("../optnetalign/tests/cg1a.net");
    Network net2("../optnetalign/tests/cg1b.net");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/cg1.sim");
    Alignment a2 = Alignment(&net1,&net2,&b);
    
    Alignment a3 = Alignment(&net1,&net2,&b);
    mt19937 g1(1244);
    a2.shuf(g1,false,false,true);
    a3.shuf(g1,false,false,true);
    Alignment a4 = Alignment(g1,0.5,a2,a3,true);
    double sumbit2 = a4.sumBLAST();
    cout<<"sumbit2: "<<sumbit2<<endl;
    cout<<"currBitscore: "<<a4.currBitscore<<endl;
    BOOST_CHECK(approxEqual(sumbit2,a4.currBitscore));
}

BOOST_AUTO_TEST_CASE( bitscores_contains_correct_records ){
	Network net1("../optnetalign/tests/cg1a.net");
    Network net2("../optnetalign/tests/cg1b.net");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/cg1.sim");
    BOOST_CHECK(b.count(net1.nodeNameToNode.at("a1")));
}

BOOST_AUTO_TEST_CASE( save_load_inverses )
{
	Network net1("../optnetalign/tests/lccstest1.net");
	Network net2("../optnetalign/tests/lccstest2.net");
	Alignment aln1(&net1,&net2,"../optnetalign/tests/lccstest.aln",nullptr);
	aln1.save("../optnetalign/tests/testSavecpp.aln");
	Alignment aln2(&net1,&net2,"../optnetalign/tests/testSavecpp.aln",nullptr);

	for(int i = 0; i < net1.nodeToNodeName.size(); i++){
		BOOST_CHECK(aln1.aln[i] == aln2.aln[i]);
	}
}

BOOST_AUTO_TEST_CASE( rand_aln_consistent )
{
	Network net1("../optnetalign/tests/lccstest1.net");
	Network net2("../optnetalign/tests/lccstest2.net");
	Alignment aln1(&net1,&net2,nullptr);
	set<node> v2Set;
	//check that alignment is valid permutation.
	for(int i = 0; i < aln1.aln.size(); i++){
		v2Set.insert(aln1.aln[i]);
		BOOST_CHECK(aln1.aln[i] >= 0 && 
			        aln1.aln[i] < net2.nodeToNodeName.size());
	}
	BOOST_CHECK(v2Set.size() == aln1.aln.size());
}

BOOST_AUTO_TEST_CASE( load_aln_consistent )
{
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,"../optnetalign/tests/cg1.aln",nullptr);
	set<node> v2Set;
	vector<int> counts(aln1.aln.size(),0);

	for(int i = 0; i < aln1.aln.size(); i++){
		BOOST_CHECK(aln1.aln[i] >= 0 &&
			        aln1.aln[i] < net2.nodeToNodeName.size());
		if(v2Set.count(aln1.aln.at(i))){
			cout<<"Already aligned! aln1.aln["<<i<<"]: "<<net2.nodeToNodeName[aln1.aln.at(i)]<<endl;
		}
		v2Set.insert(aln1.aln.at(i));
		if(counts.at(aln1.aln.at(i)) > 1){
			cout<<"node "<<aln1.aln.at(i)<<" appears twice!"<<endl;
		}
		counts[aln1.aln.at(i)]++;
	}

	for(int i = 0; i < aln1.aln.size(); i++){
		for(int j=i+1; j < aln1.aln.size(); j++){
			if(i == j){
				cout<<"found same node aligned to twice."<<endl;
				cout<<"node "<<net2.nodeToNodeName[i]<<endl;
			}
		}
	}

	for(int i = 0; i < aln1.aln.size(); i++){
		if(!v2Set.count(aln1.aln[i])){
			cout<<"missing node "<<net2.nodeToNodeName[aln1.aln[i]]<<endl;
		}
	}

	aln1.save("problem.aln");
	BOOST_CHECK(v2Set.size() == aln1.aln.size());
}

BOOST_AUTO_TEST_CASE( ics_match_1 )
{
	Network net1("../optnetalign/tests/lccstest1.net");
	Network net2("../optnetalign/tests/lccstest2.net");
	Alignment aln(&net1,&net2,"../optnetalign/tests/lccstest.aln",nullptr);
	double ics = aln.ics();
	BOOST_CHECK(approxEqual(ics,0.8));
}

BOOST_AUTO_TEST_CASE( ics_match_2 )
{
	Network net1("../optnetalign/tests/lccstest1.net");
	Network net2("../optnetalign/tests/lccstest2.net");
	Alignment aln(&net1,&net2,"../optnetalign/tests/lccstest2.aln",nullptr);
	double ics = aln.ics();
	BOOST_CHECK(approxEqual(ics,0.0));
}

BOOST_AUTO_TEST_CASE( ics_match_3 )
{
	Network net1("../optnetalign/tests/newmetrica.net");
	Network net2("../optnetalign/tests/newmetricb.net");
	Alignment aln(&net1,&net2,"../optnetalign/tests/newmetric.aln",nullptr);
	double ics = aln.ics();
	BOOST_CHECK(approxEqual(ics,0.666666666666666));
}

//todo: the following 2 tests are broken while self-loops are ignored.
/*
BOOST_AUTO_TEST_CASE( ics_match_4 )
{
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln(&net1,&net2,"../optnetalign/tests/cg1.aln",nullptr);
	double ics = aln.ics();
	BOOST_CHECK(approxEqual(ics,0.8253323173313013));
}

BOOST_AUTO_TEST_CASE( ics_match_5 )
{
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln(&net1,&net2,"../optnetalign/tests/cg1partial.aln",nullptr);
	double ics = aln.ics();
	BOOST_CHECK(approxEqual(ics,0.920399546905571));
}
*/
BOOST_AUTO_TEST_CASE( consistent_after_mutate )
{
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln(&net1,&net2,"../optnetalign/tests/cg1.aln",nullptr);
	mt19937 g1(12);
	aln.mutate(g1, 0.3);
	set<node> v2Set;
	for(int i = 0; i < aln.aln.size(); i++){
		v2Set.insert(aln.aln[i]);
		BOOST_CHECK(aln.aln[i] >=0 && 
			        aln.aln[i] < net2.nodeToNodeName.size());
	}

	BOOST_CHECK(v2Set.size() == aln.aln.size());

}

BOOST_AUTO_TEST_CASE( consistent_after_crossover )
{
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln(&net1,&net2,"../optnetalign/tests/cg1.aln",nullptr);
	Alignment aln2(&net1,&net2,"../optnetalign/tests/cg1partial.aln",nullptr);
	mt19937 g1(12);
	Alignment child(g1, 0.3, aln, aln2);
	set<node> v2Set;
	for(int i = 0; i < aln.aln.size(); i++){
		v2Set.insert(aln.aln[i]);
		BOOST_CHECK(aln.aln[i] >=0 && 
			        aln.aln[i] < net2.nodeToNodeName.size());
	}

	BOOST_CHECK(v2Set.size() == aln.aln.size());

}

BOOST_AUTO_TEST_CASE( nondominated_sort )
{
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,"../optnetalign/tests/cg1.aln",nullptr);
	aln1.fitness.push_back(0.7);
	aln1.fitness.push_back(0.8);
	Alignment aln2(&net1,&net2,nullptr);
	aln2.fitness.push_back(0.75);
	aln2.fitness.push_back(0.75);
	Alignment aln3(&net1,&net2,nullptr);
	aln3.fitness.push_back(0.8);
	aln3.fitness.push_back(0.7);
	Alignment aln4(&net1,&net2,nullptr);
	aln4.fitness.push_back(0.4);
	aln4.fitness.push_back(0.6);
	Alignment aln5(&net1,&net2,nullptr);
	aln5.fitness.push_back(0.5);
	aln5.fitness.push_back(0.5);
	Alignment aln6(&net1,&net2,nullptr);
	aln6.fitness.push_back(0.6);
	aln6.fitness.push_back(0.4);
	vector<Alignment*> toSort;
	toSort.push_back(&aln3);
	toSort.push_back(&aln2);
	toSort.push_back(&aln4);
	toSort.push_back(&aln1);
	toSort.push_back(&aln6);
	toSort.push_back(&aln5);
	vector<vector<Alignment*> > fronts = nonDominatedSort(toSort);
	for(int i = 0; i < fronts.size(); i++){
		cout<<"FRONT "<<i<<endl;
		cout<<"------------"<<endl;
		for(int j = 0; j <fronts[i].size(); j++){
			for(int k = 0; k < fronts[i][j]->fitness.size(); k++){
				cout<<fronts[i][j]->fitness[k]<<" "<<endl;
			}
			cout<<endl;
		}
	}
	//ensure we have the right number of fronts
	BOOST_CHECK(fronts.size()==2);
	//and they are the right sizes
	BOOST_CHECK(fronts[0].size() == 3);
	BOOST_CHECK(fronts[1].size() == 3);
	//and they each contain the fitnesses we think they should
	for(int i = 0; i < 3; i++){
		BOOST_CHECK(fronts[0][i]->fitness[0] > 0.6);
		BOOST_CHECK(fronts[0][i]->fitness[1] > 0.6);
	}
	for(int i = 0; i < 3; i++){
		BOOST_CHECK(fronts[1][i]->fitness[0] < 0.7);
		BOOST_CHECK(fronts[1][i]->fitness[1] < 0.7);
	}
}

BOOST_AUTO_TEST_CASE( crowding_dist_assignment )
{
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,nullptr);
	aln1.fitness.push_back(0.6);
	aln1.fitness.push_back(0.9);
	Alignment aln2(&net1,&net2,nullptr);
	aln2.fitness.push_back(0.7);
	aln2.fitness.push_back(0.8);
	Alignment aln3(&net1,&net2,nullptr);
	aln3.fitness.push_back(0.75);
	aln3.fitness.push_back(0.75);
	Alignment aln4(&net1,&net2,nullptr);
	aln4.fitness.push_back(0.8);
	aln4.fitness.push_back(0.7);
	Alignment aln5(&net1,&net2,nullptr);
	aln5.fitness.push_back(0.9);
	aln5.fitness.push_back(0.6);
	vector<Alignment*> front;
	front.push_back(&aln1);
	front.push_back(&aln2);
	front.push_back(&aln3);
	front.push_back(&aln4);
	front.push_back(&aln5);
	setCrowdingDists(front);
	BOOST_CHECK(approxEqual(1.0,front[1]->crowdDist));
	BOOST_CHECK(approxEqual(0.666666666666666,front[2]->crowdDist));
	BOOST_CHECK(approxEqual(1.0,front[3]->crowdDist));
}

BOOST_AUTO_TEST_CASE( adjacency_list_correct ){
	Network net1("../optnetalign/tests/lccstest1.net");
	BOOST_CHECK(net1.degree(net1.nodeNameToNode.at("1")) == 1);
	BOOST_CHECK(net1.degree(net1.nodeNameToNode.at("2")) == 2);
	BOOST_CHECK(net1.degree(net1.nodeNameToNode.at("3")) == 2);
	BOOST_CHECK(net1.degree(net1.nodeNameToNode.at("4")) == 3);
}

BOOST_AUTO_TEST_CASE( fast_ics_works_total_crossover ){
	cout<<endl<<"BEGIN fast_ics_works_total_crossover"<<endl;
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,nullptr);
	Alignment aln2(&net1,&net2,nullptr);
	mt19937 g1(12);
	aln1.shuf(g1,false,false,true);
	aln2.shuf(g1,false,false,true);
	Alignment child(g1,0.2,aln1,aln2,true);
	double ics = child.ics();
	double fastICS = child.fastICS();
	cout<<"ICS is "<<ics<<endl;
	cout<<"Fast ICS is "<<fastICS<<endl;
	BOOST_CHECK(approxEqual(ics,fastICS));
}

BOOST_AUTO_TEST_CASE( fast_ics_works_partial_crossover ){
	cout<<endl<<"BEGIN fast_ics_works_partial_crossover"<<endl;
	Network net1("../optnetalign/tests/selflooptest.net");
	Network net2("../optnetalign/tests/selflooptest.net");
	Alignment aln1(&net1,&net2,nullptr);
	Alignment aln2(&net1,&net2,nullptr);
	mt19937 g1(12);
	//aln1.shuf(g1,false,false,false);
	//aln2.shuf(g1,false,false,false);
	Alignment child(g1,0.2,aln1,aln2,false);
	double ics = child.ics();
	double fastICS = child.fastICS();
	cout<<"ICS is "<<ics<<endl;
	cout<<"Fast ICS is "<<fastICS<<endl;
	BOOST_CHECK(approxEqual(ics,fastICS));
}

BOOST_AUTO_TEST_CASE( fast_ics_works_total_mutation ){
	cout<<endl<<"BEGIN fast_ics_works_total_mutation"<<endl;
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,nullptr);
	mt19937 g1(12);
	aln1.shuf(g1,false,false,true);
	aln1.mutate(g1,0.1,true);
	double ics = aln1.ics();
	double fastICS = aln1.fastICS();
	cout<<"ICS is "<<ics<<endl;
	cout<<"Fast ICS is "<<fastICS<<endl;
	BOOST_CHECK(approxEqual(ics,fastICS));	
}

BOOST_AUTO_TEST_CASE( fast_ics_works_total_repeat_mutation ){
	cout<<endl<<"BEGIN fast_ics_works_total_repeat_mutation"<<endl;
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,nullptr);
	mt19937 g1(12);
	aln1.shuf(g1,false,false,true);
	for(int i = 0; i < 100; i++){
		aln1.mutate(g1,0.1,true);
	}
	double ics = aln1.ics();
	double fastICS = aln1.fastICS();
	cout<<"ICS is "<<ics<<endl;
	cout<<"Fast ICS is "<<fastICS<<endl;
	BOOST_CHECK(approxEqual(ics,fastICS));	
}

BOOST_AUTO_TEST_CASE( fast_ics_works_partial_mutation ){
	cout<<endl<<"BEGIN fast_ics_works_partial_mutation"<<endl;
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,nullptr);
	mt19937 g1(4);
	aln1.shuf(g1,false,false,false);
	cout<<"BEFORE MUTATION"<<endl;

	aln1.mutate(g1,0.1,false);
	
	cout<<"AFTER MUTATION"<<endl;
	
	vector<int> badCounts = aln1.conservedCounts;
	//reinit counts to correct them
	for(int i = 0; i < aln1.actualSize; i++){
		aln1.initConservedCount(i,aln1.aln[i], aln1.alnMask[i]);
	}
	vector<int> goodCounts = aln1.conservedCounts;
	for(int i = 0; i < goodCounts.size(); i++){
		if(badCounts[i] != goodCounts[i]){
			cout<<i<<" (AKA "<<net1.nodeToNodeName.at(i)
				<<") count should be "<<goodCounts[i]
                <<" but is actually "<<badCounts[i]<<endl;
		}
	}
	aln1.conservedCounts = badCounts;
	double ics = aln1.ics();
	double fastICS = aln1.fastICS();
	cout<<"ICS is "<<ics<<endl;
	cout<<"Fast ICS is "<<fastICS<<endl;
	BOOST_CHECK(approxEqual(ics,fastICS));	
}

BOOST_AUTO_TEST_CASE( fast_ics_works_partial_repeat_mutation ){
	cout<<endl<<"BEGIN fast_ics_works_total_repeat_mutation"<<endl;
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,nullptr);
	mt19937 g1(12);
	aln1.shuf(g1,false,false,false);
	for(int i = 0; i < 100; i++){
		aln1.mutate(g1,0.1,false);
	}
	double ics = aln1.ics();
	double fastICS = aln1.fastICS();
	cout<<"ICS is "<<ics<<endl;
	cout<<"Fast ICS is "<<fastICS<<endl;
	BOOST_CHECK(approxEqual(ics,fastICS));	
}


BOOST_AUTO_TEST_CASE( fast_ics_works_after_init_total ){
	cout<<endl<<"BEGIN fast_ics_works_after_init_total"<<endl;
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,nullptr);
	mt19937 g1(12);
	aln1.shuf(g1,false,false,true);
	double ics = aln1.ics();
	double fastICS = aln1.fastICS();
	cout<<"ICS is "<<ics<<endl;
	cout<<"Fast ICS is "<<fastICS<<endl;
	BOOST_CHECK(approxEqual(ics,fastICS));	
}


BOOST_AUTO_TEST_CASE( fast_ics_works_after_init_partial ){
	cout<<endl<<"BEGIN fast_ics_works_after_init_partial"<<endl;
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,nullptr);
	mt19937 g1(123);
	aln1.shuf(g1,false,false,false);
	double ics = aln1.ics();
	double fastICS = aln1.fastICS();
	cout<<"ICS is "<<ics<<endl;
	cout<<"Fast ICS is "<<fastICS<<endl;
	BOOST_CHECK(approxEqual(ics,fastICS));	
}

BOOST_AUTO_TEST_CASE( do_swap_consistent_ics ){
	cout<<endl<<"BEGIN do_swap_consistent_ics"<<endl;
	Network net1("../optnetalign/tests/selflooptest.net");
	Network net2("../optnetalign/tests/selflooptest.net");
	Alignment aln1(&net1,&net2,nullptr);
	mt19937 g1(124);
	aln1.shuf(g1,false,false,true);
	
	cout<<"before doSwap:"<<endl;
	cout<<"Alignment is: "<<endl;
	for(int i = 0; i < aln1.actualSize; i++){
		cout<<net1.nodeToNodeName.at(i)<<" "
		    <<net2.nodeToNodeName.at(aln1.aln[i])<<" mask: "
		    <<aln1.alnMask[i]<<endl;
	}
	cout<<endl<<"Conserved counts are: "<<endl;
	for(int i = 0; i < aln1.conservedCounts.size(); i++){
		cout<<net1.nodeToNodeName.at(i)<<" "
		    <<aln1.conservedCounts[i]<<endl;
	}
	cout<<endl;
	
	double icsBefore = aln1.ics();
	double fastICSBefore = aln1.fastICS();
	cout<<"icsBefore: "<<icsBefore<<endl;
	cout<<"fastICSBefore: "<<fastICSBefore<<endl;
	BOOST_CHECK(approxEqual(icsBefore,fastICSBefore));
	auto dist = uniform_int_distribution<int>(0,aln1.aln.size()-1);
	node x = dist(g1);
	node y = dist(g1);
	//cout<<"x is "<<net1.nodeToNodeName.at(x)<<endl;
	//cout<<"y is "<<net1.nodeToNodeName.at(y)<<endl;
	aln1.doSwap(x,y);
	
	cout<<"after doSwap:"<<endl;
	cout<<"Alignment is: "<<endl;
	for(int i = 0; i < aln1.actualSize; i++){
		cout<<net1.nodeToNodeName.at(i)<<" "
		    <<net2.nodeToNodeName.at(aln1.aln[i])<<" mask: "
		    <<aln1.alnMask[i]<<endl;
	}
	cout<<endl<<"Conserved counts are: "<<endl;
	for(int i = 0; i < aln1.conservedCounts.size(); i++){
		cout<<net1.nodeToNodeName.at(i)<<" "
		    <<aln1.conservedCounts[i]<<endl;
	}
	cout<<endl;
	
	double ics = aln1.ics();
	double fastICS = aln1.fastICS();
	cout<<"ICS is "<<ics<<endl;
	cout<<"Fast ICS is "<<fastICS<<endl;
	BOOST_CHECK(approxEqual(ics,fastICS));	
}


BOOST_AUTO_TEST_CASE( hypothetical_swap_correct ){
	cout<<endl<<"BEGIN hypothetical_swap_correct"<<endl;
	Network net1("../optnetalign/tests/lccstest1.net");
	Network net2("../optnetalign/tests/lccstest2.net");
	Alignment aln1(&net1,&net2,"../optnetalign/tests/lccstest.aln",nullptr);
	mt19937 g1(123);
	vector<string> fitnessNames;
	fitnessNames.push_back("ICS");
	aln1.computeFitness(fitnessNames);

	uniform_int_distribution<int> randIndex(0,aln1.aln.size());


	for(int i = 0; i < 100; i++){
		node x = randIndex(g1);
		node y = x;
		while(y == x){
			y = randIndex(g1);
		}
		vector<double> delta = aln1.doSwapHypothetical(x,y);

		//get sum of conservedCounts before swapping
		int sumBefore = 0;
		for(int j = 0; j < aln1.conservedCounts.size(); j++){
			sumBefore += aln1.conservedCounts[j];
		}

		aln1.doSwap(x,y);

		int sumAfter = 0;
		for(int j = 0; j < aln1.conservedCounts.size(); j++){
			sumAfter += aln1.conservedCounts[j];
		}

		int diff = sumAfter - sumBefore;
		if(diff != 2*delta[0]){
			cout<<"diff: "<<diff<<endl;
			cout<<"2*delta[0]: "<<2*delta[0]<<endl;
			cout<<"x: "<<x<<endl;
			cout<<"y: "<<y<<endl;
		}
	}


}


BOOST_AUTO_TEST_CASE( v1unaligned_consistent_after_load ){
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln(&net1, &net2, "../optnetalign/tests/cg1partial.aln", nullptr);

	for(int i = 0; i < aln.alnMask.size(); i++){
		if(!aln.alnMask[i])
			BOOST_CHECK(aln.v1Unaligned.count(i));
		else
			BOOST_CHECK(!aln.v1Unaligned.count(i));
	}
}

BOOST_AUTO_TEST_CASE( v1unaligned_consistent_after_shuf ){
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln(&net1, &net2, nullptr);
	mt19937 g(12);
	aln.shuf(g, false, false, false);

	for(int i = 0; i < aln.alnMask.size(); i++){
		if(!aln.alnMask[i])
			BOOST_CHECK(aln.v1Unaligned.count(i));
		else
			BOOST_CHECK(!aln.v1Unaligned.count(i));
	}
}

BOOST_AUTO_TEST_CASE( v1unaligned_consistent_after_onBit ){
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln(&net1, &net2, nullptr);
	mt19937 g(12);
	aln.shuf(g, false, false, false);

	int offbit = 0;
	for(; offbit < aln.alnMask.size(); offbit++){
		if(!aln.alnMask[offbit])
			break;
	}

	aln.onBit(offbit);

	for(int i = 0; i < aln.alnMask.size(); i++){
		if(!aln.alnMask[i])
			BOOST_CHECK(aln.v1Unaligned.count(i));
		else
			BOOST_CHECK(!aln.v1Unaligned.count(i));
	}
}

BOOST_AUTO_TEST_CASE( greedy_match_aln_consistent ){
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	BLASTDict b = loadBLASTInfo(&net1,&net2,"../optnetalign/tests/cg1.sim");
	Alignment aln(&net1,&net2,&b);

	//check that if a node is in v1Unaligned, then it is unaligned
	for(auto n : aln.v1Unaligned){
		BOOST_CHECK(!aln.alnMask[n]);
	}

	//check that if a node is unaligned, it is in v1Unaligned
	for(int i = 0; i < aln.alnMask.size(); i++){
		BOOST_CHECK_EQUAL(aln.v1Unaligned.count(i),!aln.alnMask[i]);
	}

	//check that the aln is a permutation
	unordered_set<node> alnv2;

	for(int i = 0; i < aln.aln.size(); i++){
		alnv2.insert(aln.aln[i]);
	}

	BOOST_CHECK_EQUAL(alnv2.size(),aln.aln.size());
}

BOOST_AUTO_TEST_CASE( self_similarity_is_one){
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	BLASTDict b = loadBLASTInfo(&net1,&net2,"../optnetalign/tests/cg1.sim");
	Alignment aln(&net1,&net2,"../optnetalign/tests/cg1.aln",&b);
	BOOST_CHECK(approxEqual(alnSimilarity(&aln,&aln),1.0));

}

BOOST_AUTO_TEST_CASE(similarity_check){
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln1(&net1,&net2,"../optnetalign/tests/cg1.aln",nullptr);
	Alignment aln2(&net1,&net2,"../optnetalign/tests/cg1partial.aln",nullptr);
	double sim = alnSimilarity(&aln1,&aln2);
	BOOST_CHECK(approxEqual(sim,0.8891585761));
}
//____________________________________________________________________________//

// EOF