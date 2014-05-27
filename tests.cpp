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
#include "Alignment.h"
#include "Network.h"
#include "blastinfo.h"

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
    Alignment a= Alignment(net1,net2,"../optnetalign/tests/small.aln");
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/small.bitScores");
    double sumbit = a.sumBLAST(net1,net2,b);
    //std::cout<<"sumbit is "<<sumbit<<std::endl;
    BOOST_CHECK(approxEqual(2.0,a.sumBLAST(net1,net2,b)));
    // reports 'error in "test2": check i == 2 failed [0 != 2]'
    //BOOST_CHECK_EQUAL( i, 2 );

    //BOOST_CHECK_EQUAL( i, 0 );
}

BOOST_AUTO_TEST_CASE( bitscoreSum_match_2 )
{
	Network net1("../optnetalign/tests/small.net");
    Network net2("../optnetalign/tests/small.net");
    Alignment a2 = Alignment(net1,net2,"../optnetalign/tests/small2.aln");
    for(int i = 0; i < a2.aln.size(); i++){
    	cout<<"node "<<i<<" is aligned to "<<a2.aln[i]<<endl;
    }
    BLASTDict b = loadBLASTInfo(&net1,&net2,
    	                        "../optnetalign/tests/small.bitScores");
    double sumbit2 = a2.sumBLAST(net1,net2,b);
    std::cout<<"sumbit2 is "<<sumbit2<<std::endl;
    BOOST_CHECK(approxEqual(sumbit2,0.0));
}

BOOST_AUTO_TEST_CASE( save_load_inverses )
{
	Network net1("../optnetalign/tests/lccstest1.net");
	Network net2("../optnetalign/tests/lccstest2.net");
	Alignment aln1(net1,net2,"../optnetalign/tests/lccstest.aln");
	aln1.save(net1,net2,"../optnetalign/tests/testSavecpp.aln");
	Alignment aln2(net1,net2,"../optnetalign/tests/testSavecpp.aln");

	for(int i = 0; i < net1.nodeToNodeName.size(); i++){
		BOOST_CHECK(aln1.aln[i] == aln2.aln[i]);
	}
}

BOOST_AUTO_TEST_CASE( rand_aln_consistent )
{
	Network net1("../optnetalign/tests/lccstest1.net");
	Network net2("../optnetalign/tests/lccstest2.net");
	Alignment aln1(net1,net2);
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
	Alignment aln1(net1,net2,"../optnetalign/tests/cg1.aln");
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

	aln1.save(net1,net2,"problem.aln");
	BOOST_CHECK(v2Set.size() == aln1.aln.size());
}

BOOST_AUTO_TEST_CASE( ics_match_1 )
{
	Network net1("../optnetalign/tests/lccstest1.net");
	Network net2("../optnetalign/tests/lccstest2.net");
	Alignment aln(net1,net2,"../optnetalign/tests/lccstest.aln");
	double ics = aln.ics(net1,net2);
	BOOST_CHECK(approxEqual(ics,0.8));
}

BOOST_AUTO_TEST_CASE( ics_match_2 )
{
	Network net1("../optnetalign/tests/lccstest1.net");
	Network net2("../optnetalign/tests/lccstest2.net");
	Alignment aln(net1,net2,"../optnetalign/tests/lccstest2.aln");
	double ics = aln.ics(net1,net2);
	BOOST_CHECK(approxEqual(ics,0.0));
}

BOOST_AUTO_TEST_CASE( ics_match_3 )
{
	Network net1("../optnetalign/tests/newmetrica.net");
	Network net2("../optnetalign/tests/newmetricb.net");
	Alignment aln(net1,net2,"../optnetalign/tests/newmetric.aln");
	double ics = aln.ics(net1,net2);
	BOOST_CHECK(approxEqual(ics,0.666666666666666));
}

BOOST_AUTO_TEST_CASE( ics_match_4 )
{
	Network net1("../optnetalign/tests/cg1a.net");
	Network net2("../optnetalign/tests/cg1b.net");
	Alignment aln(net1,net2,"../optnetalign/tests/cg1.aln");
	double ics = aln.ics(net1,net2);
	BOOST_CHECK(approxEqual(ics,0.8253323173313013));
}

//____________________________________________________________________________//

// EOF