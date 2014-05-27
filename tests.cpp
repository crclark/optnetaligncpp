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
//____________________________________________________________________________//

// EOF