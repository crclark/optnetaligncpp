debug:
	clang++ main.cpp Network.cpp blastinfo.cpp Alignment.cpp -o optnetalign -g -Wall -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_filesystem-mt -lboost_system-mt -lboost_program_options-mt

optimized:
	clang++ main.cpp Network.cpp blastinfo.cpp Alignment.cpp -o optnetalign -Os -Wall -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_filesystem-mt -lboost_system-mt -lboost_program_options-mt

test:
	clang++ tests.cpp blastinfo.cpp Alignment.cpp Network.cpp -o test -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_unit_test_framework-mt