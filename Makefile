debug:
	clang++ main.cpp Network.cpp blastinfo.cpp Alignment.cpp -o optnetalign -stdlib=libc++ -g -Wall -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt -lboost_program_options-mt -pthread

optimized:
	clang++ main.cpp Network.cpp blastinfo.cpp Alignment.cpp -o optnetalign -stdlib=libc++ -Os -Wall -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt -lboost_program_options-mt -pthread

test:
	clang++ tests.cpp blastinfo.cpp Alignment.cpp Network.cpp -o test -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_unit_test_framework-mt