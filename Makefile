all:
	clang++ main.cpp -o optnetalign -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_filesystem-mt -lboost_system-mt