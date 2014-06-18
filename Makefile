debug:
	clang++ main.cpp Network.cpp blastinfo.cpp Alignment.cpp nsga-ii.cpp localSearch.cpp goc.cpp -o optnetalign -stdlib=libc++ -g -Wall -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -I /usr/local/Cellar/tbb/4.2.4/include -L /usr/local/Cellar/tbb/4.2.4/lib -L /usr/local/Cellar/boost/1.53.0/lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt -lboost_program_options-mt -pthread -ltbb
optimized:
	icc main.cpp Network.cpp blastinfo.cpp Alignment.cpp nsga-ii.cpp localSearch.cpp goc.cpp -o optnetalign -g -prof-use -ipo -prof-dir. -use-intel-optimized-headers -tbb -par-num-threads=16 -opt-mem-layout-trans=3 -opt-subscript-in-range -ansi-alias -xHost -unroll-aggressive -opt-calloc -no-prec-div -stdlib=libc++ -O3 -Wall -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt -lboost_program_options-mt -pthread

profile:
	icpc main.cpp Network.cpp blastinfo.cpp Alignment.cpp nsga-ii.cpp localSearch.cpp goc.cpp -o optnetalign -O2 -prof-gen -prof-dir. -tbb -ldl -Wall -std=c++11 -g -stdlib=libc++ -I /usr/include/ -L /usr/lib/x86_64-linux-gnu/ -lboost_filesystem-mt -lboost_system-mt -lboost_program_options-mt -lboost_thread-mt	
ubuntu:
	icc main.cpp Network.cpp blastinfo.cpp Alignment.cpp nsga-ii.cpp localSearch.cpp goc.cpp -o optnetalign -prof-use -ipo -prof-dir. -use-intel-optimized-headers -tbb -par-num-threads=8 -opt-mem-layout-trans=3 -opt-subscript-in-range -ansi-alias -fno-alias -O3 -xHost -unroll-aggressive -opt-calloc -no-prec-div -Wall -std=c++11 -I /usr/include/ -L /usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_system -lboost_program_options -lboost_thread	

test:
	clang++ tests.cpp blastinfo.cpp Alignment.cpp Network.cpp nsga-ii.cpp localSearch.cpp goc.cpp -o test -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_unit_test_framework-mt

testubuntu:
	clang++ tests.cpp blastinfo.cpp Alignment.cpp Network.cpp nsga-ii.cpp localSearch.cpp goc.cpp -o test -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_unit_test_framework

debugubuntu:
	icc main.cpp Network.cpp blastinfo.cpp Alignment.cpp nsga-ii.cpp localSearch.cpp goc.cpp -o optnetalign -prof-gen -prof-dir/home/connor/Drobox/profiled -tbb -ldl -Wall -std=c++11 -I /usr/include/ -L /usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_system -lboost_program_options -lboost_thread	


localTest:
	icc localTest.cpp Network.cpp blastinfo.cpp Alignment.cpp nsga-ii.cpp localSearch.cpp goc.cpp -o localTest -stdlib=libc++ -O3 -Wall -tbb -std=c++11 -I /usr/local/Cellar/boost/1.53.0/include/ -L /usr/local/Cellar/boost/1.53.0/lib -lboost_thread-mt -lboost_filesystem-mt -lboost_system-mt -lboost_program_options-mt -pthread

tempubuntu:
	icc temp.cpp Network.cpp blastinfo.cpp Alignment.cpp nsga-ii.cpp localSearch.cpp goc.cpp -o temp -prof-gen -prof-dir/home/connor/Drobox/profiled -tbb -ldl -Wall -std=c++11 -I /usr/include/ -L /usr/lib/x86_64-linux-gnu/ -lboost_filesystem -lboost_system -lboost_program_options -lboost_thread	
