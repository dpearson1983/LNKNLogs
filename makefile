CXX = g++ -std=c++11
CXXFLAGS = -lgsl -lgslcblas -lm -lfftw3 -lfftw3_omp -fopenmp
CXXOPTIMIZE = -march=native -O3
DEPS = include/harppi.h include/lognormal.h include/tpods.h include/constants.h

lnknlogs: main.cpp obj/lognormal.o obj/harppi.o $(DEPS)
	$(CXX) $(CXXFLAGS) $(CXXOPTIMIZE) -o LNKNLogs main.cpp obj/harppi.o obj/lognormal.o
	
obj/harppi.o: source/harppi.cpp include/harppi.h
	$(CXX) $(CXXOPTIMIZE) -c -o obj/harppi.o source/harppi.cpp
	
obj/lognormal.o: source/lognormal.cpp include/lognormal.h
	$(CXX) $(CXXFLAGS) $(CXXOPTIMIZE) -c -o obj/lognormal.o source/lognormal.cpp
	
clean:
	rm obj/harppi.o obj/lognormal.o LNKNLogs
