CXX = g++ -std=c++11
CXXFLAGS = -lgsl -lgslcblas -lm -lfftw3 -lfftw3_omp -lCCfits -lcfitsio -fopenmp
CXXOPTIMIZE = -march=native -mtune=native -O3
DEPS = include/harppi.h include/lognormal.h include/tpods.h include/constants.h

lnknlogs: main.cpp lognormal harppi file_io cosmology $(DEPS)
	$(CXX) $(CXXFLAGS) $(CXXOPTIMIZE) -o LNKNLogs main.cpp obj/*.o
	
harppi: source/harppi.cpp include/harppi.h
	mkdir -p obj
	$(CXX) $(CXXOPTIMIZE) -c -o obj/harppi.o source/harppi.cpp
	
lognormal: source/lognormal.cpp include/lognormal.h
	mkdir -p obj
	$(CXX) $(CXXFLAGS) $(CXXOPTIMIZE) -c -o obj/lognormal.o source/lognormal.cpp
	
file_io:
	mkdir -p obj
	$(CXX) $(CXXFLAGS) $(CXXOPTIMIZE) -c -o obj/file_io.o source/file_io.cpp
	
cosmology:
	mkdir -p obj
	$(CXX) $(CXXFLAGS) $(CXXOPTIMIZE) -c -o obj/cosmology.o source/cosmology.cpp
	
clean:
	rm obj/*.o LNKNLogs
