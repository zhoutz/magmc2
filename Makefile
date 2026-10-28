warning-flags = -Wall -Wextra -Wpedantic

HDF5_PREFIX := $(shell brew --prefix hdf5)
HDF5_FLAGS := -I$(HDF5_PREFIX)/include -L$(HDF5_PREFIX)/lib -lhdf5_cpp -lhdf5
GSL_PREFIX := $(shell brew --prefix gsl)
GSL_FLAGS := -I$(GSL_PREFIX)/include -L$(GSL_PREFIX)/lib -lgsl -lgslcblas
HIGHFIVE_PREFIX := $(shell pwd)/ext/highfive
HIGHFIVE_FLAGS := -I$(HIGHFIVE_PREFIX)/include

main:
	g++-16 -std=c++23 -lstdc++exp src/main.cpp -o build/main $(HDF5_FLAGS) $(GSL_FLAGS) $(HIGHFIVE_FLAGS) $(warning-flags) -O3 -march=native 
# 	build/main

main.1:
	g++-16 -std=c++23 -pthread -lstdc++exp src/main.1.cpp -o build/main.1 $(HDF5_FLAGS) $(GSL_FLAGS) $(HIGHFIVE_FLAGS) $(warning-flags) -O3 -march=native
	build/main.1
	
