warning-flags = -Wall -Wextra -Wpedantic

HDF5_PREFIX := $(shell brew --prefix hdf5)
HDF5_FLAGS := -I$(HDF5_PREFIX)/include -L$(HDF5_PREFIX)/lib -lhdf5_cpp -lhdf5
GSL_PREFIX := $(shell brew --prefix gsl)
GSL_FLAGS := -I$(GSL_PREFIX)/include -L$(GSL_PREFIX)/lib -lgsl -lgslcblas
HIGHFIVE_PREFIX := $(shell pwd)/ext/highfive
HIGHFIVE_FLAGS := -I$(HIGHFIVE_PREFIX)/include

CXX_FLAGS := -std=c++23 -pthread -O3 -march=native 

main:
	g++-16 $(CXX_FLAGS) src/main.cpp -o build/main $(HDF5_FLAGS) $(GSL_FLAGS) $(HIGHFIVE_FLAGS) $(warning-flags)
# 	build/main

main.1:
	g++-16 $(CXX_FLAGS) src/main.1.cpp -o build/main.1 $(HDF5_FLAGS) $(GSL_FLAGS) $(HIGHFIVE_FLAGS) $(warning-flags)

.PHONY: main main.1 test-qed test-qed-smoke
test-qed:
	g++-16 $(CXX_FLAGS) -Isrc tests/test_qed.cpp -o build/test_qed $(GSL_FLAGS) $(warning-flags)
	build/test_qed full
	g++-16 $(CXX_FLAGS) -Isrc tests/test_qed_output.cpp -o build/test_qed_output $(HDF5_FLAGS) $(HIGHFIVE_FLAGS) $(warning-flags)
	build/test_qed_output

test-qed-smoke: test-qed
	build/test_qed smoke
	
