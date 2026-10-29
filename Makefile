warning-flags = -Wall -Wextra -Wpedantic

HDF5_PREFIX := $(shell brew --prefix hdf5)
HDF5_FLAGS := -I$(HDF5_PREFIX)/include -L$(HDF5_PREFIX)/lib -lhdf5_cpp -lhdf5
GSL_PREFIX := $(shell brew --prefix gsl)
GSL_FLAGS := -I$(GSL_PREFIX)/include -L$(GSL_PREFIX)/lib -lgsl -lgslcblas
HIGHFIVE_PREFIX := $(shell pwd)/ext/highfive
HIGHFIVE_FLAGS := -I$(HIGHFIVE_PREFIX)/include

CXX_FLAGS := -std=c++23 -pthread -O3 -march=native 

HEADERS := $(wildcard src/*.hpp)
.PHONY: main main.1 test diagnostics benchmark-reference

main: build/main

build/main: src/main.cpp $(HEADERS)
	mkdir -p build
	g++-16 $(CXX_FLAGS) src/main.cpp -o build/main $(HDF5_FLAGS) $(GSL_FLAGS) $(HIGHFIVE_FLAGS) $(warning-flags)

main.1: build/main.1
	build/main.1

build/main.1: src/main.1.cpp $(HEADERS)
	mkdir -p build
	g++-16 $(CXX_FLAGS) src/main.1.cpp -o build/main.1 $(HDF5_FLAGS) $(GSL_FLAGS) $(HIGHFIVE_FLAGS) $(warning-flags)

test: build/test_numerics build/main
	build/test_numerics
	python3 benchmarks/check_cli.py build/main

build/test_numerics: benchmarks/test_numerics.cpp $(HEADERS)
	mkdir -p build
	g++-16 $(CXX_FLAGS) benchmarks/test_numerics.cpp -Isrc -o $@ $(GSL_FLAGS) $(warning-flags)

diagnostics: build/diagnostic

build/diagnostic: src/main.cpp $(HEADERS)
	mkdir -p build
	g++-16 $(CXX_FLAGS) -DMAGMC_DIAGNOSTICS src/main.cpp -o $@ $(GSL_FLAGS) $(warning-flags)

benchmark-reference:
	python3 benchmarks/run.py build stage00
	python3 benchmarks/run.py build original
