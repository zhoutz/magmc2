warning-flags = -Wall -Wextra -Wpedantic
compiler = g++-16
cxxflags = -std=c++23 -O3 -march=native $(warning-flags) -I/opt/homebrew/include
ldflags = -L/opt/homebrew/lib -lstdc++exp -lgsl -lgslcblas

.PHONY: main compile test benchmark validate

main: compile
	build/main

compile:
	mkdir -p build
	$(compiler) $(cxxflags) src/main.cpp -o build/main $(ldflags)

test:
	mkdir -p build
	$(compiler) $(cxxflags) -Isrc tests/test_thinning.cpp -o build/test_thinning $(ldflags)
	build/test_thinning

benchmark:
	python3 py/validate_thinning.py build
	python3 py/validate_thinning.py benchmark

validate:
	python3 py/validate_thinning.py all
