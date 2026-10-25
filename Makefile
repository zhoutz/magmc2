main:
	g++-16 -std=c++23 -lstdc++exp src/main.cpp -o build/main -I/opt/homebrew/include -L/opt/homebrew/lib -lgsl -lgslcblas -Wall
	build/main

.PHONY: fd11_fig7
fd11_fig7: build/fd11_fig7

build/fd11_fig7: src/fd11_fig7.cpp $(wildcard src/*.hpp)
	mkdir -p build
	g++-16 -std=c++23 -O3 -DNDEBUG src/fd11_fig7.cpp -o $@ -I/opt/homebrew/include -L/opt/homebrew/lib -lstdc++exp -lgsl -lgslcblas -Wall -Wextra

.PHONY: test_resonant
test_resonant: build/test_resonant
	build/test_resonant

build/test_resonant: tests/test_resonant.cpp $(wildcard src/*.hpp)
	mkdir -p build
	g++-16 -std=c++23 -O2 tests/test_resonant.cpp -o $@ -Isrc -I/opt/homebrew/include -L/opt/homebrew/lib -lstdc++exp -lgsl -lgslcblas -Wall -Wextra
