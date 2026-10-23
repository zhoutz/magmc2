CXX = g++-16
CPPFLAGS = -I src -I/opt/homebrew/include
CXXFLAGS = -std=c++23 -O2 -Wall -Wextra -pedantic
LDFLAGS = -L/opt/homebrew/lib
LDLIBS = -lgsl -lgslcblas -lstdc++exp

.PHONY: main test
main: build/main
	./build/main

build/main: src/main.cpp $(wildcard src/*.hpp) | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

build/test_%: tests/%.cpp $(wildcard src/*.hpp) src/main.cpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDFLAGS) $(LDLIBS)

build:
	mkdir -p build

test: build/test_find_root build/test_magnus build/test_initialize build/test_evolution
	./build/test_find_root
	./build/test_magnus
	./build/test_initialize
	./build/test_evolution
