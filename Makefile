main:
	g++-16 -std=c++23 -lstdc++exp src/main.cpp -o build/main -I/opt/homebrew/include -L/opt/homebrew/lib -lgsl -lgslcblas -Wall
	build/main

.PHONY: fd11_fig7
fd11_fig7: build/fd11_fig7

build/fd11_fig7: src/fd11_fig7.cpp $(wildcard src/*.hpp)
	mkdir -p build
	g++-16 -std=c++23 -O3 -DNDEBUG -fopenmp src/fd11_fig7.cpp -o $@ -I/opt/homebrew/include -L/opt/homebrew/lib -lstdc++exp -lgsl -lgslcblas -Wall -Wextra

.PHONY: test-transport-bounds
test-transport-bounds: build/test_transport_bounds
	build/test_transport_bounds

build/test_transport_bounds: tests/transport_bounds.cpp $(wildcard src/*.hpp)
	mkdir -p build
	g++-16 -std=c++23 -O3 tests/transport_bounds.cpp -o $@ -Isrc -I/opt/homebrew/include -L/opt/homebrew/lib -lgsl -lgslcblas -Wall -Wextra

.PHONY: test-transport-cache
test-transport-cache: build/test_transport_cache
	build/test_transport_cache

build/test_transport_cache: tests/transport_cache.cpp src/transport_cache.hpp src/photon_stream.hpp
	mkdir -p build
	g++-16 -std=c++23 -O3 tests/transport_cache.cpp -o $@ -Isrc -Wall -Wextra

.PHONY: test-transport-driver
test-transport-driver: build/fd11_fig7
	python3 tests/test_transport_driver.py -v
