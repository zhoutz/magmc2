warning-flags = -Wall -Wextra -Wpedantic

.PHONY: main test-transport benchmark-transport

main:
	g++-16 -std=c++23 -lstdc++exp src/main.cpp -o build/main -I/opt/homebrew/include -L/opt/homebrew/lib -lgsl -lgslcblas -O3 -march=native $(warning-flags)
	build/main

test-transport:
	mkdir -p build/single_thread_optimization
	g++-16 -std=c++23 tests/test_transport.cpp -Isrc -o build/single_thread_optimization/test_transport -I/opt/homebrew/include -L/opt/homebrew/lib -lstdc++exp -lgsl -lgslcblas -O3 -march=native $(warning-flags)
	build/single_thread_optimization/test_transport
	python3 -m unittest discover -s py -p 'test_benchmark_transport.py'

benchmark-transport:
	python3 py/benchmark_transport.py --production --warmup --repeats 5
