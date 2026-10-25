warning-flags = -Wall -Wextra -Wpedantic

main:
	g++-16 -std=c++23 -lstdc++exp src/main.cpp -o build/main -I/opt/homebrew/include -L/opt/homebrew/lib -lgsl -lgslcblas -O3 -march=native $(warning-flags)
	build/main
