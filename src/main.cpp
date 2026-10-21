
#include "find_root.hpp"
#include <numeric>
#include <print>

double
f(double x) {
  return x - 1e-100;
}

int
main() {
  double a = 0, b = 1;
  // double root = bisect_root(f, a, b, 0, f(a), f(b));
  double root = ridders_root(f, a, b, 0, f(a), f(b));
  std::print("Root: {}\n", root);
}
