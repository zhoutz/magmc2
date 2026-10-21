// #include "doubles.hpp"

#include <numeric>
#include <print>
int
main() {
  double a = 1,b=10;
  while(true){
    double c = std::midpoint(a, b);
    if (c == a || c == b) break;
    b = c;
  }
  std::print("{} {} {}\n", a, b, std::midpoint(a, b));
}
