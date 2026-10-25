#pragma once
#include "doubles.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

// Local degree-four interpolation, checked at four independent points.
// A failed leaf always uses the original field evaluator. No state or
// scattering rates are interpolated here. Reset at every geodesic interval.
struct RotationCache {
  struct Leaf {
    double a, b;
    std::array<double3, 5> c;
    bool exact;
    double3
    at(double x) const {
      double t = 4 * (x - a) / (b - a);
      return c[0] +
             t * (c[1] + (t - 1) * (c[2] + (t - 2) * (c[3] + (t - 3) * c[4])));
    }
  };
  std::vector<Leaf> leaves;
  double tolerance = 1e-9;
  template <class F>
  void
  build(double a, double b, F const& fn, int depth = 0) {
    std::array<double3, 5> v;
    double scale = 0;
    for (int i = 0; i < 5; ++i) {
      v[i] = fn(std::lerp(a, b, i / 4.0));
      scale = std::max(scale, v[i].length());
    }
    Leaf leaf{a, b, {}, false};
    double factorial = 1;
    for (int i = 0; i < 5; ++i) {
      if (i)
        factorial *= i;
      leaf.c[i] = v[0] / factorial;
      for (int j = 0; j < 4 - i; ++j)
        v[j] = v[j + 1] - v[j];
    }
    bool valid = true;
    for (int i = 0; i < 4; ++i) {
      double x = std::lerp(a, b, (2 * i + 1) / 8.0);
      if ((fn(x) - leaf.at(x)).length() > tolerance * std::max(scale, 1e-100))
        valid = false;
    }
    double mid = std::midpoint(a, b);
    if (valid || depth == 4 || !(a < mid && mid < b)) {
      leaf.exact = !valid;
      leaves.push_back(leaf);
    } else {
      build(a, mid, fn, depth + 1);
      build(mid, b, fn, depth + 1);
    }
  }
  template <class F>
  double3
  at(double x, F const& fn) const {
    auto it = std::lower_bound(leaves.begin(), leaves.end(), x,
                               [](Leaf const& l, double t) { return l.b < t; });
    if (it == leaves.end() || x < it->a || it->exact)
      return fn(x);
    return it->at(x);
  }
};
