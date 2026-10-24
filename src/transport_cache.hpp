#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

// Local interpolation of r, mu, x_e, twist, B1, B2 and the energy factor.
// Independent interleaved points estimate error; this is not a rigorous
// uniform error bound. Split difficult intervals and use exact evaluation
// when the subdivision limit is reached. Rebuild after every geodesic step.
struct TransportCache {
  using Values = std::array<double, 7>;
  struct Leaf {
    double a, b;
    std::array<Values, 5> coeff;
    bool exact = false;

    Values at(double l) const {
      double t = 4 * (l - a) / (b - a);
      Values y{};
      for (int k = 0; k < 7; ++k)
        y[k] = coeff[0][k] + t * (coeff[1][k] + (t - 1) *
            (coeff[2][k] + (t - 2) * (coeff[3][k] + (t - 3) * coeff[4][k])));
      return y;
    }
  };

  std::vector<Leaf> leaves;
  double tolerance = 1e-7;

  template <class Sample>
  void build(double a, double b, Sample const& sample, int depth = 0) {
    std::array<Values, 5> v;
    Values scales{};
    for (int i = 0; i < 5; ++i) {
      v[i] = sample(std::lerp(a, b, i / 4.0));
      for (int k = 0; k < 7; ++k)
        scales[k] = std::max(scales[k], std::abs(v[i][k]));
    }
    scales[1] = 1; // Absolute error in the direction cosine.
    // A shared transverse scale remains regular when either component is zero.
    scales[4] = scales[5] = std::hypot(scales[4], scales[5]);
    Leaf leaf{a, b, {}};
    double factorial = 1;
    for (int i = 0; i < 5; ++i) {
      if (i) factorial *= i;
      for (int k = 0; k < 7; ++k) leaf.coeff[i][k] = v[0][k] / factorial;
      for (int j = 0; j < 4 - i; ++j)
        for (int k = 0; k < 7; ++k) v[j][k] = v[j + 1][k] - v[j][k];
    }
    bool valid = true;
    for (int i = 0; i < 4; ++i) {
      double l = std::lerp(a, b, (2 * i + 1) / 8.0);
      auto exact = sample(l), fit = leaf.at(l);
      for (int k = 0; k < 7; ++k)
        if (!std::isfinite(fit[k]) || !std::isfinite(exact[k]) ||
            std::abs(exact[k] - fit[k]) >
                tolerance * std::max(1e-100, scales[k]))
          valid = false;
    }
    double mid = std::midpoint(a, b);
    if (valid || depth == 8 || !(a < mid && mid < b)) {
      leaf.exact = !valid;
      leaves.push_back(leaf);
    } else {
      build(a, mid, sample, depth + 1);
      build(mid, b, sample, depth + 1);
    }
  }

  template <class Sample>
  Values at(double l, Sample const& sample) const {
    if (leaves.empty() || l < leaves.front().a || l > leaves.back().b)
      return sample(l);
    auto it = std::lower_bound(leaves.begin(), leaves.end(), l,
        [](auto const& leaf, double x) { return leaf.b < x; });
    return it->exact ? sample(l) : it->at(l);
  }
};
