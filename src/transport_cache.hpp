#pragma once
#include "doubles.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <vector>

// Only coefficients of the transport equations, never an interpolated
// scattering state. Accepted events always reconstruct the full exact Basis.
struct TransportCoefficients {
  double Bmag, energy_local, J, mu;
  double3 omega;
};

struct OpticalCoefficients {
  double Bmag, energy_local, J, mu;
};

class TransportCache {
  using Values = std::array<double, 7>;
  struct Segment {
    double left, right;
    std::array<Values, 5> differences;
    bool direct = false;
    template <int N = 7>
    std::array<double, N>
    at(double l) const {
      double u = 4 * ((l - left) / (right - left));
      std::array<double, N> out;
      for (int j = 0; j < N; ++j)
        out[j] = differences[0][j] +
                 u * (differences[1][j] +
                      (u - 1) * (differences[2][j] / 2 +
                                 (u - 2) * (differences[3][j] / 6 +
                                            (u - 3) * differences[4][j] / 24)));
      return out;
    }
  };
  std::vector<Segment> segments;
  double tolerance;
  bool rotation_coefficients = true;

  template <class F>
  static
  Values
  sample(F const& exact, double l) {
    auto b = exact(l);
    return {b.Bmag, b.energy_local, b.J, b.mu, b.omega.x, b.omega.y, b.omega.z};
  }

  template <class F>
  void
  add(double left, double right, F const& exact, int depth) {
    Segment s{left, right, {}};
    std::array<Values, 9> values;
    for (int k = 0; k < 9; ++k)
      values[k] = sample(exact, std::lerp(left, right, k / 8.));
    auto d = std::array<Values, 5>{values[0], values[2], values[4], values[6],
                                   values[8]};
    for (int k = 0; k < 5; ++k) {
      s.differences[k] = d[0];
      for (int i = 0; i < 4 - k; ++i)
        for (int j = 0; j < 7; ++j)
          d[i][j] = d[i + 1][j] - d[i][j];
    }
    Values scale{};
    int active = rotation_coefficients ? 7 : 4;
    for (auto const& v : values)
      for (int j = 0; j < active; ++j)
        scale[j] = std::max(scale[j], std::abs(v[j]));
    scale[3] = 1;
    if (rotation_coefficients)
      scale[4] = scale[5] = scale[6] = std::hypot(scale[4], scale[5], scale[6]);
    bool accurate = true;
    for (int k : {1, 3, 5, 7}) {
      auto estimate = s.at(std::lerp(left, right, k / 8.));
      for (int j = 0; j < active; ++j)
        accurate &= std::abs(estimate[j] - values[k][j]) <=
                    tolerance * std::max(scale[j], 1e-100);
    }
    if (!accurate && depth < 4) {
      double mid = std::midpoint(left, right);
      add(left, mid, exact, depth + 1);
      add(mid, right, exact, depth + 1);
    } else {
      s.direct = !accurate;
      segments.push_back(s);
    }
  }

public:
  explicit TransportCache(double tolerance = 1e-5) : tolerance(tolerance) {}
  template <class F>
  void
  reset(double left, double right, F const& exact, bool need_rotation = true) {
    segments.clear();
    rotation_coefficients = need_rotation;
    if (right > left)
      add(left, right, exact, 0);
  }
  template <class F>
  OpticalCoefficients
  optical_at(double l, F const& exact) const {
    auto it = std::lower_bound(
        segments.begin(), segments.end(), l,
        [](Segment const& s, double l) { return s.right < l; });
    if (it == segments.end() || l < it->left || it->direct) {
      auto b = exact(l);
      return {b.Bmag, b.energy_local, b.J, b.mu};
    }
    auto v = it->template at<4>(l);
    if (!(v[0] > 0 && v[1] > 0 && v[2] >= 0 && std::abs(v[3]) <= 1)) {
      auto b = exact(l);
      return {b.Bmag, b.energy_local, b.J, b.mu};
    }
    return {v[0],v[1],v[2],v[3]};
  }
  template <class F>
  TransportCoefficients
  at(double l, F const& exact) const {
    auto it = std::lower_bound(
        segments.begin(), segments.end(), l,
        [](Segment const& s, double l) { return s.right < l; });
    if (it == segments.end() || l < it->left || it->direct) {
      auto b = exact(l);
      return {b.Bmag, b.energy_local, b.J, b.mu, b.omega};
    }
    auto v = it->at(l);
    // Positive endpoint samples do not guarantee a positive polynomial.
    if (!(v[0] > 0 && v[1] > 0 && v[2] >= 0 && std::abs(v[3]) <= 1)) {
      auto b = exact(l);
      return {b.Bmag, b.energy_local, b.J, b.mu, b.omega};
    }
    return {v[0], v[1], v[2], v[3], {v[4], v[5], v[6]}};
  }
};
