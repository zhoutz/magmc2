#pragma once

#include "polarization.hpp"

struct Rotation {
  double w = 1;
  double3 v{0, 0, 0};

  static
  Rotation
  exponential(double3 a) {
    double angle = a.length();
    double scale =
        angle < 1e-10 ? 0.5 - angle * angle / 48 : std::sin(angle / 2) / angle;
    return {std::cos(angle / 2), scale * a};
  }

  static
  Rotation
  compose(Rotation after, Rotation before) {
    return {after.w * before.w - dot(after.v, before.v),
            after.w * before.v + before.w * after.v + cross(after.v, before.v)};
  }

  Stokes
  apply(Stokes s) const {
    auto t = 2 * cross(v, s);
    return s + w * t + cross(v, t);
  }
};

template <class F>
Rotation
magnus(F const& omega, double x, double h) {
  constexpr double g = 0.28867513459481288225;
  constexpr double a = 0.53867513459481288225, b = -0.03867513459481288225;
  auto c1 = omega(x + (0.5 - g) * h), c2 = omega(x + (0.5 + g) * h);
  auto first = Rotation::exponential(h * (a * c1.omega() + b * c2.omega()));
  auto second = Rotation::exponential(h * (b * c1.omega() + a * c2.omega()));
  return Rotation::compose(second, first);
}
