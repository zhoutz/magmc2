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






template <class DerivFunc>
struct StepperMagnus {
  static constexpr double EPS = std::numeric_limits<double>::epsilon();
  
  DerivFunc const& derivs;
  double x_old, h_old, x_max;
  // double3 y_old, y_new, y_err;
  double atol, rtol;

  StepperMagnus(DerivFunc const& derivs, double atol, double rtol)
      : derivs(derivs), atol(atol), rtol(rtol) {}


  void
  init(double x_init, double h_init, double x_end, double3 const& y_init) {
    if (!(std::isfinite(h_init) && h_init > 0)) {
      throw std::invalid_argument("h_init must be positive and finite");
    }
    x_old = x_init;
    h_old = h_init;
    x_max = x_end;
    // y_old = y_init;
    // derivs(x_old, y_old, dydx_old);
    
  }

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

  void
  do_step() {
    h_old = std::min(h_old, max_step);
    while (true) {
      if (std::abs(h_old) <= std::abs(x_old) * EPS) {
        throw std::runtime_error("stepsize underflow in StepperDopr5");
      }
      try_step();
      if (success(error()))
        break;
    }
    prepare_dense();
  }

  void
  update_old() {
    dydx_old = dydx_new;
    y_old = y_new;
    x_old += h_old;
    h_old = h_new;
  }


  void
  prepare_dense() {
    
  }

  double3
  dense_out(double x) const {
    
  }

};
