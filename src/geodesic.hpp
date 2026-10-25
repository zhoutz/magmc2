#pragma once

#include "constants.hpp"
#include "dopr5.hpp"
#include <array>
#include <cmath>

struct Geodesic {
  struct Options {
    double orbit_atol = 1e-10;
    double orbit_rtol = 1e-10;
    double escape_factor = 1e5;
    Options() : orbit_atol(1e-10), orbit_rtol(1e-10), escape_factor(1e5) {}
  } options;

  using State = std::array<double, 3>;

  double rs, R_star;
  double max_step_fraction = 0.1;
  StepperDopr5<3, Geodesic> stepper;

  Geodesic(double M_star, double R_star, Options options = {})
      : options(options), rs(M_star * schwarzschild_radius_of_sun_in_km),
        R_star(R_star), stepper(*this, options.orbit_atol, options.orbit_rtol) {
    // absorb at the surface
    stepper.add_event(
        [R_star](double, State const& y) { return y[0] - R_star; });
    // escape at escape_factor*R_star
    stepper.add_event([R_star, options](double, State const& y) {
      return y[0] - options.escape_factor * R_star;
    });
  }

  void
  init(double r0, double psi0, double alpha0) {
    stepper.init(0.0, 1e-3 * R_star, {r0, psi0, alpha0});
  }

  double
  l_old() const {
    return stepper.x_old;
  }
  double
  l_new() const {
    return stepper.x_old + stepper.h_old;
  }

  int
  step() {
    stepper.do_step(max_step_fraction * stepper.y_old[0]);
    return stepper.detect_event();
  }

  State
  interp(double l) const {
    return stepper.interp(l);
  }

  void
  operator()(double, State const& y, State& dydx) const {
    auto [r, psi, alpha] = y;
    double L = std::sqrt(1 - rs / r);
    dydx[0] = L * std::cos(alpha);
    dydx[1] = std::sin(alpha) / r;
    dydx[2] = -std::sin(alpha) / (r * L) * (1 - 3 * rs / (2 * r));
  }

  Geodesic(Geodesic const&) = delete;
  Geodesic& operator=(Geodesic const&) = delete;
  Geodesic(Geodesic&&) = delete;
  Geodesic& operator=(Geodesic&&) = delete;
};
