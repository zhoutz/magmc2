#pragma once
#include "diagnostics.hpp"

#include "polarization.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>

struct Rotation {
  double w = 1;
  double3 v{0, 0, 0};

  static
  Rotation
  exponential(double3 a) {
    double angle = a.length();
    if (!std::isfinite(angle))
      throw std::runtime_error("Non-finite rotation in StepperMagnus");
    double scale =
        angle < 1e-10 ? 0.5 - angle * angle / 48 : std::sin(angle / 2) / angle;
    return {std::cos(angle / 2), scale * a};
  }

  // Apply before first, then after.
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

// Fourth-order commutator-free Magnus integration of S' = omega(x) cross S.
// DerivFunc is a caller-owned callable returning double3, independent of S.
// Tolerances control the endpoint step-doubling estimate, NOT an independent
// dense error estimate. Fourth-order accuracy requires smooth coefficients.
// Only the latest accepted step is cached; dense queries never call derivs.
struct StepperMagnus {
  double x_old = 0, x_new = 0, x_max = 0;
  double h_old = 0, h_new = 0;
  // Integral of |omega| over the latest accepted step, in radians.
  // Estimated from the accepted Gauss samples; not the net rotation angle.
  double phase_advance = 0;
  Stokes y_old{}, y_new{};
  double atol, rtol;

  StepperMagnus(double atol = 1e-6, double rtol = 1e-6)
      : atol(atol), rtol(rtol) {
    validate_tolerances();
  }

  void
  init(double x_init, double x_end, double h_init, Stokes const& y_init) {
    if (!(std::isfinite(h_init) && h_init > 0))
      throw std::invalid_argument("h_init must be positive and finite");
    if (!std::isfinite(x_init) || !std::isfinite(x_end) || x_end < x_init ||
        !std::isfinite(x_end - x_init) || !finite(y_init) ||
        !std::isfinite(norm(y_init)))
      throw std::invalid_argument(
          "Invalid interval or Stokes in StepperMagnus");
    x_old = x_new = x_init;
    x_max = x_end;
    h_old = h_new = h_init;
    y_old = y_new = y_init;
    phase_advance = 0;
    prepared = false;
  }

  // Inspect/integrate optical depth using interp before update_old().
  // false means the initialized interval has been completely committed.
  template <class DerivFunc>
  bool
  do_step(DerivFunc const& derivs) {
    if (prepared)
      throw std::logic_error(
          "StepperMagnus requires update_old() before stepping");
    if (x_old == x_max)
      return false;

    double requested = std::min({h_old, x_max - x_old});
    while (true) {
      double end = requested == x_max - x_old
                       ? x_max
                       : std::min(x_old + requested, x_max);
      // Addition can round upward: respect the requested maximum even then.
      if (end - x_old > requested)
        end = std::nextafter(end, x_old);
      double mid = std::midpoint(x_old, end);
      if (!(x_old < mid && mid < end))
        throw std::runtime_error("Stepsize underflow in StepperMagnus");
      double h = end - x_old;
      Segment whole = sample(derivs, x_old, end);
      Segment first = sample(derivs, x_old, mid),
              second = sample(derivs, mid, end);
      first.before = y_old;
      first.after = second.before = first.rotation.apply(y_old);
      second.after = second.rotation.apply(first.after);
      auto coarse = whole.rotation.apply(y_old);
      auto delta = (1.0 / 15) * (second.after + (-1.0) * coarse);
      double scale = atol + rtol * std::max(norm(y_old), norm(second.after));
      double absolute_error = norm(delta);
      // Pure relative tolerance is meaningful for S=0: all rotations fix it.
      double err = absolute_error == 0 ? 0 : absolute_error / scale;
      if (!finite(second.after) || !finite(delta) || !std::isfinite(scale) ||
          !std::isfinite(err))
        throw std::runtime_error("Non-finite error in StepperMagnus");
      double factor =
          err == 0 ? 3 : std::clamp(0.9 * std::pow(err, -0.2), 0.2, 3.0);
      if (err > 1) {
        MC_COUNT(magnus_reject);
        requested = h * factor;
        continue;
      }

      // Reuse the four accepted Gauss samples; only three new endpoint
      // evaluations are needed because the two halves share their midpoint.
      first.omega[0] = coefficient(derivs, x_old);
      first.omega[3] = second.omega[0] = coefficient(derivs, mid);
      second.omega[3] = coefficient(derivs, end);
      segments = {first, second};
      MC_COUNT(magnus_accept);
      x_new = end;
      h_old = h;
      h_new = h * factor;
      y_new = second.after;
      phase_advance = (0.5 * (mid - x_old)) * norm(first.omega[1]) +
                      (0.5 * (mid - x_old)) * norm(first.omega[2]) +
                      (0.5 * (end - mid)) * norm(second.omega[1]) +
                      (0.5 * (end - mid)) * norm(second.omega[2]);
      prepared = true;
      return true;
    }
  }

  void
  update_old() {
    if (!prepared)
      throw std::logic_error("StepperMagnus has no accepted step to commit");
    x_old = x_new;
    y_old = y_new;
    h_old = h_new;
    prepared = false;
  }

  Stokes
  interp(double x) const {
    if (!prepared)
      throw std::logic_error(
          "StepperMagnus dense output requires an accepted step");
    if (!(x >= x_old && x <= x_new))
      throw std::out_of_range("StepperMagnus dense output outside step");
    return (x <= segments[0].end ? segments[0] : segments[1]).at(x);
  }

private:
  static constexpr double g = 0.28867513459481288225;
  static constexpr std::array<double, 4> nodes{0, 0.5 - g, 0.5 + g, 1};

  static
  bool
  finite(double3 s) {
    return std::isfinite(s.x) && std::isfinite(s.y) && std::isfinite(s.z);
  }

  static
  double
  norm(double3 s) {
    return std::hypot(s.x, s.y, s.z);
  }

  void
  validate_tolerances() const {
    if (!std::isfinite(atol) || !std::isfinite(rtol) || atol < 0 || rtol < 0 ||
        (atol == 0 && rtol == 0))
      throw std::invalid_argument("Invalid tolerances in StepperMagnus");
  }

  static
  Rotation
  magnus_rotation(double3 w1, double3 w2, double h) {
    constexpr double a = 0.53867513459481288225, b = -0.03867513459481288225;
    auto first = Rotation::exponential(h * (a * w1 + b * w2));
    auto second = Rotation::exponential(h * (b * w1 + a * w2));
    return Rotation::compose(second, first);
  }

  struct Segment {
    double begin = 0, end = 0;
    Stokes before{}, after{};
    std::array<double3, 4> omega{};
    Rotation rotation;

    double3
    interpolate(double t) const {
      double3 result{};
      for (std::size_t i = 0; i < nodes.size(); ++i) {
        if (t == nodes[i])
          return omega[i];
        double weight = 1;
        for (std::size_t j = 0; j < nodes.size(); ++j)
          if (j != i)
            weight *= (t - nodes[j]) / (nodes[i] - nodes[j]);
        result = result + weight * omega[i];
      }
      return result;
    }

    Stokes
    at(double x) const {
      if (x == begin)
        return before;
      if (x == end)
        return after;
      double t = (x - begin) / (end - begin);
      return magnus_rotation(interpolate(t * nodes[1]),
                             interpolate(t * nodes[2]), x - begin)
          .apply(before);
    }
  };

  std::array<Segment, 2> segments{};
  bool prepared = false;

  template <class DerivFunc>
  double3
  coefficient(DerivFunc const& derivs, double x) const {
    auto w = derivs(x);
    if (!finite(w))
      throw std::runtime_error("Non-finite generator in StepperMagnus");
    return w;
  }

  template <class DerivFunc>
  Segment
  sample(DerivFunc const& derivs, double begin, double end) const {
    Segment segment;
    segment.begin = begin;
    segment.end = end;
    for (int i : {1, 2})
      segment.omega[i] = coefficient(derivs, std::lerp(begin, end, nodes[i]));
    segment.rotation =
        magnus_rotation(segment.omega[1], segment.omega[2], end - begin);
    return segment;
  }
};

struct StepperMagnusWithoutDense {
  double x_old = 0, x_new = 0, x_max = 0;
  double h_old = 0, h_new = 0;
  Stokes y_old{}, y_new{};
  double atol, rtol;

  StepperMagnusWithoutDense(double atol = 1e-6, double rtol = 1e-6)
      : atol(atol), rtol(rtol) {}

  void
  init(double x_init, double x_end, double h_init, Stokes const& y_init) {
    x_old = x_new = x_init;
    x_max = x_end;
    h_old = h_new = h_init;
    y_old = y_new = y_init;
  }

  template <class DerivFunc>
  bool
  do_step(DerivFunc const& derivs) {
    if (x_old == x_max)
      return false;
    double requested = std::min(h_old, x_max - x_old);
    double old_norm = y_old.length();
    while (true) {
      double end = (requested == x_max - x_old) ? x_max : (x_old + requested);
      double mid = std::midpoint(x_old, end);
      if (!(x_old < mid && mid < end))
        throw std::runtime_error("Stepsize underflow in StepperMagnus");

      Rotation whole = sample(derivs, x_old, end);
      Rotation first = sample(derivs, x_old, mid),
               second = sample(derivs, mid, end);
      auto coarse = whole.apply(y_old);
      auto fine = second.apply(first.apply(y_old));
      double delta = (1.0 / 15) * (fine - coarse).length();
      double scale = atol + rtol * std::max(old_norm, fine.length());
      double err = delta / scale;

      double factor =
          err == 0 ? 3 : std::clamp(0.9 * std::pow(err, -0.2), 0.2, 3.0);

      if (err > 1) {
        MC_COUNT(magnus_reject);
        requested *= factor;
        continue;
      }

      MC_COUNT(magnus_accept);
      x_new = end;
      h_old = requested;
      h_new = requested * factor;
      y_new = fine;
      return true;
    }
  }

  void
  update_old() {
    x_old = x_new;
    y_old = y_new;
    h_old = h_new;
  }

private:
  static constexpr double g = 0.28867513459481288225;
  static constexpr std::array<double, 2> nodes{0.5 - g, 0.5 + g};

  static
  Rotation
  magnus_rotation(double3 w1, double3 w2, double h) {
    constexpr double a = 0.53867513459481288225, b = -0.03867513459481288225;
    auto first = Rotation::exponential(h * (a * w1 + b * w2));
    auto second = Rotation::exponential(h * (b * w1 + a * w2));
    return Rotation::compose(second, first);
  }

  template <class DerivFunc>
  Rotation
  sample(DerivFunc const& derivs, double begin, double end) const {
    Rotation rotation =
        magnus_rotation(derivs(std::lerp(begin, end, nodes[0])),
                        derivs(std::lerp(begin, end, nodes[1])), end - begin);
    return rotation;
  }
};
