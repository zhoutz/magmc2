#pragma once
#include <algorithm>
#include <cmath>
#include <limits>

// Outward-rounded interval operations used only to exclude scattering.
struct Interval {
  double lo, hi;
  Interval(double value) : lo(value), hi(value) {}
  Interval(double low, double high) : lo(low), hi(high) {}
  static
  double
  down(double x) {
    return std::nextafter(x, -std::numeric_limits<double>::infinity());
  }
  static
  double
  up(double x) {
    return std::nextafter(x, std::numeric_limits<double>::infinity());
  }
  friend
  Interval
  operator+(Interval a, Interval b) {
    return {down(a.lo + b.lo), up(a.hi + b.hi)};
  }
  friend
  Interval
  operator-(Interval a, Interval b) {
    return {down(a.lo - b.hi), up(a.hi - b.lo)};
  }
  friend
  Interval
  operator*(Interval a, Interval b) {
    double v[] = {a.lo * b.lo, a.lo * b.hi, a.hi * b.lo, a.hi * b.hi};
    for (double x : v)
      if (std::isnan(x))
        return {-INFINITY, INFINITY};
    return {down(*std::min_element(v, v + 4)), up(*std::max_element(v, v + 4))};
  }
  friend
  Interval
  operator/(Interval a, Interval b) {
    if (b.lo <= 0 && b.hi >= 0)
      return {-INFINITY, INFINITY};
    return a * Interval(down(1 / b.hi), up(1 / b.lo));
  }
  Interval
  abs() const {
    return {lo <= 0 && hi >= 0 ? 0 : std::min(std::abs(lo), std::abs(hi)),
            std::max(std::abs(lo), std::abs(hi))};
  }
  Interval
  square() const {
    auto a = abs();
    return {std::max(0.0, down(a.lo * a.lo)), up(a.hi * a.hi)};
  }
  Interval
  sqrt() const {
    return {std::max(0.0, down(std::sqrt(std::max(0.0, lo)))),
            up(std::sqrt(std::max(0.0, hi)))};
  }
  Interval
  sin() const {
    if (!std::isfinite(lo) || !std::isfinite(hi) ||
        hi - lo >= 2 * std::acos(-1.0))
      return {-1, 1};
    double low = std::min(std::sin(lo), std::sin(hi));
    double high = std::max(std::sin(lo), std::sin(hi));
    double pi = std::acos(-1.0);
    if (std::ceil((lo - pi / 2) / (2 * pi)) <=
        std::floor((hi - pi / 2) / (2 * pi)))
      high = 1;
    if (std::ceil((lo + pi / 2) / (2 * pi)) <=
        std::floor((hi + pi / 2) / (2 * pi)))
      low = -1;
    double slack = 16 * std::numeric_limits<double>::epsilon() *
                   (1 + std::max(std::abs(lo), std::abs(hi)));
    return {std::max(-1.0, down(low - slack)), std::min(1.0, up(high + slack))};
  }
  Interval
  cos() const {
    return (*this + Interval(std::acos(-1.0) / 2)).sin();
  }
};
