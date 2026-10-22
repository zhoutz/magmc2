#pragma once

#include "bessel_K1_scaled.hpp"
#include <cmath>
#include <concepts>
#include <cstdlib>
#include <stdexcept>

template <class T>
concept ParticleDistribution = requires(T b, double x) {
  { b.b_min } -> std::same_as<double>;
  { b.b_max } -> std::same_as<double>;
  { b.inv_abs_b_mean } -> std::same_as<double>;
  { b.f(x) } -> std::same_as<double>;
};

struct Boltzmann {
  static constexpr double b_min = 0, b_max = 1;
  double b0, exp_k1, a;
  double inv_abs_b_mean;

  Boltzmann(double b0) : b0(b0) {
    if (!std::isfinite(b0) || b0 <= 0 || b0 >= 1) {
      throw std::runtime_error("Invalid b0 value for Boltzmann distribution");
    }
    double s0 = std::sqrt((1 - b0) * (1 + b0));
    a = s0 * (1 + s0) / (b0 * b0);
    exp_k1 = bessel_K1_scaled(a);
    double b_mean = 1 / (a * exp_k1);
    inv_abs_b_mean = 1 / b_mean;
  }

  double
  f(double b) const {
    if (!std::isfinite(b) || b <= b_min || b >= b_max) {
      return 0;
    }
    double s = std::sqrt((1 - b) * (1 + b));
    double gm1 = b * b / (s * (1 + s));
    return std::exp(-a * gm1) / (exp_k1 * s * s * s);
  }
};

struct PowerLaw {
  double b_min, b_max, alpha_pl;
  double inv_abs_b_mean;

  PowerLaw(double b_min, double b_max, double alpha_pl)
      : b_min(b_min), b_max(b_max), alpha_pl(alpha_pl) {
    if (!std::isfinite(b_min) || b_min <= 0 || b_min >= 1) {
      throw std::runtime_error("Invalid b_min value for PowerLaw distribution");
    }
    if (!std::isfinite(b_max) || b_max <= 0 || b_max >= 1) {
      throw std::runtime_error("Invalid b_max value for PowerLaw distribution");
    }
    if (!std::isfinite(alpha_pl) || alpha_pl > 0) {
      throw std::runtime_error(
          "Invalid alpha_pl value for PowerLaw distribution");
    }
    if (b_min >= b_max) {
      throw std::runtime_error(
          "Invalid range: b_min must be less than b_max for PowerLaw "
          "distribution");
    }
  }

  double
  f(double b) const {
    if (!std::isfinite(b) || b <= b_min || b >= b_max) {
      return 0;
    }
  }
};
