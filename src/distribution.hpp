#pragma once

#include "bessel_K1_scaled.hpp"
#include "quad.hpp"
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

    // In log momentum, dP/d(log u) = u^(alpha_pl + 1) / N_PL.
    // Reference the endpoint with the largest weight to avoid overflowing
    // either the unnormalized power or N_PL. expm1 retains the alpha -> -1
    // limit, where N_PL = log(u_max / u_min).
    double p = alpha_pl + 1;
    double log_range = log_momentum_ratio(b_max, b_min);
    b_ref = p > 0 ? b_max : b_min;
    double log_weight =
        p == 0 ? -std::log(log_range)
               : std::log(std::abs(p)) -
                     std::log(-std::expm1(-std::abs(p) * log_range));
    log_weight_ref = log_weight;

    // Evaluate <beta> = integral beta(u) g(u) du in log momentum.
    // For very steep indices, resolve the concentration near u_min. Beyond
    // |p| log(u/u_min) = 80, the omitted first moment is less than 1e-17
    // of the mean when p <= -2, since beta(u)/beta(u_min) <= u/u_min.
    double span = log_range;
    if (p <= -2 && span > 80 / std::abs(p)) {
      span = 80 / std::abs(p);
    }
    double log_u_ref =
        std::log(b_ref) - 0.5 * (std::log1p(-b_ref) + std::log1p(b_ref));
    double weight = std::exp(log_weight + std::log(span));
    auto first_moment = [&](double x) {
      double offset = span * x;
      double u = std::exp(log_u_ref + (p > 0 ? -offset : offset));
      return (u / std::hypot(1.0, u)) * std::exp(-std::abs(p) * offset) *
             weight;
    };
    Quad quad;
    double b_mean = quad.qags(first_moment, 0, 1, 0, 1e-11);
    if (!std::isfinite(b_mean) || b_mean <= 0) {
      throw std::runtime_error("Invalid mean speed for PowerLaw distribution");
    }
    inv_abs_b_mean = 1 / b_mean;
  }

  double
  f(double b) const {
    if (!std::isfinite(b) || b <= b_min || b >= b_max) {
      return 0;
    }
    double log_ratio = b >= b_ref ? log_momentum_ratio(b, b_ref)
                                  : -log_momentum_ratio(b_ref, b);
    // f(beta) = gamma^3 u^alpha_pl / N_PL, with gamma^3/u
    // = 1/[beta (1-beta^2)]. The population assigns the direction separately.
    return std::exp(log_weight_ref + (alpha_pl + 1) * log_ratio - std::log(b) -
                    std::log1p(-b) - std::log1p(b));
  }

private:
  double b_ref, log_weight_ref;

  static
  double
  log_momentum_ratio(double hi, double lo) {
    // log(u_hi/u_lo), also accurate for nearly coincident cutoffs.
    double delta = hi - lo;
    double log_speed_ratio =
        hi > 2 * lo ? std::log(hi) - std::log(lo) : std::log1p(delta / lo);
    return log_speed_ratio +
           0.5 * (std::log1p(delta / (1 - hi)) - std::log1p(delta / (1 + lo)));
  }
};
