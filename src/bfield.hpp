#pragma once
#include "diagnostics.hpp"

#include "constants.hpp"
#include "doubles.hpp"
#include "hunt.hpp"
#include <algorithm>
#include <cctype>
#include <charconv>
#include <cmath>
#include <format>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct BField {
  double B_pole, R_star;
  double Delta_phi, p, A, C;
  double mu_min, mu_max;
  int mu_num;
  std::vector<double> f, fp;
  std::vector<double> tabulated_twist;
  std::vector<bool> interpolate_twist;
  double angular_upper_bound = std::numeric_limits<double>::infinity();
  double angular_lower_bound = 0;

  BField(std::string fname, double B_pole, double R_star)
      : B_pole(B_pole), R_star(R_star) {
    std::ifstream fin(fname);
    if (!fin) {
      throw std::runtime_error(
          std::format("Bfield::Bfield(): Cannot open file {}", fname));
    }
    auto fail = [&](std::size_t line_number, std::string const& reason) {
      throw std::runtime_error(
          std::format("BField: {}:{}: {}", fname, line_number, reason));
    };
    std::string line;
    if (!std::getline(fin, line))
      fail(1, "Cannot read table header");
    std::istringstream header(line);
    if (!(header >> Delta_phi >> p >> A >> C >> mu_min >> mu_max >> mu_num))
      fail(1, "Expected seven header fields");
    header >> std::ws;
    if (!header.eof())
      fail(1, "Unexpected data after header fields");
    if (!std::isfinite(Delta_phi) || !std::isfinite(p) || !std::isfinite(A) ||
        !std::isfinite(C) || !std::isfinite(mu_min) || !std::isfinite(mu_max))
      fail(1, "Header values must be finite");
    if (p <= 0)
      fail(1, "Field exponent p must be positive");
    if (mu_num < 2 || mu_min < 0 || mu_max > 1 || mu_min >= mu_max)
      fail(1,
           "Expected at least two grid points with 0 <= mu_min < mu_max <= 1");

    auto skip_space = [](char const*& cursor, char const* end) {
      while (cursor != end && std::isspace(static_cast<unsigned char>(*cursor)))
        ++cursor;
    };
    auto number = [&](char const*& cursor, char const* end, double& value) {
      skip_space(cursor, end);
      if (cursor == end)
        return false;
      // Unlike iostreams, from_chars does not consume a leading '+'.
      if (*cursor == '+') {
        ++cursor;
        if (cursor == end || *cursor == '+' || *cursor == '-')
          return false;
      }
      auto parsed = std::from_chars(cursor, end, value);
      cursor = parsed.ptr;
      return parsed.ec == std::errc{};
    };
    // Grow as rows are read, avoiding a huge allocation for a corrupt count.
    for (int i = 0; i < mu_num; i++) {
      std::size_t line_number = static_cast<std::size_t>(i) + 2;
      if (!std::getline(fin, line))
        fail(line_number, "Missing or unreadable table row");
      char const* cursor = line.data();
      char const* end = cursor + line.size();
      double fval, fpval;
      if (!number(cursor, end, fval) || !number(cursor, end, fpval))
        fail(line_number, "Expected two numeric fields: f and f'");
      skip_space(cursor, end);
      if (cursor != end)
        fail(line_number, "Unexpected data after f and f'");
      if (!std::isfinite(fval) || !std::isfinite(fpval))
        fail(line_number, "Table values must be finite");
      f.push_back(fval);
      fp.push_back(fpval);
    }
    std::size_t line_number = static_cast<std::size_t>(mu_num) + 1;
    while (std::getline(fin, line)) {
      ++line_number;
      if (line.find_first_not_of(" \t\r\n\f\v") != std::string::npos)
        fail(line_number, "Unexpected data after declared table rows");
    }
    if (fin.bad() || !fin.eof())
      fail(line_number + 1, "Error reading table file");
    for (double value : f)
      tabulated_twist.push_back(A * std::pow(value, 1 / p));
    for (int i = 0; i + 1 < mu_num; ++i) {
      double q = 1 / p, delta = f[i + 1] - f[i];
      double a = std::abs(tabulated_twist[i]),
             b = std::abs(tabulated_twist[i + 1]);
      // Linear interpolation error <= max|twist''| * delta_f^2 / 8.
      double curvature = std::abs(q * (q - 1)) *
                         std::max(a / (f[i] * f[i]), b / (f[i + 1] * f[i + 1]));
      interpolate_twist.push_back(f[i] > 0 && f[i + 1] > 0 &&
                                  curvature * delta * delta / 8 <=
                                      1e-8 * std::min(a, b));
    }
    // A conservative bound for the ORIGINAL piecewise-linear f/f' table.
    // On each cell, extrema of f(mu)/sqrt(1-mu^2) occur at the endpoints
    // or at mu=-slope/intercept. No sampled maximum is used as a bound.
    if (mu_min == 0 && mu_max == 1 && f.back() == 0 &&
        std::all_of(f.begin(), f.end(), [](double x) { return x >= 0; })) {
      double br = 0, bt = 0, fmax = 0;
      for (int i = 0; i < mu_num; ++i) {
        br = std::max(br, std::abs(fp[i]));
        fmax = std::max(fmax, f[i]);
        double mu = double(i) / (mu_num - 1);
        if (mu < 1)
          bt = std::max(bt, p * f[i] / std::sqrt((1 - mu) * (1 + mu)));
        if (i + 1 < mu_num) {
          double hi = double(i + 1) / (mu_num - 1);
          double slope = (f[i + 1] - f[i]) / (hi - mu);
          double intercept = f[i] - slope * mu;
          double stationary = -slope / intercept;
          if (stationary > mu && stationary < hi)
            bt = std::max(bt,
                          p * (intercept + slope * stationary) /
                              std::sqrt((1 - stationary) * (1 + stationary)));
        }
      }
      double twist_max = std::abs(A) * std::pow(fmax, 1 / p);
      angular_upper_bound = std::hypot(br, bt, twist_max * bt) * (1 + 1e-12);
      double minimum = std::numeric_limits<double>::infinity();
      for (int i = 0; i + 1 < mu_num; ++i) {
        double lo = double(i) / (mu_num - 1), hi = double(i + 1) / (mu_num - 1);
        double slope = (f[i + 1] - f[i]) / (hi - lo),
               intercept = f[i] - slope * lo;
        double br_min = fp[i] * fp[i + 1] <= 0
                            ? 0
                            : std::min(std::abs(fp[i]), std::abs(fp[i + 1]));
        double bt_min = std::min(
            p * f[i] / std::sqrt((1 - lo) * (1 + lo)),
            hi == 1 ? 0 : p * f[i + 1] / std::sqrt((1 - hi) * (1 + hi)));
        double critical = -slope / intercept;
        if (critical > lo && critical < hi)
          bt_min =
              std::min(bt_min, p * (intercept + slope * critical) /
                                   std::sqrt((1 - critical) * (1 + critical)));
        double twist_min =
            std::abs(A) * std::pow(std::min(f[i], f[i + 1]), 1 / p);
        minimum =
            std::min(minimum, std::hypot(br_min, bt_min, twist_min * bt_min));
      }
      angular_lower_bound = minimum * (1 - 1e-10);
    }
  }

  double
  max_strength_at(double r) const {
    return 0.5 * std::abs(B_pole) * angular_upper_bound *
           std::pow(R_star / r, 2 + p);
  }
  double
  min_strength_at(double r) const {
    return 0.5 * std::abs(B_pole) * angular_lower_bound *
           std::pow(R_star / r, 2 + p);
  }

  struct BFieldResult {
    double3 B;
    // J denotes J_physical/c, with length expressed in km (J/B is per km).
    double twist, J, Bmag;
  };

  BFieldResult
  evaluate(double r, double mu) const {
    MC_COUNT(field);
    double mu_sign = (mu >= 0) ? 1 : -1;
    mu = std::abs(mu);

    UniformHunt mu_hunt(mu_min, mu_max, mu_num);
    auto [i, a] = mu_hunt(mu);

    double fval = (1 - a) * f[i] + a * f[i + 1];
    double fpval = (1 - a) * fp[i] + a * fp[i + 1];
    fpval *= mu_sign;
    double sth = std::sqrt(1 - mu * mu);

    double Br = -fpval;
    double Bth = (sth == 0) ? 0 : p * fval / sth;
    double twist = interpolate_twist[i] ? (1 - a) * tabulated_twist[i] +
                                              a * tabulated_twist[i + 1]
                                        : A * std::pow(fval, 1 / p);
    double Bph = twist * Bth;

    double scale = 0.5 * B_pole * std::pow(R_star / r, 2 + p);
    double3 B = double3(scale * Br, scale * Bth, scale * Bph);
    double Bmag = B.length();
    double J = (p + 1) / (4 * pi * r) * twist * Bmag;

    return BFieldResult{B, twist, J, Bmag};
  }
};
