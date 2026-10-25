#pragma once

#include "constants.hpp"
#include "doubles.hpp"
#include "hunt.hpp"
#include "interval.hpp"
#include <cmath>
#include <format>
#include <fstream>
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

  // Conservative bounds for the *interpolated* angular field. Filled once.
  struct AngularBound {
    double low, high;
  };
  std::vector<AngularBound> angular_bounds;
  AngularBound global_bound{INFINITY, 0};

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

    // Grow as rows are read, avoiding a huge allocation for a corrupt count.
    for (int i = 0; i < mu_num; i++) {
      std::size_t line_number = static_cast<std::size_t>(i) + 2;
      if (!std::getline(fin, line))
        fail(line_number, "Missing or unreadable table row");
      std::istringstream row(line);
      double fval, fpval;
      if (!(row >> fval >> fpval))
        fail(line_number, "Expected two numeric fields: f and f'");
      row >> std::ws;
      if (!row.eof())
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
    for (int i = 0; i + 1 < mu_num; ++i) {
      double u = std::lerp(mu_min, mu_max, double(i) / (mu_num - 1));
      double v = std::lerp(mu_min, mu_max, double(i + 1) / (mu_num - 1));
      Interval br = Interval(fp[i], fp[i + 1]);
      if (br.lo > br.hi)
        std::swap(br.lo, br.hi);
      Interval fv(std::min(f[i], f[i + 1]), std::max(f[i], f[i + 1]));
      Interval bt = Interval(p) * fv /
                    Interval(std::sqrt(std::max(0.0, 1 - v * v)),
                             std::sqrt(std::max(0.0, 1 - u * u)));
      // At a pole f=0 and sin(theta)=0 simultaneously. Bound the ratio
      // using the linear interpolant f=k*(1-mu) on the last cell.
      if (v == 1 && f[i + 1] == 0) {
        double k = f[i] / (1 - u);
        bt = Interval(
            0, Interval::up(std::abs(p * k) * std::sqrt((1 - u) / (1 + u))));
      }
      if (fv.lo < 0 || A < 0 || !std::isfinite(bt.hi)) {
        angular_bounds.push_back({0, INFINITY});
        continue;
      }
      Interval tw(Interval::down(A * std::pow(fv.lo, 1 / p)),
                  Interval::up(A * std::pow(fv.hi, 1 / p)));
      auto mag = (br.square() + bt.square() + (tw * bt).square()).sqrt();
      // Cover the table's floating-point interpolation and libm rounding too.
      angular_bounds.push_back(
          {std::max(0.0, mag.lo * (1 - 1e-10)), mag.hi * (1 + 1e-10)});
    }
    for (auto const& bound : angular_bounds) {
      global_bound.low = std::min(global_bound.low, bound.low);
      global_bound.high = std::max(global_bound.high, bound.high);
    }
  }

  Interval
  strength_bounds(Interval radius, Interval mu) const {
    mu = mu.abs();
    if (!(radius.lo > 0) || !(R_star > 0) || !(B_pole > 0) || mu.lo < mu_min ||
        mu.hi > mu_max)
      return {0, INFINITY};
    int i =
        std::clamp(int((mu.lo - mu_min) / (mu_max - mu_min) * (mu_num - 1)) - 1,
                   0, mu_num - 2);
    int j =
        std::clamp(int((mu.hi - mu_min) / (mu_max - mu_min) * (mu_num - 1)) + 1,
                   0, mu_num - 2);
    double low = INFINITY, high = 0;
    for (; i <= j; ++i) {
      low = std::min(low, angular_bounds[i].low);
      high = std::max(high, angular_bounds[i].high);
    }
    return Interval(low, high) *
           Interval(Interval::down(0.5 * B_pole *
                                   std::pow(R_star / radius.hi, 2 + p)),
                    Interval::up(0.5 * B_pole *
                                 std::pow(R_star / radius.lo, 2 + p)));
  }

  Interval
  strength_bounds(Interval radius) const {
    if (!(radius.lo > 0) || !(R_star > 0) || !(B_pole > 0))
      return {0, INFINITY};
    return Interval(global_bound.low, global_bound.high) *
           Interval(Interval::down(0.5 * B_pole *
                                   std::pow(R_star / radius.hi, 2 + p)),
                    Interval::up(0.5 * B_pole *
                                 std::pow(R_star / radius.lo, 2 + p)));
  }

  struct BFieldResult {
    double3 B;
    // J denotes J_physical/c, with length expressed in km (J/B is per km).
    double twist, J, Bmag;
  };

  BFieldResult
  evaluate(double r, double mu) const {
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
    double twist = A * std::pow(fval, 1 / p);
    double Bph = twist * Bth;

    double scale = 0.5 * B_pole * std::pow(R_star / r, 2 + p);
    double3 B = double3(scale * Br, scale * Bth, scale * Bph);
    double Bmag = B.length();
    double J = (p + 1) / (4 * pi * r) * twist * Bmag;

    return BFieldResult{B, twist, J, Bmag};
  }
};
