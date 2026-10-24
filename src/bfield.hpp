#pragma once

#include "doubles.hpp"
#include "hunt.hpp"
#include <algorithm>
#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>

struct BField {
  double B_pole, R_star;
  double Delta_phi, p, A, C;
  double mu_min, mu_max;
  int mu_num;
  std::vector<double> f, fp;
  // Bounds on |B| / [0.5*B_pole*(R_star/r)^(2+p)], over all mu.
  double angular_min = 0, angular_max = std::numeric_limits<double>::infinity();

  BField(std::string fname, double B_pole, double R_star)
      : B_pole(B_pole), R_star(R_star) {
    std::ifstream fin(fname);
    if (!fin) {
      throw std::runtime_error(
          std::format("Bfield::Bfield(): Cannot open file {}", fname));
    }
    fin >> Delta_phi >> p >> A >> C >> mu_min >> mu_max >> mu_num;
    f.resize(mu_num);
    fp.resize(mu_num);
    for (int i = 0; i < mu_num; i++) {
      fin >> f[i] >> fp[i];
    }
    prepare_magnitude_bounds();
  }

  void
  prepare_magnitude_bounds() {
    angular_min = 0;
    angular_max = std::numeric_limits<double>::infinity();
    // Unsupported tables retain uninformative bounds and use full transport.
    if (!(p > 0 && A >= 0 && mu_min == 0 && mu_max == 1 && mu_num >= 2) ||
        f.back() != 0)
      return;
    double lo = std::numeric_limits<double>::infinity(), hi = 0;
    for (int i = 0; i + 1 < mu_num; ++i) {
      if (!(std::isfinite(f[i]) && std::isfinite(f[i + 1]) &&
            std::isfinite(fp[i]) && std::isfinite(fp[i + 1]) &&
            f[i] >= 0 && f[i + 1] >= 0))
        return;
      double u = double(i) / (mu_num - 1), v = double(i + 1) / (mu_num - 1);
      double fmin = std::min(f[i], f[i + 1]), fmax = std::max(f[i], f[i + 1]);
      double brmin = fp[i] * fp[i + 1] <= 0
                         ? 0 : std::min(std::abs(fp[i]), std::abs(fp[i + 1]));
      double brmax = std::max(std::abs(fp[i]), std::abs(fp[i + 1]));
      double btmin = p * fmin / std::sqrt((1 - u) * (1 + u));
      // In the last cell f(mu)=f(u)*(1-mu)/(1-u), so f/sin(theta)
      // decreases to zero. Else use separate numerator/denominator bounds.
      double btmax = p * fmax / std::sqrt(
          i + 2 == mu_num ? (1 - u) * (1 + u) : (1 - v) * (1 + v));
      double bpmax = A * std::pow(fmax, 1 / p) * btmax;
      lo = std::min(lo, std::hypot(brmin, btmin));
      hi = std::max(hi, std::hypot(brmax, btmax, bpmax));
    }
    // Pad interval arithmetic and the Cartesian basis roundoff generously.
    angular_min = lo * (1 - 1e-10);
    angular_max = hi * (1 + 1e-10);
  }

  std::pair<double3, double>
  B_and_twist(double r, double mu) const {
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

    return std::make_pair(double3(scale * Br, scale * Bth, scale * Bph), twist);
  }
};
