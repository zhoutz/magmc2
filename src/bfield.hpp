#pragma once

#include "constants.hpp"
#include "doubles.hpp"
#include "hunt.hpp"
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
