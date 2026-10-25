#pragma once

#include "constants.hpp"
#include "doubles.hpp"
#include "hunt.hpp"
#include <cmath>
#include <format>
#include <fstream>
#include <iostream>
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
    fin >> Delta_phi >> p >> A >> C >> mu_min >> mu_max >> mu_num;
    f.resize(mu_num);
    fp.resize(mu_num);
    for (int i = 0; i < mu_num; i++) {
      fin >> f[i] >> fp[i];
    }
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
