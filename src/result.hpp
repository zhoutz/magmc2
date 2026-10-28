#pragma once

#include <vector>
// #include "../ext/highfive/include/highfive/highfive.hpp"
#include "highfive/highfive.hpp"

struct ResultPhoton {
  int id, n_scattered;
  double energy_inf, muk;
  double Q, U, V;
};

struct ResultPhotonGroup {
  std::vector<ResultPhoton> photons;
  int n_absorbed, n_escaped;
  int seed;
  ResultPhotonGroup(int seed) : n_absorbed(0), n_escaped(0), seed(seed) {}
  void
  add_absorbed() {
    ++n_absorbed;
  }
  void
  add_escaped(ResultPhoton const& photon) {
    ++n_escaped;
    photons.push_back(photon);
  }
};
