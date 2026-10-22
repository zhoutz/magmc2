#include "bfield.hpp"
#include "geodesic.hpp"
#include "initialize.hpp"
#include "ran.hpp"

int
main() {
  double M_star = 1.4;
  double R_star = 10;
  double B_pole = 1e14;
  double T_inf = 0.4;
  BField bfield("table/bfield_t10.txt", B_pole, R_star);
  Ran ran(42);
  Geodesic geodesic(M_star, R_star);
  for (int i = 0; i < 10; ++i) {
    auto photon = sample_initial_photon_fd11(ran, R_star, T_inf);
    geodesic.init(photon.r0, photon.psi0, photon.alpha0);
  }
}

