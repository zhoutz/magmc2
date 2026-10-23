#include "bfield.hpp"
#include "geodesic.hpp"
#include "initialize.hpp"
#include "ran.hpp"

double eta_couple = 1e-3;

struct Basis {
  double3 r_hat, psi_hat, k_hat, theta_hat, phi_hat;
  double3 s1, s2;
  double3 B;
  Basis(double r, double psi, double alpha, double3 n, double3 e1, double3 e2,
        BField const& bfield) {
    r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
    psi_hat = cross(n, r_hat);
    k_hat = std::cos(alpha) * r_hat + std::sin(alpha) * psi_hat;
    s1 = cross(n, k_hat);
    s2 = n;
    double mu = std::clamp(r_hat.z, -1.0, 1.0),
           rho = std::hypot(r_hat.x, r_hat.y);
    theta_hat = rho > 0 ? double3{r_hat.x * mu / rho, r_hat.y * mu / rho, -rho}
                        : double3{std::copysign(1.0, mu), 0, 0};
    phi_hat =
        rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0} : double3{0, 1, 0};
    auto [B_vec, _] = bfield.B_and_twist(r, mu);
    B = B_vec.x * r_hat + B_vec.y * theta_hat + B_vec.z * phi_hat;
  }
};

template <PhotonInitializer Init>
struct PhotonEvolution {
  BField const& bfield;
  Init const& initializer;
  Ran ran;
  Geodesic geodesic;

  double3 n, e1, e2;
  double omega_inf;
  Polarization pol;

  PhotonEvolution(double M_star, double R_star, int seed, BField const& bfield,
                  Init const& init)
      : bfield(bfield), initializer(init), ran(seed), geodesic(M_star, R_star) {
  }

  void
  init() {
    InitializedPhoton photon = initializer.sample_initial_photon(ran);
    geodesic.init(photon.r0, photon.psi0, photon.alpha0);
    n = photon.n;
    e1 = photon.e1;
    e2 = photon.e2;
    omega_inf = photon.omega_inf;
    pol = photon.pol;
  }

  void
  evolve() {
    bool absorbed = geodesic.step();
  }
};

int
main() {
  double M_star = 1.4;
  double R_star = 10;
  double B_pole = 1e14;
  double T_inf = 0.4;
  BField bfield("table/bfield_t10.txt", B_pole, R_star);
  InitFD11 initializer(R_star, T_inf);
  PhotonEvolution photon_evolution(M_star, R_star, 42, bfield, initializer);
  photon_evolution.init();
  photon_evolution.evolve();
}
