#include "bfield.hpp"
#include "component.hpp"
#include "find_root.hpp"
#include "geodesic.hpp"
#include "initialize.hpp"
#include "magnus.hpp"
#include "ran.hpp"

template <PhotonInitializer Init, typename Comp>
struct PhotonEvolution {
  struct Options {
    double eta_couple = 1e-3;
  } options;
  BField const& bfield;
  Init const& initializer;
  Comp const& comp;
  Ran ran;
  Geodesic geodesic;

  double3 n, e1, e2;
  double omega_inf;

  bool is_adiabatic;
  AdiabaticPolarization apol;
  Stokes stokes;

  PhotonEvolution(double M_star, double R_star, int seed, BField const& bfield,
                  Init const& init, Comp const& comp)
      : PhotonEvolution(M_star, R_star, seed, bfield, init, comp, Options{}) {}

  PhotonEvolution(double M_star, double R_star, int seed, BField const& bfield,
                  Init const& init, Comp const& comp, Options options)
      : options(options), bfield(bfield), initializer(init), comp(comp),
        ran(seed), geodesic(M_star, R_star) {}

  void
  init() {
    InitializedPhoton photon = initializer.sample_initial_photon(ran);
    geodesic.init(photon.r0, photon.psi0, photon.alpha0);
    n = photon.n;
    e1 = photon.e1;
    e2 = photon.e2;
    omega_inf = photon.omega_inf;
    is_adiabatic = true;
    apol = photon.pol;
  }

  struct Basis {
    double3 r_hat, psi_hat, k_hat, theta_hat, phi_hat;
    double3 s1, s2;
    double3 B;
    double Bperp2, r;
    double c, s; // cos(2chi_B), sin(2chi_B)
    double kappa;
    Basis(double r, double psi, double alpha, double3 n, double3 e1, double3 e2,
          BField const& bfield, double omega_inf, double rs) {
      r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
      psi_hat = cross(n, r_hat);
      k_hat = std::cos(alpha) * r_hat + std::sin(alpha) * psi_hat;
      s1 = cross(n, k_hat);
      s2 = n;
      double mu = std::clamp(r_hat.z, -1.0, 1.0),
             rho = std::hypot(r_hat.x, r_hat.y);
      theta_hat = rho > 0
                      ? double3{r_hat.x * mu / rho, r_hat.y * mu / rho, -rho}
                      : double3{std::copysign(1.0, mu), 0, 0};
      phi_hat = rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0}
                        : double3{0, 1, 0};
      auto [B_vec, _] = bfield.B_and_twist(r, mu);
      B = B_vec.x * r_hat + B_vec.y * theta_hat + B_vec.z * phi_hat;
      double B1 = dot(B, s1), B2 = dot(B, s2);
      Bperp2 = B1 * B1 + B2 * B2;
      this->r = r;
      c = (B1 * B1 - B2 * B2) / Bperp2;
      s = 2 * B1 * B2 / Bperp2;
      kappa = 2.013919e-19 * omega_inf / std::sqrt(1 - rs / r) * Bperp2;
    }
  };

  Basis
  get_basis(double l) const {
    auto [r, psi, alpha] = geodesic.interp(l);
    return Basis(r, psi, alpha, n, e1, e2, bfield, omega_inf, geodesic.rs);
  }

  enum class EvolveResult { Absorbed, Escaped, Scattered };

  struct EscapedData {
    double omega_inf, muk;
    double Q, U, V;
  } escaped_data;

  struct ScatteredData {
    Geodesic::YVector r_psi_alpha;
    double beta;
  } scattered_data;

  EvolveResult
  evolve() {
    while (true) {
      bool absorbed = geodesic.step();
      if (is_adiabatic) {
        auto f_couple = [&](double l) {
          Basis basis = get_basis(l);
          return options.eta_couple * basis.kappa * basis.r - 1;
        };
        if (f_couple(geodesic.l_old()) < 0) {
          double sign = (apol == AdiabaticPolarization::O) ? 1 : -1;
          Basis basis = get_basis(geodesic.l_old());
          stokes = {sign * basis.c, sign * basis.s, 0};
          is_adiabatic = false;
        }
      }

      if (is_adiabatic) {
      } else {
        StepperMagnus stokes_stepper;
        stokes_stepper.init(geodesic.l_old(), geodesic.l_new(),
                            0.1 * geodesic.l_old(), stokes);
        auto omega = [&](double l) {
          Basis basis = get_basis(l);
          return double3{basis.kappa * basis.c, basis.kappa * basis.s, 0};
        };
        while (stokes_stepper.do_step(omega)) {
          // 在 [x_old, x_new] 内积分散射光深。
          // 任意采样点 l：
          //   geometry = geodesic.interp(l);
          //   S = stokes_stepper.dense_out(l);
          // 若发生散射，在此定位、读取状态并退出。
          stokes_stepper.update_old();
        }
        // 保存末态和 h_new，再提交当前测地线步。
      }

      if (absorbed)
        return EvolveResult::Absorbed;
    }
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
  auto powerlaw_electron = CurrentComponent{.dist = PowerLaw2(0.2, 0.8),
                                            .charge_sign = -1,
                                            .current_fraction = 0.5,
                                            .charge_mass_ratio = 1.0};
  auto powerlaw_proton =
      CurrentComponent{.dist = PowerLaw2(0.2, 0.8),
                       .charge_sign = 1,
                       .current_fraction = 0.5,
                       .charge_mass_ratio = 1.0 / 1836.15267343};
  PhotonEvolution photon_evolution(
      M_star, R_star, 42, bfield, initializer,
      std::make_tuple(powerlaw_electron, powerlaw_proton));
  photon_evolution.init();
  photon_evolution.evolve();
}
