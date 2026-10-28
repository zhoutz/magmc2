#include "bfield.hpp"
#include "component.hpp"
#include "evolution.hpp"
#include "resonant.hpp"
#include <print>
#include <tuple>

int
main() {
  double M_star = 1.4;
  double R_star = 10;
  double B_pole = 1e14;
  double T_inf = 0.4;
  BField bfield("table/bfield_t10.txt", B_pole, R_star);
  InitFD11 initializer(R_star, T_inf);

  auto powerlaw_electron =
      CurrentComponent{.dist = NegateDist(PowerLaw2(0.2, 0.8)),
                       .charge_sign = -1,
                       .current_fraction = 0.5,
                       .charge_mass_ratio = 1.0};

  auto powerlaw_positron = CurrentComponent{.dist = PowerLaw2(0.2, 0.8),
                                            .charge_sign = 1,
                                            .current_fraction = 0.5,
                                            .charge_mass_ratio = 1.0};

  ResonantScattering resonant(
      std::make_tuple(powerlaw_electron, powerlaw_positron));

  int N = 10000;
  for (int i = 0; i < N; ++i) {
    PhotonEvolution photon(M_star, R_star, 42 + i, bfield, initializer,
                           resonant);
    photon.init();
    while (true) {
      auto result = photon.evolve();
      if (result == EvolveResult::Scattered) {
        // std::println("Scattered: beta={}", photon.scatter_data.beta);
        photon.perform_scattering();
      } else if (result == EvolveResult::Escaped) {
        auto const& e = photon.escape_data;
        std::println("Escaped: energy_inf={} muk={} Q={} U={} V={}",
                     e.energy_inf, e.muk, e.Q, e.U, e.V);
        break;
      } else {
        std::println("Absorbed");
        break;
      }
    }
  }
}
