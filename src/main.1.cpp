#include "bfield.hpp"
#include "component.hpp"
#include "evolution.hpp"
#include "resonant.hpp"
#include "result.hpp"
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

  int N_group = 10;
  int group_size = 100;
  OutputHDF5 output("data/output.h5");
  for (int i = 0; i < N_group; ++i) {
    int seed = 42 + i;
    ResultPhotonGroup group(seed);
    PhotonEvolution photon(M_star, R_star, seed, bfield, initializer, resonant);
    for (int j = 0; j < group_size; ++j) {
      photon.init();
      int n_scattered = 0;
      while (true) {
        auto result = photon.evolve();
        if (result == EvolveResult::Scattered) {
          ++n_scattered;
          photon.perform_scattering();
        } else if (result == EvolveResult::Escaped) {
          auto const& e = photon.escape_data;
          group.add_escaped(ResultPhoton{.id = j,
                                         .n_scattered = n_scattered,
                                         .energy_inf = e.energy_inf,
                                         .muk = e.muk,
                                         .Q = e.Q,
                                         .U = e.U,
                                         .V = e.V});
          break;
        } else {
          group.add_absorbed();
          break;
        }
      }
    }
    output.write_photon_group("seed_" + std::to_string(seed), group);
  }
}
