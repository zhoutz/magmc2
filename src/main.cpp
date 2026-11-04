#include "bfield.hpp"
#include "component.hpp"
#include "evolution.hpp"
#include "qed_scattering.hpp"
#include <print>
#include <tuple>

int
main(int argc, char** argv) try {
  int N = 10000, first_seed = 42;
  bool trace = false;
  double outer_factor = 1e5, tolerance_scale = 1;
  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "--trace")
      trace = true;
    else if (arg == "--outer-factor" && i + 1 < argc)
      outer_factor = std::stod(argv[++i]);
    else if (arg == "--seed" && i + 1 < argc)
      first_seed = std::stoi(argv[++i]);
    else if (arg == "--tolerance-scale" && i + 1 < argc)
      tolerance_scale = std::stod(argv[++i]);
    else if (i == 1 && !arg.starts_with("--"))
      N = std::stoi(arg);
    else
      throw std::invalid_argument(
          "Usage: main [N] [--trace] [--seed N] [--outer-factor R/Rstar] "
          "[--tolerance-scale factor]");
  }
  if (N < 0 || !std::isfinite(outer_factor) || outer_factor <= 1 ||
      !std::isfinite(tolerance_scale) || tolerance_scale <= 0)
    throw std::invalid_argument(
        "Invalid photon count, outer boundary or tolerance scale");
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

  qed::Options scattering_options;
  scattering_options.width_rtol *= tolerance_scale;
  scattering_options.angular_rtol *= tolerance_scale;
  scattering_options.speed_rtol *= tolerance_scale;
  scattering_options.path_rtol *= tolerance_scale;
  scattering_options.rate_atol *= tolerance_scale;
  QEDScattering resonant(std::make_tuple(powerlaw_electron, powerlaw_positron),
                         scattering_options);
  using Evolution = PhotonEvolution<InitFD11, decltype(resonant)>;
  Evolution::Options options;
  options.geodesic.escape_factor = outer_factor;
  if (trace)
    options.progress = [](double l, double r) {
      std::println(stderr, "Flight: l={} km r={} km", l, r);
    };

  int failed = 0, domain_failed = 0;
  for (int i = 0; i < N; ++i) {
    try {
      Evolution photon(M_star, R_star, first_seed + i, bfield, initializer,
                       resonant, options);
      photon.init();
      int n_scattered = 0;
      while (true) {
        auto result = photon.evolve();
        if (result == EvolveResult::Scattered) {
          ++n_scattered;
          photon.perform_scattering();
        } else if (result == EvolveResult::Escaped) {
          auto const& e = photon.escape_data;
          std::println(
              "Escaped: energy_inf={} muk={} Q={} U={} V={} n_scattered={}",
              e.energy_inf, e.muk, e.Q, e.U, e.V, n_scattered);
          break;
        } else {
          std::println("Absorbed");
          break;
        }
      }
    } catch (std::exception const& error) {
      ++failed;
      auto qed_error = dynamic_cast<qed::Error const*>(&error);
      domain_failed += qed_error && qed_error->kind == qed::Failure::PairDomain;
      std::println(stderr, "Failed: seed={} photon={} category={} {}",
                   first_seed + i, i,
                   qed_error ? qed::failure_name(qed_error->kind) : "transport",
                   error.what());
    }
  }
  std::println("Completed: photons={} failed={} domain_failed={}", N, failed,
               domain_failed);
} catch (std::exception const& error) {
  std::println(stderr, "Simulation setup failed: {}", error.what());
  return 1;
}
