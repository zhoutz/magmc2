#include "bfield.hpp"
#include "component.hpp"
#include "evolution.hpp"
#include "resonant.hpp"
#include <charconv>
#include <cstdio>
#include <limits>
#include <print>
#include <string_view>
#include <tuple>

int
main(int argc, char** argv) {
  int N = 1000, seed = 42;
  bool thinning = false, diagnostic = false, no_cache = false;
  double accuracy = 1, step_scale = 1;
  try {
    for (int i = 1; i < argc; ++i) {
      std::string_view arg(argv[i]);
      if (arg == "--diagnostics") {
        diagnostic = true;
        continue;
      }
      if (arg == "--no-cache") {
        no_cache = true;
        continue;
      }
      auto eq = arg.find('=');
      auto key = arg.substr(0, eq);
      std::string_view value;
      if (eq != arg.npos)
        value = arg.substr(eq + 1);
      else if (i + 1 < argc)
        value = argv[++i];
      else
        throw std::invalid_argument("Missing argument value");
      if (key == "--transport") {
        if (value != "reference" && value != "thinning")
          throw std::invalid_argument("Unknown transport mode");
        thinning = value == "thinning";
      } else if (key == "--photons" || key == "--seed") {
        int n;
        auto [end, ec] =
            std::from_chars(value.data(), value.data() + value.size(), n);
        if (ec != std::errc{} || end != value.data() + value.size() || n < 0)
          throw std::invalid_argument("Invalid nonnegative integer");
        (key == "--photons" ? N : seed) = n;
      } else if (key == "--accuracy" || key == "--step-scale") {
        double x;
        auto [end, ec] =
            std::from_chars(value.data(), value.data() + value.size(), x);
        if (ec != std::errc{} || end != value.data() + value.size() ||
            !(x > 0 && x <= 1))
          throw std::invalid_argument("Scale must be in (0,1]");
        (key == "--accuracy" ? accuracy : step_scale) = x;
      } else
        throw std::invalid_argument("Unknown option");
    }
    if (N > 0 && seed > std::numeric_limits<int>::max() - (N - 1))
      throw std::invalid_argument("Seed range overflows int");
  } catch (std::exception const& e) {
    std::println(stderr, "{}", e.what());
    return 2;
  }
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
  resonant.options.quad_atol *= accuracy;
  resonant.options.quad_rtol *= accuracy;

  for (int i = 0; i < N; ++i) {
    PhotonEvolution photon(M_star, R_star, seed + i, bfield, initializer,
                           resonant);
    photon.options.thinning = thinning;
    photon.options.profile = diagnostic;
    photon.options.polarization_tolerance *= accuracy;
    photon.options.cache_tolerance *= accuracy;
    if (no_cache)
      photon.options.cache_tolerance = 0;
    photon.geodesic.max_step_fraction *= step_scale;
    int scatterings = 0;
    photon.init();
    while (true) {
      auto result = photon.evolve();
      if (result == decltype(photon)::EvolveResult::Scattered) {
        // std::println("Scattered: beta={}", photon.scatter_data.beta);
        photon.perform_scattering();
        ++scatterings;
      } else if (result == decltype(photon)::EvolveResult::Escaped) {
        auto const& e = photon.escape_data;
        std::println("Escaped: energy_inf={} muk={} Q={} U={} V={}",
                     e.energy_inf, e.muk, e.Q, e.U, e.V);
        break;
      } else {
        std::println("Absorbed");
        break;
      }
    }
    if (diagnostic)
      std::println(stderr, "{} {} {} {} {} {} {}", seed + i, scatterings,
                   photon.candidates, photon.rejections, photon.magnus_steps,
                   photon.envelope_seconds, photon.polarization_seconds);
  }
}
