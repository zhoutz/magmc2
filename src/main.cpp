#include "bfield.hpp"
#include "component.hpp"
#include "diagnostics.hpp"
#include "evolution.hpp"
#include "resonant.hpp"
#include <charconv>
#include <chrono>
#include <climits>
#include <cstdio>
#include <print>
#include <stdexcept>
#include <string_view>
#include <tuple>

int
main(int argc, char** argv) try {
  int N = 10000, seed = 42;
  std::string_view output = "text", accuracy_mode = "fast";
  bool summary = false;
  auto integer = [](std::string_view value) {
    int result;
    auto [end, error] =
        std::from_chars(value.data(), value.data() + value.size(), result);
    if (error != std::errc{} || end != value.data() + value.size() ||
        result < 0)
      throw std::invalid_argument("Expected a non-negative integer");
    return result;
  };
  for (int i = 1; i < argc; ++i) {
    std::string_view arg = argv[i];
    if (arg == "--help") {
      std::println("Usage: main [--count N] [--seed N] [--output "
                   "text|csv|none] [--summary] [--accuracy fast|strict]");
      return 0;
    }
    if (arg == "--summary") {
      summary = true;
      continue;
    }
    if (i + 1 == argc)
      throw std::invalid_argument("Missing option value");
    std::string_view value = argv[++i];
    if (arg == "--count")
      N = integer(value);
    else if (arg == "--seed")
      seed = integer(value);
    else if (arg == "--output")
      output = value;
    else if (arg == "--accuracy")
      accuracy_mode = value;
    else
      throw std::invalid_argument("Unknown option");
  }
  if (N == 0 || seed > INT_MAX - (N - 1))
    throw std::invalid_argument("Invalid photon count or seed range");
  if (output != "text" && output != "csv" && output != "none")
    throw std::invalid_argument("Output must be text, csv, or none");
  if (accuracy_mode != "fast" && accuracy_mode != "strict")
    throw std::invalid_argument("Accuracy must be fast or strict");
  auto start = std::chrono::steady_clock::now();
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

  using Scattering = ResonantScattering<decltype(powerlaw_electron),
                                        decltype(powerlaw_positron)>;
  Scattering::Options scattering_accuracy;
  if (accuracy_mode == "strict") {
    scattering_accuracy.n_knots = 4;
    scattering_accuracy.fast_quadrature = false;
    scattering_accuracy.quad_atol = scattering_accuracy.quad_rtol = 1e-8;
    std::fill(bfield.interpolate_twist.begin(), bfield.interpolate_twist.end(),
              false);
  }
  Scattering resonant(std::make_tuple(powerlaw_electron, powerlaw_positron),
                      scattering_accuracy);

  using Evolution = PhotonEvolution<InitFD11, decltype(resonant)>;
  Evolution::Options accuracy;
  if (accuracy_mode == "strict") {
    accuracy.orbit.orbit_atol = accuracy.orbit.orbit_rtol = 1e-10;
    accuracy.orbit.max_step_fraction = .1;
    accuracy.eta_couple = 1e-3;
    accuracy.magnus_atol = accuracy.magnus_rtol = 1e-6;
    accuracy.cache_enabled = false;
  }
  int escaped = 0, absorbed = 0;
  long long scattered = 0;
  if (output == "csv")
    std::println("seed,escaped,n_scattered,energy_inf,muk,Q,U,V");
  auto simulation_start = std::chrono::steady_clock::now();
  PhotonEvolution photon(M_star, R_star, seed, bfield, initializer, resonant,
                         accuracy);
  for (int i = 0; i < N; ++i) {
    photon.ran = Ran(seed + i);
    photon.init();
    int n_scattered = 0;
    while (true) {
      auto result = photon.evolve();
      if (result == EvolveResult::Scattered) {
        ++n_scattered;
        ++scattered;
        photon.perform_scattering();
      } else if (result == EvolveResult::Escaped) {
        auto const& e = photon.escape_data;
        ++escaped;
        if (output == "csv")
          std::println("{},{},{},{},{},{},{},{}", seed + i, 1, n_scattered,
                       e.energy_inf, e.muk, e.Q, e.U, e.V);
        if (output == "text")
          std::println("Escaped: energy_inf={} muk={} Q={} U={} V={}",
                       e.energy_inf, e.muk, e.Q, e.U, e.V);
        break;
      } else {
        ++absorbed;
        if (output == "csv")
          std::println("{},{},{},0,0,0,0,0", seed + i, 0, n_scattered);
        if (output == "text")
          std::println("Absorbed");
        break;
      }
    }
  }
  auto finish = std::chrono::steady_clock::now();
  if (summary)
    std::fprintf(
        stderr,
        "{\"count\":%d,\"seed\":%d,\"escaped\":%d,\"absorbed\":%d,"
        "\"scattered\":%lld,\"total_seconds\":%.9g,\"loop_seconds\":%.9g}\n",
        N, seed, escaped, absorbed, scattered,
        std::chrono::duration<double>(finish - start).count(),
        std::chrono::duration<double>(finish - simulation_start).count());
#ifdef MAGMC_DIAGNOSTICS
  diagnostic::report();
#endif
} catch (std::exception const& error) {
  std::fprintf(stderr, "Simulation failed: %s\n", error.what());
  return 1;
}
