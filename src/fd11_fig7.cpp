#include "bfield.hpp"
#include "component.hpp"
#include "evolution.hpp"
#include "photon_stream.hpp"
#include <charconv>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <omp.h>
#include <print>
#include <string_view>

// PhotonEvolution evaluates dist.f(charge_sign * beta). FD11's electrons
// travel along +B (north to south), independently of their negative charge.
// Reflect only the distribution's argument, not the field or the sky axes.
struct FD11Electrons {
  PowerLaw2 positive{0.2, std::sqrt(0.75)};
  double b_min = -positive.b_max, b_max = -positive.b_min;
  double inv_abs_b_mean = positive.inv_abs_b_mean;
  double f(double x) const { return positive.f(-x); }
};

namespace {
struct Arguments {
  unsigned long long photons = 1000000;
  int seed = 110701;
  std::string output = "data/fd11_fig7_photons.txt";
  std::string field = "table/bfield_t10.txt";
  bool tight = false;
  bool fast = false, photon_streams = false;
  int threads = 1;
};

template <class T>
T integer(std::string_view text) {
  T value{};
  auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
  if (error != std::errc{} || end != text.data() + text.size() || value <= 0)
    throw std::invalid_argument("Expected a positive integer: " + std::string(text));
  return value;
}

Arguments parse(int argc, char** argv) {
  Arguments args;
  for (int i = 1; i < argc; ++i) {
    std::string_view key(argv[i]);
    if (key == "--help") {
      std::println("Usage: fd11_fig7 [--photons N] [--seed N] [--output FILE]"
                   " [--field FILE] [--tight] [--fast] [--threads N] [--photon-streams]\n"
                   "N counts emitted photons, including those later absorbed.\n"
                   "--tight tightens coupling, escape and quadrature tolerances.\n"
                   "--fast enables validated local interpolation and cheaper optical-depth settings.\n"
                   "--threads N uses N workers and stable per-photon random streams.\n"
                   "--photon-streams uses the same random streams with one worker for comparisons.\n"
                   "Run from the repository root or supply --field.");
      std::exit(0);
    }
    if (key == "--tight") {
      args.tight = true;
      continue;
    }
    if (key == "--fast") { args.fast = true; continue; }
    if (key == "--photon-streams") { args.photon_streams = true; continue; }
    if (i + 1 == argc)
      throw std::invalid_argument("Missing value for " + std::string(key));
    std::string_view value(argv[++i]);
    if (key == "--photons") args.photons = integer<unsigned long long>(value);
    else if (key == "--seed") args.seed = integer<int>(value);
    else if (key == "--output") args.output = value;
    else if (key == "--field") args.field = value;
    else if (key == "--threads") args.threads = integer<int>(value);
    else throw std::invalid_argument("Unknown argument: " + std::string(key));
  }
  if (args.threads > 256) throw std::invalid_argument("At most 256 workers supported");
  if (args.threads > 1) args.photon_streams = true;
  return args;
}
} // namespace

int main(int argc, char** argv) {
  unsigned long long emitted = 0, escaped = 0, absorbed = 0, scatterings = 0;
  std::ofstream out;
  try {
    auto args = parse(argc, argv);
    constexpr double R = 10, temperature = 0.4, B = 1e14;
    constexpr double mass = R / (3 * schwarzschild_radius_of_sun_in_km);
    BField field(args.field, B, R);
    if (!std::isfinite(field.Delta_phi) || std::abs(field.Delta_phi - 1) > 1e-12)
      throw std::runtime_error("FD11 Figure 7 requires Delta_phi = 1 rad");
    InitFD11 initializer(R, temperature);
    CurrentComponent electrons{.dist = FD11Electrons{}, .charge_sign = -1,
                               .current_fraction = 1.0, .charge_mass_ratio = 1.0};
    using Evolution = PhotonEvolution<InitFD11, decltype(electrons)>;
    Evolution::Options options;
    if (args.fast) {
      options.quad_atol = 1e-6;
      options.quad_rtol = 1e-5;
      options.n_knots = 2;
      options.cache_tolerance = 1e-7;
      options.event_rtol = 1e-7;
    }
    if (args.tight) {
      options.eta_couple = 3e-4;
      options.escape_tolerance = 3e-4;
      options.quad_atol = options.quad_rtol = 1e-9;
      options.n_scan = 16;
      if (args.fast) {
        options.quad_atol = 1e-8;
        options.quad_rtol = 1e-7;
        options.n_knots = 8;
        options.cache_tolerance = 1e-9;
        options.event_rtol = 1e-9;
      }
    }
    // Construct on the main thread: Quad sets GSL's global error handler.
    // Each worker owns its RNG, integrators and mutable interpolation cache.
    std::vector<std::unique_ptr<Evolution>> workers;
    for (int i = 0; i < args.threads; ++i)
      workers.push_back(std::make_unique<Evolution>(
          mass, R, args.seed, field, initializer, electrons, options));
    for (double beta : workers.front()->knots)
      if (!(beta >= 0.2 && beta <= std::sqrt(0.75)))
        throw std::runtime_error("Incorrect electron velocity direction");

    if (std::filesystem::exists(args.output))
      throw std::runtime_error("Output already exists: " + args.output);
    auto parent = std::filesystem::path(args.output).parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    out.exceptions(std::ios::failbit | std::ios::badbit);
    out.open(args.output);
    out << std::setprecision(17)
        << "# format = fd11_fig7_photons_v1\n"
        << "# model = t10_g20_c90_Eu\n"
        << "# requested = " << args.photons << "\n# seed = " << args.seed
        << "\n# R_km = " << R << "\n# M_solar = " << mass
        << "\n# R_over_rs = 3\n# B_pole_G = " << B
        << "\n# T_infinity_keV = " << temperature
        << "\n# Delta_phi_rad = " << field.Delta_phi
        << "\n# field_p = " << field.p << "\n# field_A = " << field.A
        << "\n# beta_min = 0.2\n# beta_max = " << std::sqrt(0.75)
        << "\n# gamma_max = 2\n# momentum_exponent = -2"
        << "\n# charge_sign = -1\n# velocity_direction = along_positive_B"
        << "\n# current_fraction = 1\n# charge_mass_ratio = 1"
        << "\n# cap_half_angle_rad = " << initializer.options.cap_half_angle
        << "\n# north_fraction = 0.5\n# e_mode_fraction = 1"
        << "\n# eta_couple = " << options.eta_couple
        << "\n# escape_tolerance = " << options.escape_tolerance
        << "\n# quad_atol = " << options.quad_atol
        << "\n# quad_rtol = " << options.quad_rtol
        << "\n# n_knots = " << options.n_knots << "\n# n_scan = " << options.n_scan
        << "\n# geodesic_atol = 1e-10\n# geodesic_rtol = 1e-10"
        << "\n# magnus_atol = 1e-6\n# magnus_rtol = 1e-6"
        << "\n# stokes_basis = projected_magnetic_axis_right_handed_sky"
        << "\n# photon_id energy_keV mu Q U V n_scatter\n";
    if (args.fast || args.photon_streams)
      out << "# numeric_mode = " << (args.fast ? "fast_v1" : "reference")
          << "\n# random_stream = " << (args.photon_streams ? "photon_id_v1" : "legacy")
          << "\n# cache_tolerance = " << options.cache_tolerance
          << "\n# event_rtol = " << options.event_rtol << '\n';
    out.flush();
    auto start = std::chrono::steady_clock::now(), last_report = start;
    struct Record {
      Evolution::EscapedData data{};
      unsigned long long scatterings = 0;
      bool escaped = false;
    };
    constexpr unsigned long long block_size = 2048;
    std::vector<Record> records(block_size);
    std::vector<std::exception_ptr> errors(args.threads);
    for (unsigned long long begin = 0; begin < args.photons;) {
      int count = static_cast<int>(std::min(block_size, args.photons - begin));
      #pragma omp parallel for num_threads(args.threads) schedule(dynamic, 4)
      for (int j = 0; j < count; ++j) {
        int worker = omp_get_thread_num();
        if (errors[worker]) continue;
        try {
          auto& photon = *workers[worker];
          auto& record = records[j];
          record = {};
          if (args.photon_streams)
            photon.ran.set_seed(photon_seed(args.seed, begin + j));
          photon.init();
          while (true) {
            auto result = photon.evolve();
            if (result == Evolution::EvolveResult::Scattered) {
              if (!(photon.scattered_data.beta > 0.2 &&
                    photon.scattered_data.beta < std::sqrt(0.75)))
                throw std::runtime_error("Scattered electron outside FD11 velocity range");
              ++record.scatterings;
              photon.perform_scattering();
            } else if (result == Evolution::EvolveResult::Escaped) {
              auto const& e = photon.escaped_data;
              double norm = std::hypot(e.Q, e.U, e.V);
              if (!(std::isfinite(e.omega_inf) && e.omega_inf > 0 &&
                    std::isfinite(e.muk) && std::abs(e.muk) <= 1 &&
                    std::isfinite(norm) && std::abs(norm - 1) < 1e-7))
                throw std::runtime_error("Nonphysical escaped photon");
              record.data = e;
              record.escaped = true;
              break;
            } else {
              break;
            }
          }
        } catch (...) {
          errors[worker] = std::current_exception();
        }
      }
      for (auto const& error : errors)
        if (error) std::rethrow_exception(error);
      // A single writer emits records in ID order, independent of scheduling.
      for (int j = 0; j < count; ++j) {
        auto const& record = records[j];
        scatterings += record.scatterings;
        if (record.escaped) {
          auto const& e = record.data;
          out << begin + j << ' ' << e.omega_inf << ' ' << e.muk << ' ' << e.Q
              << ' ' << e.U << ' ' << e.V << ' ' << record.scatterings << '\n';
          ++escaped;
        } else {
          ++absorbed;
        }
        ++emitted;
      }
      begin += count;
      auto now = std::chrono::steady_clock::now();
      if (now - last_report > std::chrono::seconds(30) || emitted == args.photons) {
        double seconds = std::chrono::duration<double>(now - start).count();
        std::println(stderr, "seed={} emitted={}/{} escaped={} absorbed={} elapsed={:.1f}s rate={:.1f}/s",
                     args.seed, emitted, args.photons, escaped, absorbed, seconds,
                     emitted / seconds);
        out.flush();
        last_report = now;
      }
    }
    out << "# emitted = " << emitted << "\n# escaped = " << escaped
        << "\n# absorbed = " << absorbed << "\n# scatterings = " << scatterings
        << "\n# complete = true\n";
    out.close();
    return 0;
  } catch (std::exception const& error) {
    std::println(stderr, "fd11_fig7 failed after {} complete photons (escaped={}, absorbed={}): {}",
                 emitted, escaped, absorbed, error.what());
    return 1;
  }
}
