#include "component.hpp"
#include "evolution.hpp"
#include "resonant.hpp"
#include <atomic>
#include <charconv>
#include <chrono>
#include <condition_variable>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <locale>
#include <memory>
#include <mutex>
#include <print>
#include <sstream>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

namespace {
using Count = std::uint64_t;
constexpr double radius = 10, temperature = 0.4, polar_field = 1e14;
constexpr Count block_size = 64;

struct Arguments {
  char setup = 0;
  unsigned threads = 0;
  Count seed = 0, photons = 76800000;
  std::string output;
  bool single = false, untwisted = false, help = false;
};

template <class T>
T
integer(std::string_view value, bool allow_zero = false) {
  T result{};
  auto [end, error] =
      std::from_chars(value.data(), value.data() + value.size(), result);
  if (error != std::errc{} || end != value.data() + value.size() ||
      (!allow_zero && result == 0))
    throw std::invalid_argument("Invalid integer: " + std::string(value));
  return result;
}

Arguments
parse(int argc, char** argv) {
  Arguments args;
  std::vector<std::string_view> seen;
  for (int i = 1; i < argc; ++i) {
    std::string_view key(argv[i]);
    if (std::find(seen.begin(), seen.end(), key) != seen.end())
      throw std::invalid_argument("Repeated argument: " + std::string(key));
    seen.push_back(key);
    if (key == "--help") {
      args.help = true;
      continue;
    }
    if (key == "--single-scattering" || key == "--untwisted") {
      (key == "--untwisted" ? args.untwisted : args.single) = true;
      continue;
    }
    if (key != "--setup" && key != "--threads" && key != "--seed" &&
        key != "--output" && key != "--photons")
      throw std::invalid_argument("Unknown argument: " + std::string(key));
    if (++i == argc)
      throw std::invalid_argument("Missing value for " + std::string(key));
    std::string_view value(argv[i]);
    if (key == "--setup") {
      if (value.size() != 1 || value[0] < 'a' || value[0] > 'j')
        throw std::invalid_argument("--setup must be a single letter a-j");
      args.setup = value[0];
    } else if (key == "--threads") {
      args.threads = integer<unsigned>(value);
    } else if (key == "--seed") {
      args.seed = integer<Count>(value, true);
    } else if (key == "--photons") {
      args.photons = integer<Count>(value);
    } else {
      args.output = value;
    }
  }
  if (args.help)
    return args;
  if (!args.setup || !args.threads || args.output.empty() ||
      std::find(seen.begin(), seen.end(), "--seed") == seen.end())
    throw std::invalid_argument(
        "Required: --setup a-j --threads N --seed N --output FILE");
  if ((args.single || args.untwisted) && args.setup != 'a')
    throw std::invalid_argument("Extra models require --setup a");
  if (args.single && args.untwisted)
    throw std::invalid_argument(
        "--single-scattering and --untwisted are mutually exclusive");
  return args;
}

struct Model {
  char const* name;
  int twist_tenths;
  double gamma_max;
  bool bidirectional = false;
  double e_fraction = 1, cap_degrees = 90, north_fraction = 0.5;
  bool planar = false;
};

constexpr Model models[] = {
    {"t10_g20_c90_Eu", 10, 2},
    {"t10_g20_c90_Eb", 10, 2, true},
    {"t08_g24_c90_Eu", 8, 2.4},
    {"t12_g18_c90_Eu", 12, 1.8},
    {"t14_g17_c90_Eu", 14, 1.65},
    {"t16_g16_c90_Eu", 16, 1.6},
    {"t10_g20_c90_Ou", 10, 2, false, 0},
    {"t10_g20_c5b_Eu", 10, 2, false, 1, 5},
    {"t10_g20_c5s_Eu", 10, 2, false, 1, 5, 0},
    {"t10_g20_c5p_Eu", 10, 2, false, 1, 5, 0.5, true},
};

// SplitMix64's bijective finalizer, with explicitly modulo-2^64 arithmetic.
// A photon's stream depends on its ID, never on scheduling or thread count.
Count
photon_seed(Count seed, Count id) {
  Count z = seed + id + UINT64_C(0x9e3779b97f4a7c15);
  z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
  z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
  return z ^ (z >> 31);
}

using Electrons = CurrentComponent<PowerLaw2>;
using Positrons = CurrentComponent<NegateDist<PowerLaw2>>;
using Physics = ResonantScattering<Electrons, Positrons>;

// Each worker owns its switch. Turning off scattering does not bypass
// geodesic evolution, absorption, or the propagation of Stokes parameters.
struct Scattering {
  using ScatterData = Physics::ScatterData;
  Physics const& physics;
  std::atomic<bool> const& cancelled;
  bool enabled = true;

  template <class BasisAt, class PolarizationAt>
  std::optional<ScatterData>
  integrate_optical_depth(double left, double right, BasisAt const& basis,
                          PolarizationAt const& pol, double& tau, Ran& ran,
                          Quad& quad) const {
    if (cancelled.load(std::memory_order_relaxed))
      throw std::runtime_error("Calculation cancelled after another failure");
    if (!enabled)
      return std::nullopt;
    return physics.integrate_optical_depth(left, right, basis, pol, tau, ran,
                                           quad);
  }

  bool
  resonance_escape(double h, Basis const& old, Basis const& now) const {
    return !enabled || physics.resonance_escape(h, old, now);
  }

  InitializedPhoton
  perform_scattering(ScatterData data, double energy, Ran& ran) const {
    return physics.perform_scattering(data, energy, ran);
  }
};

using Evolution = PhotonEvolution<InitFD11, Scattering>;

struct Worker {
  Scattering scattering;
  Evolution photon;

  Worker(double mass, BField const& field, InitFD11 const& init,
         Physics const& physics, std::atomic<bool> const& cancelled)
      : scattering{physics, cancelled},
        photon(mass, radius, 0, field, init, scattering) {
    if (!photon.quad.workspace)
      throw std::runtime_error("Cannot allocate quadrature workspace");
  }
};

struct Totals {
  Count emitted = 0, escaped = 0, absorbed = 0, scatterings = 0;
  void
  add(Totals const& b) {
    emitted += b.emitted;
    escaped += b.escaped;
    absorbed += b.absorbed;
    scatterings += b.scatterings;
  }
};

int
mode(AdiabaticPolarization pol) {
  return pol == AdiabaticPolarization::E ? 1 : 0;
}

void
simulate(Worker& worker, Arguments const& args, Count id, std::ostream& rows,
         Totals& totals) {
  auto& photon = worker.photon;
  worker.scattering.enabled = !args.untwisted;
  photon.ran.set_seed(photon_seed(args.seed, id));
  photon.init();
  int seed_mode = mode(photon.pol.apol), last_mode = -1;
  Count scatterings = 0;
  while (true) {
    auto result = photon.evolve();
    if (result == Evolution::EvolveResult::Scattered) {
      photon.perform_scattering();
      last_mode = mode(photon.pol.apol);
      ++scatterings;
      if (args.single)
        worker.scattering.enabled = false;
    } else if (result == Evolution::EvolveResult::Escaped) {
      auto const& e = photon.escape_data;
      double norm = std::hypot(e.Q, e.U, e.V);
      if (!(std::isfinite(e.energy_inf) && e.energy_inf > 0 &&
            std::isfinite(e.muk) && std::abs(e.muk) <= 1 &&
            std::isfinite(norm) && std::abs(norm - 1) <= 1e-6))
        throw std::runtime_error("Nonphysical escaped photon");
      rows << id << ' ' << e.energy_inf << ' ' << e.muk << ' ' << e.Q << ' '
           << e.U << ' ' << e.V << ' ' << scatterings << ' ' << seed_mode << ' '
           << last_mode << '\n';
      ++totals.escaped;
      break;
    } else {
      ++totals.absorbed;
      break;
    }
  }
  ++totals.emitted;
  totals.scatterings += scatterings;
}

void
header(std::ostream& out, Arguments const& args, Model const& model,
       BField const& field, std::string const& field_path,
       Worker const& worker) {
  auto const& p = worker.photon;
  auto const& init = p.initializer;
  auto const& physics = worker.scattering.physics;
  StepperMagnus magnus;
  out << std::setprecision(17)
      << "# format = fd11_photons_v1\n# setup = " << args.setup
      << "\n# table1_model = " << model.name << "\n# variant = "
      << (args.single      ? "single_scattering"
          : args.untwisted ? "untwisted"
                           : "standard")
      << "\n# requested = " << args.photons << "\n# seed = " << args.seed
      << "\n# threads = " << args.threads
      << "\n# seed_mapping = splitmix64_finalizer(seed + photon_id + "
         "0x9e3779b97f4a7c15)"
      << "\n# photon_id_origin = 0\n# block_size = " << block_size
      << "\n# R_km = " << radius
      << "\n# M_solar = " << p.geodesic.rs / schwarzschild_radius_of_sun_in_km
      << "\n# rs_km = " << p.geodesic.rs << "\n# spacetime = "
      << (model.planar ? "flat" : "Schwarzschild_R_over_rs_3")
      << "\n# B_pole_G = " << polar_field
      << "\n# T_infinity_keV = " << temperature
      << "\n# field_file = " << field_path
      << "\n# Delta_phi_rad = " << field.Delta_phi
      << "\n# field_p = " << field.p << "\n# field_A = " << field.A
      << "\n# field_C = " << field.C << "\n# gamma_max = " << model.gamma_max
      << "\n# beta_min = 0.2\n# beta_max = "
      << std::get<0>(physics.components).dist.b_max
      << "\n# momentum_exponent = -2\n# electron_charge_sign = -1"
      << "\n# electron_velocity = positive_B\n# electron_current_fraction = "
      << std::get<0>(physics.components).current_fraction
      << "\n# positron_charge_sign = 1\n# positron_velocity = negative_B"
      << "\n# positron_current_fraction = "
      << std::get<1>(physics.components).current_fraction
      << "\n# charge_mass_ratio = 1"
      << "\n# scattering_enabled_initially = " << !args.untwisted
      << "\n# cap_half_angle_rad = " << init.options.cap_half_angle
      << "\n# north_fraction = " << init.options.north_fraction
      << "\n# e_mode_fraction = " << init.options.e_mode_fraction
      << "\n# eta_couple = " << p.options.eta_couple
      << "\n# delta_stokes = " << p.options.delta_stokes
      << "\n# quad_atol = " << physics.options.quad_atol
      << "\n# quad_rtol = " << physics.options.quad_rtol
      << "\n# quad_limit = " << Quad::limit
      << "\n# n_knots = " << physics.options.n_knots
      << "\n# n_scan = " << physics.options.n_scan
      << "\n# n_od_sample = " << physics.options.n_od_sample
      << "\n# orbit_atol = " << p.geodesic.options.orbit_atol
      << "\n# orbit_rtol = " << p.geodesic.options.orbit_rtol
      << "\n# escape_factor = " << p.geodesic.options.escape_factor
      << "\n# magnus_atol = " << magnus.atol
      << "\n# magnus_rtol = " << magnus.rtol
      << "\n# escape_direction_quad_atol = 1e-11"
      << "\n# escape_direction_quad_rtol = 1e-11"
      << "\n# stokes_basis = projected_positive_magnetic_axis_right_handed_sky"
      << "\n# stokes_I_per_photon = 1\n# mode_codes = O:0 E:1 no_scattering:-1"
      << "\n# photon_id energy_inf_keV mu_k Q U V n_scatter seed_mode "
         "last_scatter_mode\n";
}

void
run(Arguments const& args) {
  Model const& model = models[args.setup - 'a'];
  int twist = args.untwisted ? 0 : model.twist_tenths;
  std::string field_path = std::format("table/bfield_t{:02}.txt", twist);
  if (!std::filesystem::exists(field_path))
    throw std::runtime_error(std::format(
        "Missing {}. Run from the repository root; generate it with: "
        "python3 py/gen_bfield.py {}",
        field_path, twist / 10.0));
  BField field(field_path, polar_field, radius);
  if (std::abs(field.Delta_phi - twist / 10.0) > 1e-12 || field.mu_min != 0 ||
      field.mu_max != 1)
    throw std::runtime_error("Magnetic table twist/grid does not match model");
  if (args.untwisted &&
      (field.A != 0 || field.C != 0 || std::abs(field.p - 1) > 1e-12))
    throw std::runtime_error("Untwisted table must have p=1, A=C=0");
  if (std::filesystem::exists(args.output))
    throw std::runtime_error("Output already exists: " + args.output);

  InitFD11 init(radius, temperature,
                {.cap_half_angle = model.cap_degrees * pi / 180,
                 .north_fraction = model.north_fraction,
                 .e_mode_fraction = model.e_fraction});
  PowerLaw2 distribution(
      0.2, std::sqrt(1 - 1 / (model.gamma_max * model.gamma_max)));
  Physics physics(std::make_tuple(
      Electrons{distribution, -1, model.bidirectional ? 0.5 : 1.0, 1.0},
      Positrons{NegateDist(distribution), 1, model.bidirectional ? 0.5 : 0.0,
                1.0}));
  double mass =
      model.planar ? 0 : radius / (3 * schwarzschild_radius_of_sun_in_km);
  std::atomic<bool> cancelled{false};
  std::vector<std::unique_ptr<Worker>> workers;
  // Quad's constructor changes GSL's global error handler. Construct ALL
  // workspaces on the main thread before any worker starts using GSL.
  for (unsigned i = 0; i < args.threads; ++i)
    workers.push_back(
        std::make_unique<Worker>(mass, field, init, physics, cancelled));

  auto parent = std::filesystem::path(args.output).parent_path();
  if (!parent.empty())
    std::filesystem::create_directories(parent);
  std::ofstream out;
  out.exceptions(std::ios::failbit | std::ios::badbit);
  out.imbue(std::locale::classic());
  // C++23 exclusive creation also protects against a concurrent creator or
  // a dangling symlink, beyond the earlier friendly existence diagnostic.
  out.open(args.output, std::ios::out | std::ios::noreplace);
  header(out, args, model, field, field_path, *workers.front());
  out.flush();

  std::mutex mutex;
  std::condition_variable ready;
  std::exception_ptr failure;
  Count next_id = 0, next_write = 0;
  Totals totals;
  std::atomic<Count> completed{0};
  auto start = std::chrono::steady_clock::now(), last_report = start;
  auto report = [&](bool force = false) {
    // Caller holds mutex. Completed may include blocks not yet written.
    auto now = std::chrono::steady_clock::now();
    if (force || now - last_report >= std::chrono::seconds(10)) {
      double seconds = std::chrono::duration<double>(now - start).count();
      auto done = completed.load(std::memory_order_relaxed);
      std::println(stderr,
                   "setup={} completed={}/{} written={} escaped={} absorbed={} "
                   "elapsed={:.1f}s rate={:.1f}/s",
                   args.setup, done, args.photons, totals.emitted,
                   totals.escaped, totals.absorbed, seconds,
                   seconds > 0 ? done / seconds : 0);
      last_report = now;
    }
  };
  auto fail = [&](std::exception_ptr error) {
    {
      std::lock_guard lock(mutex);
      if (!failure)
        failure = error;
      cancelled.store(true, std::memory_order_relaxed);
    }
    ready.notify_all();
  };
  auto work = [&](unsigned index) {
    Count active_id = 0;
    try {
      while (true) {
        Count begin, end;
        {
          std::lock_guard lock(mutex);
          if (cancelled.load(std::memory_order_relaxed) ||
              next_id == args.photons)
            return;
          begin = next_id;
          end = begin + std::min(block_size, args.photons - begin);
          next_id = end;
        }
        std::ostringstream rows;
        rows.exceptions(std::ios::failbit | std::ios::badbit);
        rows.imbue(std::locale::classic());
        rows << std::setprecision(17);
        Totals block;
        for (active_id = begin; active_id < end; ++active_id) {
          if (cancelled.load(std::memory_order_relaxed))
            return;
          simulate(*workers[index], args, active_id, rows, block);
          completed.fetch_add(1, std::memory_order_relaxed);
          std::lock_guard lock(mutex);
          report();
        }
        std::unique_lock lock(mutex);
        ready.wait(lock, [&] {
          return cancelled.load(std::memory_order_relaxed) ||
                 next_write == begin;
        });
        if (cancelled.load(std::memory_order_relaxed))
          return;
        out << rows.view();
        out.flush();
        totals.add(block);
        next_write = end;
        lock.unlock();
        ready.notify_all();
      }
    } catch (std::exception const& error) {
      fail(std::make_exception_ptr(std::runtime_error(
          std::format("worker={} photon_id={} (or block write): {}", index,
                      active_id, error.what()))));
    } catch (...) {
      fail(std::current_exception());
    }
  };

  std::println(stderr,
               "FD11 setup={} variant={} photons={} threads={} output={}",
               args.setup,
               args.single      ? "single_scattering"
               : args.untwisted ? "untwisted"
                                : "standard",
               args.photons, args.threads, args.output);
  // N-1 background workers plus the main thread: no extra writer thread.
  std::vector<std::jthread> threads;
  try {
    for (unsigned i = 1; i < args.threads; ++i)
      threads.emplace_back(work, i);
    work(0);
  } catch (...) {
    fail(std::current_exception());
  }
  for (auto& thread : threads)
    thread.join();
  if (failure)
    std::rethrow_exception(failure);
  if (totals.emitted != args.photons ||
      totals.emitted != totals.escaped + totals.absorbed)
    throw std::runtime_error("Photon accounting mismatch");
  out << "# emitted = " << totals.emitted << "\n# escaped = " << totals.escaped
      << "\n# absorbed = " << totals.absorbed
      << "\n# scatterings = " << totals.scatterings << '\n';
  out.flush();
  report(true);
  auto marker_offset = static_cast<std::uintmax_t>(out.tellp());
  try {
    out << "# complete = true\n";
    out.close();
  } catch (...) {
    // A short write can leave the entire word "true" but lose its newline.
    // Close the buffer before removing even a partially written success marker.
    out.exceptions(std::ios::goodbit);
    if (out.is_open())
      out.close();
    std::error_code error;
    std::filesystem::resize_file(args.output, marker_offset, error);
    if (error)
      throw std::runtime_error(std::format(
          "Output finalization failed; could not remove completion marker: {}. "
          "Treat {} as incomplete.",
          error.message(), args.output));
    throw;
  }
}
} // namespace

int
main(int argc, char** argv) {
  try {
    auto args = parse(argc, argv);
    if (args.help) {
      std::println(
          "Usage: fd11 --setup a-j --threads N --seed N --output FILE\n"
          "            [--photons N] [--single-scattering | --untwisted]\n"
          "--photons defaults to 76800000 emitted photons (including "
          "absorbed).\n"
          "--seed is an unsigned 64-bit integer, including zero.\n"
          "--threads counts all compute threads, including the main thread.\n"
          "Extra variants require --setup a; Table 1 maps top-to-bottom to "
          "a-j.\n"
          "Run from the repository root. Existing outputs are never "
          "overwritten.\n"
          "See docs/fd11.md for model parameters, columns and all twelve "
          "runs.");
      return 0;
    }
    run(args);
    return 0;
  } catch (std::exception const& error) {
    std::println(stderr, "fd11 failed: {}", error.what());
    return 1;
  }
}
