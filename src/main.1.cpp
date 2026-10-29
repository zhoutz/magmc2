#include "bfield.hpp"
#include "component.hpp"
#include "evolution.hpp"
#include "resonant.hpp"
#include "result.hpp"
#include <algorithm>
#include <atomic>
#include <exception>
#include <format>
#include <fstream>
#include <iostream>
#include <mutex>
#include <optional>
#include <print>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

int
main() try {
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

  int start_seed = 42;
  int N_group = 1536;
  int group_size = 100000;
  int n_threads = 23;
  std::string error_log_file = "data/error.log";

  OutputHDF5 output("data/output.h5");
  std::ofstream error_log;
  error_log.exceptions(std::ios::failbit | std::ios::badbit);
  error_log.open(error_log_file, std::ios::out | std::ios::app);

  std::vector<int> pending_seeds;
  for (int i = 0; i < N_group; ++i) {
    int seed = start_seed + i;
    // Existing datasets are assumed complete, including failed photons.
    if (!output.file.exist(std::format("seed_{}", seed)))
      pending_seeds.push_back(seed);
  }
  if (pending_seeds.empty()) {
    std::println("All photon groups already exist in output file");
    return 0;
  }

  std::mutex output_mutex, log_mutex, exception_mutex;
  std::atomic<std::size_t> next_task{0};
  std::atomic<bool> stop{false};
  std::exception_ptr first_exception;
  auto record_failure = [&] {
    std::lock_guard lock(exception_mutex);
    if (!first_exception)
      first_exception = std::current_exception();
    stop.store(true);
  };

  auto run_group = [&](int seed) {
    ResultPhotonGroup group(seed);
    PhotonEvolution photon(M_star, R_star, seed, bfield, initializer, resonant);
    for (int j = 0; j < group_size; ++j) {
      std::optional<ResultPhoton> escaped;
      try {
        photon.init();
        int n_scattered = 0;
        while (true) {
          auto result = photon.evolve();
          if (result == EvolveResult::Scattered) {
            ++n_scattered;
            photon.perform_scattering();
          } else if (result == EvolveResult::Escaped) {
            auto const& e = photon.escape_data;
            escaped = ResultPhoton{.id = j,
                                   .n_scattered = n_scattered,
                                   .energy_inf = e.energy_inf,
                                   .muk = e.muk,
                                   .Q = e.Q,
                                   .U = e.U,
                                   .V = e.V};
            break;
          } else {
            break;
          }
        }
      } catch (std::exception const& error) {
        group.add_failed();
        std::string message = error.what();
        std::lock_guard lock(log_mutex);
        error_log << seed << '\t' << j << '\t' << message << '\n';
        error_log.flush();
        continue;
      }
      // Storage failures are task failures, not photon failures.
      if (escaped)
        group.add_escaped(*escaped);
      else
        group.add_absorbed();
    }
    std::lock_guard lock(output_mutex);
    output.write_photon_group(std::format("seed_{}", seed), group);
    output.file.flush();
    std::println("Finished processing seed {}", seed);
  };

  std::vector<std::jthread> workers;
  auto worker_count =
      std::min(static_cast<std::size_t>(n_threads), pending_seeds.size());
  workers.reserve(worker_count);
  try {
    for (std::size_t t = 0; t < worker_count; ++t) {
      workers.emplace_back([&] {
        try {
          while (!stop.load()) {
            auto index = next_task.fetch_add(1);
            if (index >= pending_seeds.size() || stop.load())
              break;
            run_group(pending_seeds[index]);
          }
        } catch (...) {
          record_failure();
        }
      });
    }
  } catch (...) {
    record_failure();
  }
  for (auto& worker : workers)
    worker.join();
  if (first_exception)
    std::rethrow_exception(first_exception);
  error_log.close();
  return 0;
} catch (std::exception const& error) {
  std::cerr << "Simulation failed: " << error.what() << '\n';
  return 1;
} catch (...) {
  std::cerr << "Simulation failed: unknown exception\n";
  return 1;
}
