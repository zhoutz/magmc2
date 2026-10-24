#include "component.hpp"
#include "evolution.hpp"
#include <iostream>

void require(bool condition, char const* message) {
  if (!condition) throw std::runtime_error(message);
}

int main() {
  BField field("table/bfield_t10.txt", 1e14, 10);
  require(field.angular_min > 0 && std::isfinite(field.angular_max),
          "Expected finite bounds for the production table");
  for (int i = 0; i <= 200000; ++i) {
    double mu = -1 + 2.0 * i / 200000;
    double magnitude = field.B_and_twist(10, mu).first.length() / (0.5e14);
    require(magnitude >= field.angular_min && magnitude <= field.angular_max,
            "Field magnitude outside bounds");
  }
  auto unsupported = field;
  unsupported.p = 0;
  unsupported.prepare_magnitude_bounds();
  require(unsupported.angular_min == 0 && std::isinf(unsupported.angular_max),
          "Unsupported fields must disable exclusion");

  InitFD11 initializer(10, 0.4);
  Ran random(601028);
  unsigned long long excluded = 0, inward = 0, turning = 0;
  for (int sign : {-1, 1}) {
    auto component = CurrentComponent{PowerLaw2(0.2, std::sqrt(0.75)), sign, 1, 1};
    PhotonEvolution photon(1.4, 10, 42, field, initializer, component);
    for (int ray = 0; ray < 120; ++ray) {
      auto initial = initializer.sample_initial_photon(random);
      initial.r0 = 10 * std::exp(random.U(0, 4));
      initial.alpha0 = random.U(0, pi);
      initial.omega_inf = std::exp(random.U(std::log(0.04), std::log(40.0)));
      photon.init(initial);
      for (int step = 0; step < 160; ++step) {
        bool absorbed = photon.geodesic.step();
        auto [rmin, rmax] = photon.geodesic.radius_bounds();
        bool possible = photon.may_resonate();
        excluded += !possible;
        inward += std::cos(photon.geodesic.stepper.y_old[2]) < 0;
        turning += std::cos(photon.geodesic.stepper.y_old[2]) < 0 &&
                   std::cos(photon.geodesic.stepper.y_new[2]) >= 0;
        for (int node = 0; node <= 24; ++node) {
          double l = std::lerp(photon.geodesic.l_old(), photon.geodesic.l_new(),
                               double(node) / 24);
          auto b = photon.get_basis(l);
          require(b.r >= rmin && b.r <= rmax, "Dense radius outside bounds");
          if (!possible) {
            for (auto mode : {AdiabaticPolarization::O, AdiabaticPolarization::E})
              require(photon.total_rate(photon.rates(b, Polarization(mode))) == 0,
                      "Excluded interval contains a resonant channel");
          }
        }
        if (absorbed || photon.geodesic.stepper.y_new[0] > 1e5) break;
        photon.geodesic.stepper.update_old();
      }
    }
    photon.init();
    photon.geodesic.step();
    auto b = photon.get_basis(photon.geodesic.l_old());
    b.x_e = 1e6;
    int calls = 0;
    photon.rates_at(b, [&] { ++calls; return Polarization{}; });
    require(calls == 0, "Empty channels must not evaluate polarization");
    b.x_e = 1.1;
    b.mu = 0;
    auto lazy = photon.rates_at(b, [&] { ++calls; return Polarization{}; });
    auto eager = photon.rates(b, Polarization{});
    require(calls == 1 && photon.total_rate(lazy) > 0,
            "Active channels must evaluate polarization exactly once");
    for (int i = 0; i < 2; ++i)
      require(lazy[i].weight == eager[i].weight && lazy[i].beta == eager[i].beta,
              "Lazy and eager rates differ");
    photon.options.resonance_bounds = false;
    require(photon.may_resonate(), "Disabled bounds must use full transport");
    photon.options.resonance_bounds = true;
    photon.comp.dist.b_max = 1;
    require(photon.may_resonate(), "Unbounded Lorentz factors need full transport");
  }
  require(excluded > 100 && inward > 100 && turning > 0, "Missing interval coverage");
  std::cout << "Transport bounds passed: excluded=" << excluded
            << " inward_steps=" << inward << " turning_steps=" << turning << '\n';
}
