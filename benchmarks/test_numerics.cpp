#include "bfield.hpp"
#include "component.hpp"
#include "evolution.hpp"
#include "resonant.hpp"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <print>

int
main() {
  BField field("table/bfield_t10.txt", 1e14, 10);
  Ran random(1701);
  for (int i = 0; i < 100000; ++i) {
    double r = 10 * std::exp(10 * random.U());
    double mu = 2 * random.U() - 1;
    double magnitude = field.evaluate(r, mu).Bmag;
    assert(magnitude >= field.min_strength_at(r));
    assert(magnitude <= field.max_strength_at(r));
  }
  for (double mu : {-1., 0., 1.}) {
    auto b = field.evaluate(10, mu);
    assert(std::isfinite(b.Bmag) && std::isfinite(b.J));
  }
  // The accelerated parser must retain validation and the original numbers.
  {
    std::ifstream input("table/bfield_t10.txt");
    std::string line;
    std::getline(input, line);
    for (int i = 0; i < field.mu_num; ++i) {
      std::getline(input, line);
      std::istringstream original(line);
      double f, fp;
      original >> f >> fp;
      assert(f == field.f[i] && fp == field.fp[i]);
    }
    std::filesystem::create_directories("build/test-fixtures");
    auto filename = "build/test-fixtures/bfield.txt";
    for (std::string row :
         {"nan 0", "1e999 0", "1 0 extra", "+-1 0", "", "1"}) {
      {
        std::ofstream out(filename);
        out << "0 1 0 0 0 1 2\n" << row << "\n0 -2\n";
      }
      bool rejected = false;
      try {
        BField invalid(filename, 1e14, 10);
      } catch (std::runtime_error const&) {
        rejected = true;
      }
      assert(rejected);
    }
    {
      std::ofstream out(filename);
      out << "0 1 0 0 0 1 2\n +1.0 \t -0.0 \n0 -2\n";
    }
    BField valid(filename, 1e14, 10);
    assert(valid.f[0] == 1 && valid.fp[1] == -2);
  }
  CurrentComponent electron{.dist = NegateDist(PowerLaw2(.2, .8)),
                            .charge_sign = -1,
                            .current_fraction = .5,
                            .charge_mass_ratio = 1.};
  CurrentComponent positron{.dist = PowerLaw2(.2, .8),
                            .charge_sign = 1,
                            .current_fraction = .5,
                            .charge_mass_ratio = 1.};
  ResonantScattering scatter(std::make_tuple(electron, positron));
  for (int i = 0; i < 10000; ++i) {
    double beta = 1.6 * random.U() - .8, mu = 2 * random.U() - 1;
    double x = (1 - beta * mu) / std::sqrt((1 - beta) * (1 + beta));
    double D = scatter.discriminant(x, mu);
    if (D <= 0)
      continue;
    auto roots = scatter.resonance_roots(x, mu, D);
    assert(std::min(std::abs(roots[0].beta - beta),
                    std::abs(roots[1].beta - beta)) < 1e-7);
    for (auto root : roots)
      assert(std::abs(x * std::sqrt((1 - root.beta) * (1 + root.beta)) +
                      root.beta * mu - 1) < 1e-10);
    // Compare the cancellation-based specialization with the generic formula.
    double speed = .2 + .6 * random.U();
    double t = (1 - speed) * (1 + speed);
    double expected = positron.dist.inv_abs_b_mean * positron.dist.f(speed) *
                      x * x * t * std::sqrt(t) / std::sqrt(D);
    double fast = positron.dist.resonant_weight(speed, x, std::sqrt(D));
    assert(std::abs(fast - expected) < 1e-12 * expected);
  }
  auto basis = [](double l) {
    Basis b{};
    b.Bmag = (0.8 + 0.2 * l) / B_to_energy;
    b.energy_local = 1;
    b.mu = .5;
    b.J = b.Bmag / (4 * pi * pi);
    b.c = 1;
    b.s = 0;
    return b;
  };
  auto b = basis(0);
  // Shared-root caching must also handle unequal charge/mass ratios and the
  // non-specialized distribution path.
  auto other_mass = positron;
  other_mass.charge_mass_ratio = .5;
  ResonantScattering mixed(std::make_tuple(electron, other_mass));
  ResonantScattering first_only(std::make_tuple(electron));
  ResonantScattering second_only(std::make_tuple(other_mass));
  CurrentComponent thermal{.dist = Boltzmann(.4), .charge_sign = 1,
                           .current_fraction = .5, .charge_mass_ratio = 1};
  ResonantScattering generic(std::make_tuple(thermal));
  for (int i = 0; i < 100; ++i) {
    auto geometry = basis(i/50.);
    auto both = mixed.rates_enveloped(geometry);
    auto first = first_only.rates_enveloped(geometry);
    auto second = second_only.rates_enveloped(geometry);
    for (int j = 0; j < 2; ++j) {
      assert(both[j].weight == first[j].weight);
      assert(both[j+2].weight == second[j].weight);
    }
    for (auto rate : generic.rates_enveloped(geometry))
      assert(std::isfinite(rate.weight) && rate.weight >= 0);
  }
  b.Bmag = std::sqrt(.75) / B_to_energy;
  // Coalescence: the two resonance roots meet at beta=mu.
  auto roots = scatter.resonance_roots(std::sqrt(.75), .5, 0);
  assert(std::abs(roots[0].beta - .5) < 1e-14);
  assert(std::abs(roots[1].beta - .5) < 1e-14);
  b.Bmag *= 1 + 1e-10;
  for (auto rate : scatter.rates_enveloped(b))
    assert(std::isfinite(rate.weight) && rate.weight >= 0);
  assert(positron.dist.f(.2) == 0 && positron.dist.f(.8) == 0);
  assert(!scatter.can_resonate(0, .3));
  assert(!scatter.can_resonate(3.1, 4));
  assert(scatter.can_resonate(.9, 1.1));
  auto zero = electron;
  zero.current_fraction = 0;
  ResonantScattering empty(std::make_tuple(zero));
  Quad quad;
  for (int degree = 0; degree <= 13; ++degree) {
    auto polynomial = [=](double x) { return std::pow(x, degree); };
    auto estimate = gauss_kronrod15(polynomial, 0, 1);
    assert(std::abs(estimate.integral - 1. / (degree + 1)) < 1e-14);
  }
  auto singular = [](double x) { return 1 / std::sqrt(x); };
  assert(gauss_kronrod15(singular, 0, 1).error > 1e-4);
  assert(std::abs(quad.qags(singular, 0, 1, 1e-10, 1e-10) - 2) < 1e-9);
  double tau = -1;
  assert(!empty.next_candidate(0, 1, basis, tau, quad));
  assert(tau == -1);
  // Narrow velocity support must still yield integration cuts/candidates.
  auto narrow = positron;
  narrow.dist = PowerLaw2(.49, .51);
  ResonantScattering narrow_scatter(std::make_tuple(narrow));
  auto cuts = narrow_scatter.integration_cuts(0, 1, basis);
  assert(cuts.size() > 2);
  auto actual_integral = [&](double l) {
    return narrow_scatter.total_rate(narrow_scatter.rates_enveloped(basis(l)));
  };
  double mass = 0;
  for (std::size_t i = 1; i < cuts.size(); ++i)
    mass += quad.qags(actual_integral, cuts[i - 1], cuts[i], 1e-9, 1e-9);
  if (mass > 0) {
    tau = -.5 * mass;
    auto candidate = narrow_scatter.next_candidate(0, 1, basis, tau, quad);
    assert(candidate && *candidate > 0 && *candidate < 1);
  }
  // Known constant-coefficient Stokes rotation, including norm conservation.
  StepperMagnusWithoutDense magnus;
  magnus.init(0, 1, 1, {0, 1, 0});
  auto rotation = [](double) { return double3{0, 0, 2}; };
  while (magnus.do_step(rotation))
    magnus.update_old();
  assert((magnus.y_new - double3{-std::sin(2.), std::cos(2.), 0}).length() <
         1e-12);
  assert(std::abs(magnus.y_new.length() - 1) < 1e-12);
  // A conserved null-geodesic impact parameter tests the relaxed orbit solver.
  for (double alpha : {.1, 1.2, 2.9}) {
    Geodesic orbit(1.4, 10);
    orbit.init(11, 0, alpha);
    auto impact = [&](auto y) {
      return y[0] * std::sin(y[2]) / std::sqrt(1 - orbit.rs / y[0]);
    };
    double initial = impact(orbit.stepper.y_old);
    int event = -1;
    for (int i = 0; i < 1000; ++i) {
      event = orbit.step();
      assert(std::abs(impact(orbit.stepper.y_new) - initial) < 1e-3 * initial);
      if (event >= 0)
        break;
      orbit.stepper.update_old();
    }
    assert(event >= 0);
  }
  InitFD11 initializer(10, .4);
  PhotonEvolution photon(1.4, 10, 42, field, initializer, scatter);
  // Check the interpolation on actual trajectories between validation nodes.
  for (int i = 0; i < 8; ++i) {
    photon.init();
    for (int step = 0; step < 20; ++step) {
      int event = photon.geodesic.step();
      double left = photon.geodesic.l_old(), right = photon.geodesic.l_new();
      auto exact = [&](double l) { return photon.basis(l); };
      TransportCache cache(1e-3);
      cache.reset(left, right, exact);
      for (int j = 0; j <= 31; ++j) {
        double l = std::lerp(left, right, j / 31.);
        auto direct = exact(l);
        auto cached = cache.at(l, exact);
        assert(std::abs(direct.Bmag - cached.Bmag) <= .01 * direct.Bmag);
        assert(std::abs(direct.J - cached.J) <=
               .01 * std::max(direct.J, 1e-100));
        assert(std::abs(direct.mu - cached.mu) <= .01);
      }
      if (event >= 0)
        break;
      photon.geodesic.stepper.update_old();
    }
  }
  int scatter_count = 0;
  for (int i = 0; i < 128; ++i) {
    photon.init();
    bool done = false;
    for (int j = 0; j < 10000; ++j) {
      auto result = photon.evolve();
      if (result == EvolveResult::Scattered) {
        ++scatter_count;
        photon.perform_scattering();
        continue;
      }
      if (result == EvolveResult::Escaped) {
        auto e = photon.escape_data;
        assert(e.energy_inf > 0 && std::isfinite(e.energy_inf));
        assert(std::abs(e.muk) <= 1);
        assert(std::abs(double3{e.Q, e.U, e.V}.length() - 1) < 1e-7);
      }
      done = true;
      break;
    }
    assert(done);
  }
  assert(scatter_count > 0);
  // Pole and inward trajectories exercise the singular coordinate and surface.
  for (double sign : {-1., 1.}) {
    InitializedPhoton p{.n = {1, 0, 0},
                        .e1 = {0, 0, sign},
                        .e2 = {0, -sign, 0},
                        .energy_inf = 1,
                        .r0 = 10.01,
                        .psi0 = 0,
                        .alpha0 = pi,
                        .pol = AdiabaticPolarization::E};
    photon.init(p);
    assert(photon.evolve() == EvolveResult::Absorbed);
  }
  std::println("Numerical checks passed");
}
