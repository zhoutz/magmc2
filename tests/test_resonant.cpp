#include "resonant.hpp"
#include "component.hpp"
#include <cassert>
#include <iostream>

namespace {
// A symmetric distribution makes both resonance branches observable.
struct UniformDistribution {
  double b_min = -0.9, b_max = 0.9, inv_abs_b_mean = 1 / 0.45;
  double f(double beta) const {
    return beta > b_min && beta < b_max ? 1 / 1.8 : 0;
  }
};

// Local inputs deliberately independent of PhotonEvolution.
struct Basis {
  double r = 20, Bmag = 1, mu = 0, twist = 0.5, c = 1, s = 0;
  double3 r_hat{1, 0, 0}, B{0, 0, 1};
};

bool close(double a, double b) {
  return std::abs(a - b) <= 1e-11 * std::max({1.0, std::abs(a), std::abs(b)});
}
} // namespace

int main() {
  BField field("table/bfield_t10.txt", 1e14, 10);
  CurrentComponent component{UniformDistribution{}, -1, 1.0, 1.0};
  ResonantScattering scattering(field, component);
  using Scattering = decltype(scattering);
  Scattering::Options options;
  Basis b;
  Polarization pol(AdiabaticPolarization::E);
  // With rs=0 and Bmag=1, these energies give x=2, 1/2 and 1 exactly.
  double energy = B_to_omega / 2;
  auto rates = scattering.rates(b, pol, energy, 0);
  assert(rates[0].weight > 0 && rates[1].weight > 0);
  assert(close(rates[0].beta, -std::sqrt(0.75)));
  assert(close(rates[1].beta, std::sqrt(0.75)));
  assert(close(rates[0].weight, rates[1].weight));

  scattering.comp.current_fraction = 0;
  assert(Scattering::total_rate(scattering.rates(b, pol, energy, 0)) == 0);
  scattering.comp.current_fraction = 1;
  b.twist = 0;
  assert(Scattering::total_rate(scattering.rates(b, pol, energy, 0)) == 0);
  b.twist = 0.5;
  assert(Scattering::discriminant(0.5, 0) < 0);
  assert(Scattering::total_rate(scattering.rates(b, pol, 2 * B_to_omega, 0)) == 0);
  assert(Scattering::discriminant(1, 0) == 0);
  assert(Scattering::total_rate(scattering.rates(b, pol, B_to_omega, 0)) == 0);
  auto near = scattering.rates(b, pol, B_to_omega / (1 + 1e-8), 0);
  assert(std::isfinite(Scattering::total_rate(near)));
  assert(Scattering::total_rate(near) > 0);

  bool rejected = false;
  Ran ran(42), untouched(42);
  try {
    Scattering::select_rate({}, ran);
  } catch (std::runtime_error const&) {
    rejected = true;
  }
  assert(rejected);
  assert(ran.uint64() == untouched.uint64());

  Quad quad;
  auto basis_at = [&](double) { return b; };
  auto pol_at = [&](double) { return pol; };
  double rate = Scattering::total_rate(rates);
  double tau = -2 * rate;
  auto event = scattering.integrate_optical_depth(
      0, 1, basis_at, pol_at, energy, 0, tau, ran, quad, options);
  assert(!event && close(tau, -rate));
  assert(ran.uint64() == untouched.uint64());

  // Constant rate: the event must occur at l=3/4, in the second segment.
  // Integration preserves that segment's starting threshold until reinit.
  pol.set_stokes({0, 0, 1});
  rate = Scattering::total_rate(scattering.rates(b, pol, energy, 0));
  tau = -0.75 * rate;
  event = scattering.integrate_optical_depth(
      0, 1, basis_at, pol_at, energy, 0, tau, ran, quad, options);
  assert(event && close(event->l, 0.75));
  assert(close(tau, -0.25 * rate));
  assert(!event->incoming.is_adiabatic && event->incoming.stokes.z == 1);
  assert(close(std::abs(event->beta), std::sqrt(0.75)));

  for (double beta : {-0.8, 0.0, 0.8}) {
    for (int i = 0; i < 100; ++i) {
      auto photon = scattering.perform_scattering(b, beta, energy, ran);
      assert(std::isfinite(photon.omega_inf) && photon.omega_inf > 0);
      assert(photon.r0 == b.r && photon.psi0 == 0);
      assert(std::isfinite(photon.alpha0));
      assert(close(photon.n.length(), 1) && close(photon.e1.length(), 1));
      assert(close(photon.e2.length(), 1));
      assert(close(dot(photon.n, photon.e1), 0));
      assert(close(dot(cross(photon.e1, photon.e2), photon.n), 1));
      double3 k = std::cos(photon.alpha0) * photon.e1 +
                  std::sin(photon.alpha0) * photon.e2;
      assert(close(k.length(), 1));
      // Elastic scattering in the particle rest frame.
      assert(close(photon.omega_inf * (1 - beta * dot(k, b.B)),
                   energy * (1 - beta * b.mu)));
      assert(photon.pol == AdiabaticPolarization::E ||
             photon.pol == AdiabaticPolarization::O);
    }
  }
  std::cout << "Resonant scattering checks passed\n";
}
