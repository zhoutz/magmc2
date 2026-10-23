#define MAGMC_NO_MAIN
#include "../src/main.cpp"
#include <iostream>
#include <string_view>

namespace {
void
require(bool ok, std::string_view message) {
  if (!ok)
    throw std::runtime_error(std::string(message));
}
void
near(double actual, double expected, double atol, std::string_view message) {
  if (!std::isfinite(actual) || std::abs(actual - expected) > atol)
    throw std::runtime_error(
        std::format("{}: {} != {} (+/- {})", message, actual, expected, atol));
}
void
near_vector(double3 a, double3 b, double atol, std::string_view message) {
  near((a + (-1.0) * b).length(), 0, atol, message);
}
using Component = CurrentComponent<PowerLaw2>;
using Components = std::tuple<Component, Component>;
using Evolution = PhotonEvolution<InitFD11, Components>;
Components
components(double e_fraction = 0.5, double p_fraction = 0.5,
           double p_mass_ratio = 1 / 1836.15267343) {
  return {Component{PowerLaw2(0.2, 0.8), -1, e_fraction, 1},
          Component{PowerLaw2(0.2, 0.8), 1, p_fraction, p_mass_ratio}};
}
InitializedPhoton
radial(double radius, double z, double energy = 1) {
  double3 r{std::sqrt(1 - z * z), 0, z}, n{0, 1, 0};
  return {n, r, cross(n, r), energy, radius, 0, 0, AdiabaticPolarization::E};
}

void
test_rates(BField const& field, InitFD11 const& source) {
  Evolution mixed(1.4, 10, 42, field, source, components(0.3, 0.7, 0.9));
  mixed.init();
  auto b = mixed.basis({100, 0, 0});
  b.mu = 0;
  b.x_e = 1.15;
  b.c = 1;
  b.s = 0;
  b.twist = 0.5;
  Stokes S{-0.3, 0.4, std::sqrt(0.75)};
  auto weights = mixed.rates(b, S, false);
  double expected = 0;
  for (auto const& rate : weights) {
    if (rate.weight == 0)
      continue;
    require(rate.component == 0 ? rate.beta < 0 : rate.beta > 0,
            "Population signed velocity support");
    auto c =
        rate.component == 0 ? std::get<0>(mixed.comp) : std::get<1>(mixed.comp);
    double x = b.x_e * c.charge_mass_ratio, beta = rate.beta;
    near(x, (1 - beta * b.mu) / std::sqrt(1 - beta * beta), 2e-15,
         "Root obeys unsquared resonance equation");
    double mup = (b.mu - beta) / (1 - beta * b.mu);
    double overlap = ((1 + mup * mup) + (mup * mup - 1) * S.x +
                      2 * c.charge_sign * mup * S.z) /
                     4;
    // Independent direct Jacobian from the note.
    double w = c.current_fraction * (field.p + 1) * pi * b.twist / b.r *
               c.dist.inv_abs_b_mean * c.dist.f(c.charge_sign * beta) *
               overlap * std::pow(1 - beta * b.mu, 2) * (1 - beta * beta) /
               std::abs(b.mu - beta);
    near(rate.weight, w, 1e-14, "Full branch contribution");
    expected += w;
  }
  require(
      weights[0].weight > 0 && weights[3].weight > 0,
      "Both charged populations contribute in synthetic equal-mass mixture");
  near(Evolution::total_rate(weights), expected, 1e-14, "Mixture additivity");
  for (double mup : {-1., -0.7, 0., 0.4, 1.}) {
    for (int sign : {-1, 1}) {
      near(Evolution::overlap(-1, 0, mup, sign), 0.5, 1e-15, "E overlap");
      near(Evolution::overlap(1, 0, mup, sign), mup * mup / 2, 1e-15,
           "O overlap");
    }
    near(Evolution::overlap(0, 1, mup, 1) - Evolution::overlap(0, 1, mup, -1),
         mup, 1e-15, "Charge reversal circular term");
  }
  auto original = mixed.comp;
  std::get<1>(mixed.comp).current_fraction = 0;
  auto electron = mixed.rates(b, S, false);
  near(Evolution::total_rate(electron), weights[0].weight + weights[1].weight,
       1e-15, "Zero current removes exactly its opacity");
  mixed.comp = original;

  // At mu=.5, x=.9 both roots are positive and inside the proton support.
  b.mu = 0.5;
  b.x_e = 1;
  auto two_roots = mixed.rates(b, {0, 0, 0}, false);
  require(two_roots[2].weight > 0 && two_roots[3].weight > 0,
          "Two retained roots of one population");
  for (auto const& rates : {weights, two_roots}) {
    std::array<int, 4> hits{};
    constexpr int count = 60000;
    for (int i = 0; i < count; ++i) {
      auto r = mixed.select_rate(rates);
      ++hits[2 * r.component + r.root];
    }
    double total = Evolution::total_rate(rates);
    for (int j = 0; j < 4; ++j) {
      double p = rates[j].weight / total;
      near(double(hits[j]) / count, p,
           6 * std::sqrt(p * (1 - p) / count) + 1e-12,
           "Joint component/root selection probability");
    }
  }

  // Actual proton mass: proton resonance must not reuse electron x.
  Evolution physical(1.4, 10, 9, field, source, components());
  physical.init();
  b.mu = 0;
  b.x_e = 1.15 * 1836.15267343;
  auto proton = physical.rates(b, {-1, 0, 0}, false);
  require(proton[3].weight > 0 && proton[0].weight == 0 &&
              proton[1].weight == 0,
          "Mass-specific resonance surfaces");
}

void
test_scattering_kinematics(BField const& field, InitFD11 const& source) {
  Evolution photon(1.4, 10, 88, field, source, components());
  auto initial = radial(50, 0.3, 2);
  constexpr int count = 30000;
  int e_modes = 0;
  double sum = 0, square = 0, azimuth = 0;
  for (int i = 0; i < count; ++i) {
    photon.init(initial);
    auto before = photon.basis({50, 0, 0});
    double beta = (i % 2 ? 1 : -1) * 0.6;
    photon.scattered_data = {{50, 0, 0}, beta};
    photon.scattering_pending = true;
    double target = photon.tau_target;
    photon.perform_scattering();
    auto after = photon.basis(photon.geodesic.stepper.y_old);
    near_vector(after.r_hat, before.r_hat, 2e-15,
                "Collision preserves position");
    near(after.r, before.r, 0, "Collision preserves radius");
    near(photon.omega_inf * (1 - beta * after.mu),
         initial.omega_inf * (1 - beta * before.mu), 3e-14,
         "Rest-frame photon energy conserved");
    near_vector(cross(photon.e1, photon.e2), photon.n, 2e-15,
                "Right-handed orbit");
    near(photon.n.length(), 1, 2e-15, "Unit plane normal");
    near(dot(photon.n, photon.e1), 0, 2e-15, "Orthogonal plane basis");
    require(photon.tau == 0 && photon.tau_target > 0 &&
                photon.tau_target != target,
            "Collision starts a fresh optical-depth clock");
    double mup = (after.mu - beta) / (1 - beta * after.mu);
    sum += mup;
    square += mup * mup;
    e_modes += photon.apol == AdiabaticPolarization::E;
    auto transverse = to_unit(cross(before.B, before.r_hat));
    azimuth +=
        dot(after.k_hat, transverse) / std::sqrt(1 - after.mu * after.mu);
  }
  near(sum / count, 0, 6 * std::sqrt(0.4 / count), "Rest-frame mean mu");
  near(square / count, 0.4, 6 * std::sqrt((9.0 / 35 - 0.16) / count),
       "Rest-frame angular second moment");
  near(double(e_modes) / count, 0.75, 6 * std::sqrt(0.75 * 0.25 / count),
       "Outgoing E-mode probability integrated over direction");
  near(azimuth / count, 0, 6 * std::sqrt(0.5 / count), "Uniform azimuth");
}

void
test_events(BField const& field, InitFD11 const& source) {
  Evolution photon(1.4, 10, 42, field, source, components());
  photon.init();
  double target = photon.tau_target;
  require(photon.evolve() == Evolution::EvolveResult::Scattered,
          "Seed 42 scatters");
  require(photon.scattered_data.beta < 0 &&
              photon.scattered_data.component == 0,
          "Seed 42 selects the electron component");
  near(photon.tau, target, 0, "Target is not resampled during propagation");
  auto event = photon.scattered_data;
  // Replay from exactly the same initial RNG seed; independently integrate on
  // many smaller intervals, stopping at the first run's event position.
  Evolution replay(1.4, 10, 42, field, source, components());
  replay.init();
  Quad independent;
  double depth = 0;
  while (replay.geodesic.l_old() < event.l) {
    replay.geodesic.step();
    double a = replay.geodesic.l_old(),
           end = std::min(event.l, replay.geodesic.l_new());
    auto cuts = replay.integration_cuts(a, end);
    auto rate = [&](double l) {
      return Evolution::total_rate(replay.rates(replay.get_basis(l), {}, true));
    };
    for (size_t j = 0; j + 1 < cuts.size(); ++j)
      for (int k = 0; k < 8; ++k) {
        double x = std::lerp(cuts[j], cuts[j + 1], k / 8.0);
        double y = std::lerp(cuts[j], cuts[j + 1], (k + 1) / 8.0);
        depth += independent.qags(rate, x, y, 1e-11, 1e-10);
      }
    if (end == event.l)
      break;
    replay.geodesic.stepper.update_old();
  }
  near(depth, target, 3e-8, "Independent scattering optical-depth residual");
  Evolution::Options fine;
  fine.n_knots = fine.n_scan = 16;
  fine.quad_atol = fine.quad_rtol = 1e-10;
  fine.root_atol = 1e-11;
  Evolution refined(1.4, 10, 42, field, source, components(), fine);
  refined.init();
  require(refined.evolve() == Evolution::EvolveResult::Scattered,
          "Refined scattering");
  near(refined.scattered_data.l, event.l, 2e-7, "Event position convergence");
  photon.perform_scattering();
  require(!photon.scattering_pending && !photon.finished,
          "Ready after collision");

  Evolution absorb(1.4, 10, 7, field, source, components(0, 0));
  auto inward = radial(12, 0.2);
  inward.alpha0 = pi;
  absorb.init(inward);
  require(absorb.evolve() == Evolution::EvolveResult::Absorbed,
          "Inward surface absorption");
  near(absorb.geodesic.stepper.y_new[0], 10, 2e-12,
       "Surface event root is inside step");
  inward.r0 = 10;
  absorb.init(inward);
  require(absorb.evolve() == Evolution::EvolveResult::Absorbed,
          "Inward birth at surface");
}

void
test_singular_depth(BField const& field, InitFD11 const& source) {
  // A flat radial ray has constant mu and x(r) proportional to r^(-2-p).
  // Choose mu inside the positive speed support so the E opacity diverges
  // integrably where D -> 0+.
  Evolution photon(0, 10, 9, field, source, components(0, 1));
  double z = 0.1;
  for (; z < 0.9; z += 0.01) {
    photon.init(radial(20, z));
    auto b = photon.basis({20, 0, 0});
    if (b.mu > 0.45 && b.mu < 0.55)
      break;
  }
  require(z < 0.9, "Find radial coalescence fixture");
  auto b0 = photon.basis({20, 0, 0});
  double ratio = std::get<1>(photon.comp).charge_mass_ratio;
  auto initial = radial(20, z, b0.x_e * ratio / 1.05);
  photon.init(initial);
  b0 = photon.basis({20, 0, 0});
  double edge =
      20 * std::pow(1.05 / std::sqrt(1 - b0.mu * b0.mu), 1 / (2 + field.p));
  double length = edge - 20 + 0.5;
  photon.geodesic.stepper.init(0, length, {20, 0, 0});
  photon.geodesic.stepper.do_step(length);
  // This test isolates adiabatic optical depth even if the energy is low.
  photon.is_adiabatic = true;
  photon.tau_target = 1e9;
  auto pol = [&](double l) { return photon.mode_at(l); };
  require(!photon.integrate_scattering(0, length, pol), "No artificial event");
  double coarse = photon.tau;
  require(coarse > 0, "Nonzero depth approaching coalescence");
  auto cuts = photon.integration_cuts(0, length);
  double closest = length;
  for (double cut : cuts)
    closest = std::min(closest, std::abs(cut - (edge - 20)));
  near(closest, 0, 2e-13, "Locate D=0 boundary");
  require(cuts.size() > 4, "Power-law support and internal knot crossings");
  photon.options.n_scan = 16;
  photon.knots[1].clear();
  for (int j = 0; j < 16; ++j)
    photon.knots[1].push_back(std::lerp(0.2, 0.8, j / 15.0));
  photon.tau = 0;
  photon.options.quad_atol = photon.options.quad_rtol = 1e-9;
  photon.integrate_scattering(0, length, pol);
  near(photon.tau, coarse, 2e-7, "Singular integral knot/scan convergence");

  // Independent radial formula and t=sqrt(edge-r) substitution. Midpoint
  // quadrature avoids evaluating the exact singularity or velocity endpoints.
  auto direct = [&](double r) {
    double x = 1.05 * std::pow(20 / r, 2 + field.p), mu = b0.mu;
    double D = x * x + mu * mu - 1;
    if (D <= 0)
      return 0.0;
    auto const& c = std::get<1>(photon.comp);
    double value = 0;
    for (int sign : {-1, 1}) {
      double beta = (mu + sign * x * std::sqrt(D)) / (x * x + mu * mu);
      value += c.dist.f(beta) * 0.5 * std::pow(1 - beta * mu, 2) *
               (1 - beta * beta) / std::abs(mu - beta);
    }
    return (field.p + 1) * pi * b0.twist * c.dist.inv_abs_b_mean / r * value;
  };
  constexpr int count = 300000;
  double tmax = std::sqrt(edge - 20), h = tmax / count, reference = 0;
  for (int i = 0; i < count; ++i) {
    double t = (i + 0.5) * h;
    reference += 2 * t * direct(edge - t * t) * h;
  }
  near(photon.tau, reference, 2e-4,
       "Independent transformed singular quadrature");
}

void
test_escape_and_polarization(BField const& field, InitFD11 const& source) {
  Evolution photon(1.4, 10, 5, field, source, components());
  auto initial = radial(10, 0.3);
  initial.alpha0 = 0.4;
  photon.init(initial);
  photon.tau_target =
      1e9; // Deterministically bypass scattering, retain all opacity.
  require(photon.evolve() == Evolution::EvolveResult::Escaped,
          "Bent photon escapes");
  auto escaped = photon.escaped_data;
  require(photon.resonance_escape(escaped.r_psi_alpha),
          "All species escape test");
  require(!photon.is_adiabatic && std::isfinite(photon.coupling_l),
          "Coupling detected");
  near(photon.stokes.length(), 1, 1e-11, "Pure Stokes norm after evolution");
  near(std::hypot(escaped.Q, escaped.U, escaped.V), 1, 1e-11,
       "Sky Stokes norm");
  near(escaped.k_inf.length(), 1, 2e-15, "Unit asymptotic direction");

  // Continue birefringence beyond the nominal escape point; compare the final
  // sky result. This evaluates the local stopping prescription, not just its
  // code.
  Stokes stopped = photon.stokes;
  double3 old_n = photon.n;
  double radius = escaped.r_psi_alpha[0];
  photon.geodesic.init(radius, escaped.r_psi_alpha[1], escaped.r_psi_alpha[2]);
  while (photon.geodesic.stepper.y_old[0] < 4 * radius) {
    photon.geodesic.step();
    StepperMagnus stepper(1e-9, 1e-9);
    double a = photon.geodesic.l_old(), b = photon.geodesic.l_new();
    stepper.init(a, b, 0.01 * (b - a), photon.stokes);
    auto omega = [&](double l) { return photon.get_basis(l).omega; };
    while (stepper.do_step(omega, b - a, 0.125))
      stepper.update_old();
    photon.stokes = stepper.y_old;
    photon.geodesic.stepper.update_old();
  }
  photon.record_escape(photon.geodesic.l_old(), photon.geodesic.stepper.y_old);
  auto continued = photon.escaped_data;
  near_vector(continued.k_inf, escaped.k_inf, 2e-9,
              "Asymptotic bending consistency");
  double drift = (photon.stokes + (-1.0) * stopped).length();
  require(drift < 2e-3, "Representative escape continuation drift");
  near_vector({continued.Q, continued.U, continued.V},
              {escaped.Q, escaped.U, escaped.V}, 2e-3,
              "Observer polarization continuation");
  std::println("Escape r={} coupling_l={} continuation |delta S|={}", radius,
               photon.coupling_l, drift);

  // Global-axis fallback and screen-rotation covariance; no magnetic field
  // input.
  Stokes S{0.3, -0.4, std::sqrt(0.75)};
  auto polar = Evolution::observer_stokes(S, {0, 0, 1}, {0, 1, 0});
  near_vector(polar, S, 1e-15, "Global +x fallback at +z view");
  auto k = escaped.k_inf;
  double angle = 0.43, c = std::cos(angle), s = std::sin(angle);
  double3 old_s1 = cross(old_n, k), new_n = (-s) * old_s1 + c * old_n;
  Stokes rotated{std::cos(2 * angle) * S.x + std::sin(2 * angle) * S.y,
                 -std::sin(2 * angle) * S.x + std::cos(2 * angle) * S.y, S.z};
  near_vector(Evolution::observer_stokes(S, k, old_n),
              Evolution::observer_stokes(rotated, k, new_n), 1e-14,
              "Observer output invariant under transported-screen rotation");

  // Constant-field direction along a radial ray preserves its exact normal
  // mode.
  Evolution radial_ray(1.4, 10, 5, field, source, components());
  radial_ray.init(radial(10, 0.3));
  radial_ray.tau_target = 1e9;
  require(radial_ray.evolve() == Evolution::EvolveResult::Escaped,
          "Radial escape");
  auto radial_basis = radial_ray.basis(radial_ray.escaped_data.r_psi_alpha);
  near_vector(radial_ray.stokes, radial_ray.eigenmode(radial_basis), 1e-10,
              "Radial normal mode retained");
  // Proton already outside resonance is insufficient while electrons can
  // resonate.
  auto y = radial_ray.escaped_data.r_psi_alpha;
  y[0] = 20;
  require(!radial_ray.resonance_escape(y),
          "All modeled species must pass escape");
}

void
test_stokes_collision_and_coupling(BField const& field,
                                   InitFD11 const& source) {
  Evolution photon(1.4, 10, 17, field, source, components());
  auto initial = radial(500, -0.3);
  photon.init(initial);
  initial.omega_inf = photon.basis({500, 0, 0}).x_e / 1.05;
  photon.init(initial);
  require(!photon.is_adiabatic, "Initial resonance in Stokes region");
  photon.stokes = {0, 1, 0};
  photon.tau_target = 0.005;
  require(photon.evolve() == Evolution::EvolveResult::Scattered,
          "Collision while explicitly transporting Stokes");
  auto event = photon.scattered_data;
  require(!event.incoming_adiabatic && std::abs(event.incoming_stokes.z) > 1e-8,
          "Scattering reads evolved, not initial, polarization");
  near(event.incoming_stokes.length(), 1, 2e-13, "Dense event Stokes norm");
  Evolution::Options options;
  options.n_knots = options.n_scan = 16;
  options.max_phase = 0.125;
  options.stokes_atol = options.stokes_rtol = 1e-9;
  options.quad_atol = options.quad_rtol = 1e-10;
  Evolution fine(1.4, 10, 17, field, source, components(), options);
  fine.init(initial);
  fine.stokes = {0, 1, 0};
  fine.tau_target = 0.005;
  require(fine.evolve() == Evolution::EvolveResult::Scattered,
          "Refined Stokes collision");
  near(fine.scattered_data.l, event.l, 1e-6,
       "Stokes collision position converges");
  near_vector(fine.scattered_data.incoming_stokes, event.incoming_stokes, 1e-8,
              "Stokes collision polarization converges");
  photon.perform_scattering();
  require(!photon.is_adiabatic,
          "Scattered state already beyond coupling threshold");
  near(photon.stokes.z, 0, 0, "Collision resets circular polarization");
  near_vector(photon.stokes,
              photon.eigenmode(photon.basis(photon.geodesic.stepper.y_old)),
              1e-15, "New mode in new orbit screen");

  Evolution crossing(1.4, 10, 17, field, source, components(0, 0));
  crossing.init(radial(10, 0.3));
  require(crossing.evolve() == Evolution::EvolveResult::Escaped,
          "Locate coupling before escape");
  // For a radial ray the coupling function can be evaluated independently at
  // radius; bisection's recorded path must agree with its geodesic mapping.
  double r0 = 10, r1 = crossing.escaped_data.r_psi_alpha[0];
  auto f = [&](double r) {
    auto b = crossing.basis({r, 0, 0});
    return crossing.options.eta_couple * b.kappa * b.r - 1;
  };
  double r_couple = bisect_root(f, r0, r1, f(r0), f(r1), 1e-10);
  Quad q;
  auto radial_length = [&](double r) {
    return 1 / std::sqrt(1 - crossing.geodesic.rs / r);
  };
  near(crossing.coupling_l, q.qags(radial_length, 10, r_couple, 1e-9, 1e-10),
       2e-7, "Coupling located inside orbit step, not at its endpoint");

  Evolution aligned(1.4, 10, 5, field, source, components());
  aligned.init(radial(1000, 1));
  auto b = aligned.basis({1000, 0, 0});
  near(b.Bperp2, 0, 0, "Exactly aligned fixture");
  near_vector(b.omega, {}, 0, "Regular generator at B_perp=0");
  near_vector(aligned.stokes, {-1, 0, 0}, 0,
              "Specified birth screen for degenerate E mode");
}

void
test_repeated_events(BField const& field, InitFD11 const& source) {
  int escaped = 0, absorbed = 0, scattered = 0, proton = 0;
  for (int seed = 1; seed <= 32; ++seed) {
    Evolution photon(1.4, 10, seed, field, source, components());
    photon.init();
    for (int event = 0;; ++event) {
      require(event < 1000, "Repeated event termination");
      auto result = photon.evolve();
      if (result == Evolution::EvolveResult::Escaped) {
        ++escaped;
        break;
      }
      if (result == Evolution::EvolveResult::Absorbed) {
        ++absorbed;
        break;
      }
      ++scattered;
      proton += photon.scattered_data.component == 1;
      photon.perform_scattering();
    }
  }
  require(escaped > 0 && absorbed > 0 && scattered > 0 && proton > 0,
          "Regression sample exercises both species and terminal events");
  std::println(
      "32 photons: escaped={} absorbed={} collisions={} proton_collisions={}",
      escaped, absorbed, scattered, proton);
}

void
test_magnus_phase() {
  StepperMagnus stepper;
  auto omega = [](double x) -> double3 { return {10 + 200 * x, 0, 0}; };
  stepper.init(0, 1, 1, {1, 0, 0}); // exact eigenmode: embedded error is zero
  int count = 0;
  while (stepper.do_step(omega, 0.1, 0.25)) {
    require((10 + 200 * stepper.x_new) * stepper.h_old <= 0.25 * (1 + 1e-14),
            "Phase cap checks growth inside accepted step");
    stepper.update_old();
    ++count;
  }
  require(count > 400, "Eigenmode still resolves phase");
}
} // namespace

int
main() try {
  BField field("table/bfield_t10.txt", 1e14, 10);
  InitFD11 source(10, 0.4);
  test_rates(field, source);
  std::println("Rates and joint selection passed");
  test_scattering_kinematics(field, source);
  std::println("Scattering kinematics and sampling passed");
  test_events(field, source);
  std::println("Scattering events and absorption passed");
  test_singular_depth(field, source);
  std::println("Singular optical depth passed");
  test_escape_and_polarization(field, source);
  test_stokes_collision_and_coupling(field, source);
  test_repeated_events(field, source);
  test_magnus_phase();
  std::println("All photon evolution checks passed");
} catch (std::exception const& e) {
  std::println(stderr, "FAIL: {}", e.what());
  return 1;
}
