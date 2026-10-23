#include "bfield.hpp"
#include "component.hpp"
#include "find_root.hpp"
#include "geodesic.hpp"
#include "initialize.hpp"
#include "magnus.hpp"
#include "ran.hpp"
#include <array>
#include <limits>
#include <optional>
#include <print>
#include <tuple>
#include <vector>

// Local lengths in km, photon energies at infinity in keV, fields in gauss.
// V = +2 Im(A1 conj(A2)); s1 = n cross k, s2 = n.
template <PhotonInitializer Init, typename Comp>
struct PhotonEvolution {
  struct Options {
    double eta_couple = 1e-3, escape_tolerance = 1e-3;
    double quad_atol = 1e-8, quad_rtol = 1e-8, root_atol = 1e-10;
    double stokes_atol = 1e-6, stokes_rtol = 1e-6, max_phase = 0.25;
    int n_knots = 8, n_scan = 8;
  } options;
  BField const& bfield;
  Init const& initializer;
  Comp comp; // Own the tuple: callers may pass std::make_tuple(...).
  Ran ran;
  Geodesic geodesic;
  Quad quad;
  using State = Geodesic::YVector;
  static constexpr size_t component_count = std::tuple_size_v<Comp>;
  std::array<std::vector<double>, component_count> knots;

  double3 n{}, e1{}, e2{};
  double omega_inf = 0, tau = 0, tau_target = 0, stokes_h = 0;
  double coupling_l = std::numeric_limits<double>::infinity();
  bool is_adiabatic = true;
  AdiabaticPolarization apol = AdiabaticPolarization::E;
  Stokes stokes{};
  bool initialized = false, finished = false, scattering_pending = false;

  PhotonEvolution(double M_star, double R_star, int seed, BField const& field,
                  Init const& init, Comp const& components)
      : PhotonEvolution(M_star, R_star, seed, field, init, components,
                        Options{}) {}

  PhotonEvolution(double M_star, double R_star, int seed, BField const& field,
                  Init const& init, Comp const& components, Options opts)
      : options(opts), bfield(field), initializer(init), comp(components),
        ran(seed), geodesic(M_star, R_star) {
    for (double v :
         {options.eta_couple, options.escape_tolerance, options.quad_atol,
          options.quad_rtol, options.root_atol, options.stokes_atol,
          options.stokes_rtol, options.max_phase})
      if (!(std::isfinite(v) && v > 0))
        throw std::invalid_argument(
            "Evolution tolerances must be positive and finite");
    if (options.n_knots < 2 || options.n_scan < 2)
      throw std::invalid_argument(
          "At least two knots and scan intervals required");
    for_components([&](size_t i, auto const& c) {
      if ((c.charge_sign != -1 && c.charge_sign != 1) ||
          !std::isfinite(c.current_fraction) || c.current_fraction < 0 ||
          !(std::isfinite(c.charge_mass_ratio) && c.charge_mass_ratio > 0))
        throw std::invalid_argument("Invalid current component");
      for (int j = 0; j < options.n_knots; ++j)
        knots[i].push_back(c.charge_sign *
                           std::lerp(c.dist.b_min, c.dist.b_max,
                                     double(j) / (options.n_knots - 1)));
    });
  }

  template <class F>
  void
  for_components(F const& f) const {
    size_t i = 0;
    std::apply([&](auto const&... c) { (f(i++, c), ...); }, comp);
  }

  struct Basis {
    double3 r_hat, psi_hat, k_hat, s1, s2, B, omega;
    double Bperp2, r, c, s, kappa, Bmag, mu, x_e, twist;
    Basis(State y, double3 n, double3 e1, double3 e2, BField const& field,
          double energy, double rs) {
      auto [radius, psi, alpha] = y;
      r = radius;
      r_hat = std::cos(psi) * e1 + std::sin(psi) * e2;
      psi_hat = cross(n, r_hat);
      k_hat = std::cos(alpha) * r_hat + std::sin(alpha) * psi_hat;
      s1 = cross(n, k_hat);
      s2 = n;
      double muz = std::clamp(r_hat.z, -1.0, 1.0);
      double rho = std::hypot(r_hat.x, r_hat.y);
      double3 theta_hat =
          rho > 0 ? double3{r_hat.x * muz / rho, r_hat.y * muz / rho, -rho}
                  : double3{std::copysign(1.0, muz), 0, 0};
      double3 phi_hat = rho > 0 ? double3{-r_hat.y / rho, r_hat.x / rho, 0}
                                : double3{0, 1, 0};
      auto [B_vec, t] = field.B_and_twist(r, muz);
      twist = t;
      B = B_vec.x * r_hat + B_vec.y * theta_hat + B_vec.z * phi_hat;
      Bmag = B.length();
      mu = Bmag > 0 ? std::clamp(dot(k_hat, B) / Bmag, -1.0, 1.0) : 0;
      double L = std::sqrt(1 - rs / r);
      x_e = B_to_omega * Bmag * L / energy;
      double B1 = dot(B, s1), B2 = dot(B, s2);
      Bperp2 = B1 * B1 + B2 * B2;
      double diff = B1 * B1 - B2 * B2, twice = 2 * B1 * B2;
      // At birth exactly along B, choose O along s1. For an incoming label
      // at a degeneracy, mode_at() supplies the continuous incoming limit.
      c = Bperp2 > 0 ? diff / Bperp2 : 1;
      s = Bperp2 > 0 ? twice / Bperp2 : 0;
      double K = 2.013919e-19 * energy / L;
      kappa = K * Bperp2;
      omega = {K * diff, K * twice, 0};
    }
  };

  Basis
  basis(State y) const {
    return Basis(y, n, e1, e2, bfield, omega_inf, geodesic.rs);
  }
  Basis
  get_basis(double l) const {
    return basis(geodesic.interp(l));
  }

  Stokes
  eigenmode(Basis const& b) const {
    double sign = apol == AdiabaticPolarization::O ? 1 : -1;
    return {sign * b.c, sign * b.s, 0};
  }
  Stokes
  mode_at(double l) const {
    auto b = get_basis(l);
    if (b.Bperp2 == 0 && l > geodesic.l_old()) {
      double d = std::min(l - geodesic.l_old(), 1e-7 * b.r);
      b = get_basis(l - d);
    }
    if (b.Bperp2 == 0)
      return stokes; // Incoming limiting mode, or the specified birth screen.
    return eigenmode(b);
  }

  void
  init() {
    init(initializer.sample_initial_photon(ran));
  }
  void
  init(InitializedPhoton const& photon) {
    geodesic.init(photon.r0, photon.psi0, photon.alpha0);
    n = photon.n;
    e1 = photon.e1;
    e2 = photon.e2;
    omega_inf = photon.omega_inf;
    apol = photon.pol;
    stokes_h = 0;
    tau = 0;
    tau_target = -std::log(ran.U());
    auto b = basis(geodesic.stepper.y_old);
    stokes = eigenmode(b);
    is_adiabatic = options.eta_couple * b.kappa * b.r > 1;
    coupling_l = is_adiabatic ? std::numeric_limits<double>::infinity() : 0;
    initialized = true;
    finished = scattering_pending = false;
  }

  static
  double
  discriminant(double x, double mu) {
    return std::fma(x, x, (mu - 1) * (mu + 1));
  }

  // Sorted roots of (x^2+mu^2) beta^2 - 2 mu beta + (1-x^2)=0.
  // Compute one root without cancellation, the other by their product.
  static
  std::array<double, 2>
  resonance_roots(double x, double mu, double D) {
    double a = x * x + mu * mu;
    double q = mu + std::copysign(x * std::sqrt(D), mu);
    std::array<double, 2> roots{q / a, ((1 - x) * (1 + x)) / q};
    if (roots[1] < roots[0])
      std::swap(roots[0], roots[1]);
    return roots;
  }

  static
  double
  overlap(double qb, double v, double mup, int charge_sign) {
    double w = 0.25 * ((1 + mup * mup) + (mup * mup - 1) * qb) +
               0.5 * charge_sign * mup * v;
    if (!std::isfinite(w) || w < -1e-10 || w > 1 + 1e-10)
      throw std::runtime_error("Nonphysical polarization overlap");
    return std::clamp(w, 0.0, 1.0);
  }

  struct Rate {
    size_t component = 0;
    int root = 0; // 0 = lower beta root, 1 = upper beta root
    double beta = 0, weight = 0;
  };
  using Rates = std::array<Rate, 2 * component_count>;

  Rates
  rates(Basis const& b, Stokes S, bool adiabatic) const {
    Rates result{};
    double qb = adiabatic ? (apol == AdiabaticPolarization::O ? 1 : -1)
                          : b.c * S.x + b.s * S.y;
    double v = adiabatic ? 0 : S.z;
    for_components([&](size_t i, auto const& c) {
      if (c.current_fraction == 0 || b.twist == 0)
        return;
      double x = c.charge_mass_ratio * b.x_e, D = discriminant(x, b.mu);
      if (!(D > 0) || !(x > 0))
        return; // Exact coalescence is an integration boundary, not a sample.
      auto roots = resonance_roots(x, b.mu, D);
      for (int j = 0; j < 2; ++j) {
        double beta = roots[j];
        double f = c.dist.f(c.charge_sign * beta);
        if (!(beta > -1 && beta < 1) || f == 0)
          continue;
        double t = (1 - beta) * (1 + beta);
        double mup = std::clamp((b.mu - beta) / (1 - beta * b.mu), -1.0, 1.0);
        // Equivalent to (1-beta*mu)^2(1-beta^2)/abs(mu-beta),
        // but avoids cancellation at coalescing roots.
        double jacobian = x * x * t * std::sqrt(t) / std::sqrt(D);
        double w = c.current_fraction * (bfield.p + 1) * pi * b.twist *
                   c.dist.inv_abs_b_mean / b.r * f * jacobian *
                   overlap(qb, v, mup, c.charge_sign);
        if (!std::isfinite(w) || w < 0)
          throw std::runtime_error("Invalid differential optical depth");
        result[2 * i + j] = {i, j, beta, w};
      }
    });
    return result;
  }
  static
  double
  total_rate(Rates const& rates) {
    double sum = 0;
    for (auto const& r : rates)
      sum += r.weight;
    return sum;
  }
  Rate
  select_rate(Rates const& rates) {
    double sum = total_rate(rates);
    if (!(sum > 0 && std::isfinite(sum)))
      throw std::runtime_error("No positive rate at scattering event");
    double target = ran.U() * sum, cumulative = 0;
    Rate last{};
    for (auto const& r : rates) {
      if (r.weight == 0)
        continue;
      last = r;
      cumulative += r.weight;
      if (target < cumulative)
        return r;
    }
    return last; // Only roundoff in the final cumulative sum.
  }

  // The union of EVERY component's resonance boundary and velocity knots.
  std::vector<double>
  integration_cuts(double left, double right) const {
    // Also split the join of the two Magnus dense-output polynomials.
    std::vector<double> cuts{left, std::midpoint(left, right), right};
    auto scan = [&](auto const& f) {
      double a = left, fa = f(a);
      for (int j = 1; j <= options.n_scan; ++j) {
        double b = std::lerp(left, right, double(j) / options.n_scan),
               fb = f(b);
        if (fa == 0)
          cuts.push_back(a);
        if (fb == 0)
          cuts.push_back(b);
        if ((fa < 0 && fb > 0) || (fa > 0 && fb < 0))
          cuts.push_back(bisect_root(f, a, b, fa, fb));
        a = b;
        fa = fb;
      }
    };
    for_components([&](size_t i, auto const& c) {
      if (c.current_fraction == 0)
        return;
      scan([&](double l) {
        auto b = get_basis(l);
        return discriminant(c.charge_mass_ratio * b.x_e, b.mu);
      });
      for (double beta : knots[i]) {
        double inv_gamma = std::sqrt((1 - beta) * (1 + beta));
        scan([&](double l) {
          auto b = get_basis(l);
          return c.charge_mass_ratio * b.x_e * inv_gamma + beta * b.mu - 1;
        });
      }
    });
    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    return cuts;
  }

  enum class EvolveResult { Absorbed, Escaped, Scattered };
  struct EscapedData {
    double omega_inf = 0, muk = 0, Q = 0, U = 0, V = 0;
    double l = 0;
    State r_psi_alpha{};
    double3 k_inf{};
  } escaped_data;
  struct ScatteredData {
    State r_psi_alpha{};
    double beta = 0;
    size_t component = 0;
    int root = 0;
    double l = 0;
    Stokes incoming_stokes{};
    bool incoming_adiabatic = true;
    AdiabaticPolarization incoming_mode = AdiabaticPolarization::E;
  } scattered_data;

  template <class PolarizationAt>
  bool
  integrate_scattering(double left, double right, PolarizationAt const& pol) {
    auto cuts = integration_cuts(left, right);
    auto integrand = [&](double l) {
      return total_rate(rates(get_basis(l), pol(l), is_adiabatic));
    };
    for (size_t j = 0; j + 1 < cuts.size(); ++j) {
      double a = cuts[j], b = cuts[j + 1];
      if (!(a < b))
        continue;
      auto integral = [&](double end) {
        return quad.qags(integrand, a, end, options.quad_atol,
                         options.quad_rtol);
      };
      double dtau = integral(b);
      if (dtau < 0)
        throw std::runtime_error("Negative integrated optical depth");
      if (dtau > 0 && tau + dtau >= tau_target) {
        double remaining = tau_target - tau;
        auto residual = [&](double l) { return integral(l) - remaining; };
        double l = bisect_root(residual, a, b, -remaining, dtau - remaining,
                               options.root_atol);
        Stokes incoming = pol(l);
        Rate event = select_rate(rates(get_basis(l), incoming, is_adiabatic));
        scattered_data = {geodesic.interp(l), event.beta, event.component,
                          event.root,         l,          incoming,
                          is_adiabatic,       apol};
        stokes = incoming;
        tau = tau_target;
        scattering_pending = true;
        return true;
      }
      tau += dtau;
    }
    return false;
  }

  // Find the FIRST coupling crossing, including one inside a geodesic step.
  std::optional<double>
  coupling_point(double left, double right) const {
    auto f = [&](double l) {
      auto b = get_basis(l);
      return options.eta_couple * b.kappa * b.r - 1;
    };
    double a = left, fa = f(a);
    if (fa <= 0)
      return a;
    for (int j = 1; j <= options.n_scan; ++j) {
      double b = std::lerp(left, right, double(j) / options.n_scan), fb = f(b);
      if (fb <= 0)
        return bisect_root(f, a, b, fa, fb, options.root_atol);
      a = b;
      fa = fb;
    }
    return {};
  }

  // Symmetric trial propagation includes bending and rotating local bases.
  // No dense-output extrapolation and no mutation of the stored geodesic.
  State
  trial_state(State y, double h) const {
    auto displaced = [](State a, State const& k, double scale) {
      for (size_t i = 0; i < a.size(); ++i)
        a[i] += scale * k[i];
      return a;
    };
    State k1, k2, k3, k4;
    geodesic(0, y, k1);
    geodesic(0, displaced(y, k1, h / 2), k2);
    geodesic(0, displaced(y, k2, h / 2), k3);
    geodesic(0, displaced(y, k3, h), k4);
    for (size_t i = 0; i < y.size(); ++i)
      y[i] += h / 6 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]);
    return y;
  }

  bool
  resonance_escape(State y) const {
    if (std::cos(y[2]) <= 0)
      return false; // Infinity reconstruction is only for the outgoing branch.
    auto b = basis(y);
    double h = std::min(1e-4 * b.r, 0.1 * (b.r - geodesic.rs));
    auto minus = basis(trial_state(y, -h)), plus = basis(trial_state(y, h));
    double dx = (plus.x_e - minus.x_e) / (2 * h);
    double dmu = (plus.mu - minus.mu) / (2 * h);
    bool escaped = true;
    for_components([&](size_t, auto const& c) {
      if (c.current_fraction == 0 || bfield.A == 0)
        return;
      double x = c.charge_mass_ratio * b.x_e;
      double tx = c.charge_mass_ratio * dx;
      escaped = escaped && discriminant(x, b.mu) < 0 && tx < 0 &&
                std::abs(b.mu - x * dmu / tx) < 1;
    });
    return escaped;
  }

  static
  Stokes
  observer_stokes(Stokes S, double3 k, double3 normal) {
    double3 axis{0, 0, 1};
    double3 u = axis + (-dot(axis, k)) * k;
    if (u.length() < 1e-12) {
      axis = {1, 0, 0};
      u = axis + (-dot(axis, k)) * k;
    }
    u = to_unit(u);
    double c = dot(u, cross(normal, k)), s = dot(u, normal);
    double norm = std::hypot(c, s);
    c /= norm;
    s /= norm;
    double c2 = c * c - s * s, s2 = 2 * c * s;
    return {c2 * S.x + s2 * S.y, -s2 * S.x + c2 * S.y, S.z};
  }

  void
  record_escape(double l, State y) {
    auto [r, psi, alpha] = y;
    double impact = r * std::sin(alpha) / std::sqrt(1 - geodesic.rs / r);
    // u = r/rho maps [r,infinity) to [1,0], with a regular integrand at 0.
    auto integrand = [&](double u) {
      double t = impact / r;
      return t / std::sqrt(1 - t * t * u * u * (1 - geodesic.rs * u / r));
    };
    double psi_inf = psi + quad.qags(integrand, 0, 1, 1e-11, 1e-11);
    double3 k = to_unit(std::cos(psi_inf) * e1 + std::sin(psi_inf) * e2);
    Stokes sky = observer_stokes(stokes, k, n);
    escaped_data = {
        omega_inf, std::clamp(k.z, -1.0, 1.0), sky.x, sky.y, sky.z, l, y, k};
    finished = true;
  }

  EvolveResult
  evolve() {
    if (!initialized || finished || scattering_pending)
      throw std::logic_error(
          "Initialize a photon or perform its pending scattering first");
    while (true) {
      // An inward photon initialized exactly on the surface is absorbed there.
      if (geodesic.stepper.y_old[0] <= geodesic.R_star &&
          std::cos(geodesic.stepper.y_old[2]) < 0) {
        finished = true;
        return EvolveResult::Absorbed;
      }
      bool absorbed = geodesic.step();
      double left = geodesic.l_old(), right = geodesic.l_new();
      if (is_adiabatic) {
        auto crossing = coupling_point(left, right);
        double end = crossing.value_or(right);
        if (integrate_scattering(left, end,
                                 [&](double l) { return mode_at(l); }))
          return EvolveResult::Scattered;
        if (crossing) {
          stokes = mode_at(*crossing);
          is_adiabatic = false;
          coupling_l = left = *crossing;
        }
      }
      if (!is_adiabatic && left < right) {
        StepperMagnus stepper(options.stokes_atol, options.stokes_rtol);
        stepper.init(left, right,
                     stokes_h > 0 ? stokes_h : 0.1 * (right - left), stokes);
        auto omega = [&](double l) { return get_basis(l).omega; };
        while (stepper.do_step(omega, right - left, options.max_phase)) {
          // Integrate each component's knots BEFORE committing this dense step.
          if (integrate_scattering(stepper.x_old, stepper.x_new, [&](double l) {
                return stepper.dense_out(l);
              }))
            return EvolveResult::Scattered;
          stokes = stepper.y_new;
          stokes_h = stepper.h_new;
          auto y = geodesic.interp(stepper.x_new);
          double change = (stepper.y_new + (-1.0) * stepper.y_old).length() *
                          y[0] / (stepper.x_new - stepper.x_old);
          if (!(absorbed && stepper.x_new == right) &&
              change < options.escape_tolerance && resonance_escape(y)) {
            record_escape(stepper.x_new, y);
            return EvolveResult::Escaped;
          }
          stepper.update_old();
        }
      }
      if (absorbed) {
        finished = true;
        return EvolveResult::Absorbed;
      }
      if (is_adiabatic)
        stokes = mode_at(right);
      geodesic.stepper.update_old();
    }
  }

  static
  double
  sample_mup(Ran& ran) {
    double const u = ran.U();
    if (u < 0.75)
      return (8.0 / 3.0) * u - 1.0;
    double const y = 8.0 * u - 7.0;
    double m = std::max(std::abs(y), ran.U());
    m = std::max(m, ran.U());
    return std::copysign(m, y);
  }

  void
  perform_scattering() {
    if (!scattering_pending)
      throw std::logic_error("No pending scattering event");
    auto b = basis(scattered_data.r_psi_alpha);
    double beta = scattered_data.beta, mup = sample_mup(ran);
    double mu_out = std::clamp((mup + beta) / (1 + beta * mup), -1.0, 1.0);
    double3 bhat = b.B / b.Bmag;
    double3 k = mu_out * bhat +
                std::sqrt((1 - mu_out) * (1 + mu_out)) * ran.unit_perp_to(bhat);
    k = to_unit(k);
    double3 normal = cross(b.r_hat, k);
    double cross_len = normal.length();
    normal = cross_len > 1e-14 ? normal / cross_len : ran.unit_perp_to(b.r_hat);
    double energy = omega_inf * (1 - beta * b.mu) / (1 - beta * mu_out);
    auto mode = ran.U() < 1 / (1 + mup * mup) ? AdiabaticPolarization::E
                                              : AdiabaticPolarization::O;
    init({normal, b.r_hat, cross(normal, b.r_hat), energy, b.r, 0,
          std::atan2(cross_len, dot(k, b.r_hat)), mode});
  }
};

#ifndef MAGMC_NO_MAIN
int
main() try {
  double M_star = 1.4, R_star = 10, B_pole = 1e14, T_inf = 0.4;
  BField bfield("table/bfield_t10.txt", B_pole, R_star);
  InitFD11 initializer(R_star, T_inf);
  auto powerlaw_electron = CurrentComponent{.dist = PowerLaw2(0.2, 0.8),
                                            .charge_sign = -1,
                                            .current_fraction = 0.5,
                                            .charge_mass_ratio = 1.0};
  auto powerlaw_proton =
      CurrentComponent{.dist = PowerLaw2(0.2, 0.8),
                       .charge_sign = 1,
                       .current_fraction = 0.5,
                       .charge_mass_ratio = 1.0 / 1836.15267343};
  PhotonEvolution photon(M_star, R_star, 42, bfield, initializer,
                         std::make_tuple(powerlaw_electron, powerlaw_proton));
  photon.init();
  auto result = photon.evolve();
  if (result == decltype(photon)::EvolveResult::Scattered) {
    std::println("Scattered: component={} beta={} l={} tau={}",
                 photon.scattered_data.component, photon.scattered_data.beta,
                 photon.scattered_data.l, photon.tau);
    photon.perform_scattering();
  } else if (result == decltype(photon)::EvolveResult::Escaped) {
    auto const& e = photon.escaped_data;
    std::println("Escaped: omega_inf={} muk={} Q={} U={} V={}", e.omega_inf,
                 e.muk, e.Q, e.U, e.V);
  } else {
    std::println("Absorbed");
  }
} catch (std::exception const& e) {
  std::println(stderr, "magmc2: {}", e.what());
  return 1;
}
#endif
