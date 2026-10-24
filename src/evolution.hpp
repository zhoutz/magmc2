#pragma once

#include "bfield.hpp"
#include "find_root.hpp"
#include "geodesic.hpp"
#include "initialize.hpp"
#include "magnus.hpp"
#include "polarization.hpp"
#include "quad.hpp"
#include "ran.hpp"
#include "transport_cache.hpp"

template <PhotonInitializer Init, typename Comp>
struct PhotonEvolution {
  struct Options {
    double eta_couple = 1e-3, escape_tolerance = 1e-3;
    double quad_atol = 1e-8, quad_rtol = 1e-8;
    int n_knots = 8, n_scan = 8;
    bool resonance_bounds = true;
    double cache_tolerance = 0, event_rtol = 0; // Zero retains exact/reference paths.
  } options;

  BField const& bfield;
  Init const& initializer;
  Comp comp;
  Ran ran;
  Geodesic geodesic;
  Quad quad;

  double3 n, e1, e2;
  double omega_inf, tau;
  Polarization pol;

  std::vector<double> knots;

  PhotonEvolution(double M_star, double R_star, int seed, BField const& bfield,
                  Init const& init, Comp const& comp)
      : PhotonEvolution(M_star, R_star, seed, bfield, init, comp, Options{}) {}

  PhotonEvolution(double M_star, double R_star, int seed, BField const& bfield,
                  Init const& init, Comp const& comp, Options options)
      : options(options), bfield(bfield), initializer(init), comp(comp),
        ran(seed), geodesic(M_star, R_star) {
    for (int i = 0; i < options.n_knots; ++i)
      knots.push_back(comp.charge_sign *
                      std::lerp(comp.dist.b_min, comp.dist.b_max,
                                double(i) / (options.n_knots - 1)));
  }

  struct Transport {
    double r, c, s, kappa, mu, x_e, twist;
    double3 omega;
  };
  TransportCache cache;
  struct Basis {
    double3 r_hat, psi_hat, k_hat, s1, s2, B, omega;
    double Bperp2, r, c, s, kappa, Bmag, mu, x_e, twist;
    Basis(double radius, double psi, double alpha, double3 n, double3 e1,
          double3 e2, BField const& field, double energy, double rs) {
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
      c = Bperp2 > 0 ? diff / Bperp2 : 1;
      s = Bperp2 > 0 ? twice / Bperp2 : 0;
      double K = 2.013919e-19 * energy / L;
      kappa = K * Bperp2;
      omega = {K * diff, K * twice, 0};
    }
  };

  Basis
  basis(Geodesic::YVector y) const {
    auto [r, psi, alpha] = y;
    return Basis(r, psi, alpha, n, e1, e2, bfield, omega_inf, geodesic.rs);
  }

  Basis
  get_basis(double l) const {
    auto y = geodesic.interp(l);
    return basis(y);
  }

  TransportCache::Values sample_transport(double l) const {
    auto b = get_basis(l);
    return {b.r, b.mu, b.x_e, b.twist, dot(b.B, b.s1), dot(b.B, b.s2),
            2.013919e-19 * omega_inf / std::sqrt(1 - geodesic.rs / b.r)};
  }

  Transport transport_exact(double l) const {
    auto b = get_basis(l);
    return {b.r, b.c, b.s, b.kappa, b.mu, b.x_e, b.twist, b.omega};
  }

  Transport transport(double l) const {
    if (options.cache_tolerance == 0) return transport_exact(l);
    auto v = cache.at(l, [&](double x) { return sample_transport(x); });
    auto [r, mu, x, twist, b1, b2, K] = v;
    if (!(r > geodesic.rs && x > 0 && twist >= 0)) return transport_exact(l);
    // Reconstruct the rotation axis from transverse field components rather
    // than interpolating c/s independently; c*c+s*s remains one.
    double perp = b1 * b1 + b2 * b2, diff = b1 * b1 - b2 * b2, twice = 2 * b1 * b2;
    return {r, perp > 0 ? diff / perp : 1, perp > 0 ? twice / perp : 0, K * perp,
            std::clamp(mu, -1.0, 1.0), x, twist, {K * diff, K * twice, 0}};
  }

  void
  init() {
    init(initializer.sample_initial_photon(ran));
  }
  void
  init(InitializedPhoton const& photon) {
    cache.leaves.clear();
    geodesic.init(photon.r0, photon.psi0, photon.alpha0);
    n = photon.n;
    e1 = photon.e1;
    e2 = photon.e2;
    omega_inf = photon.omega_inf;
    tau = std::log(ran.U());
    pol.set_adiabatic(photon.pol);
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
    int iroot = 0; // 0 = lower beta root, 1 = upper beta root
    double beta = 0, weight = 0;
  };
  using Rates = std::array<Rate, 2>;

  Rates
  rates(Basis const& b, Polarization const& pol) const {
    return rates_at(b, [&] { return pol; });
  }

  // Dense polarization output requires rotations. Most quadrature samples
  // have no resonant particles, so evaluate it only for an active channel.
  template <class Geometry, class PolarizationAt>
  Rates
  rates_at(Geometry const& b, PolarizationAt const& pol_at) const {
    Rates result{};
    if (comp.current_fraction == 0 || b.twist == 0)
      return result;
    double x = comp.charge_mass_ratio * b.x_e, D = discriminant(x, b.mu);
    // Exact coalescence is an integration boundary, not a sample.
    if (!(D > 0) || !(x > 0))
      return result;
    auto roots = resonance_roots(x, b.mu, D);
    Stokes sb{};
    bool have_stokes = false;
    for (int j = 0; j < 2; ++j) {
      double beta = roots[j];
      double f = comp.dist.f(comp.charge_sign * beta);
      if (!(beta > -1 && beta < 1) || f == 0)
        continue;
      if (!have_stokes) {
        sb = pol_at().in_B(b.c, b.s);
        have_stokes = true;
      }
      double t = (1 - beta) * (1 + beta);
      double mup = std::clamp((b.mu - beta) / (1 - beta * b.mu), -1.0, 1.0);
      // Equivalent to (1-beta*mu)^2(1-beta^2)/abs(mu-beta),
      // but avoids cancellation at coalescing roots.
      double jacobian = x * x * t * std::sqrt(t) / std::sqrt(D);
      double w = comp.current_fraction * (bfield.p + 1) * pi * b.twist *
                 comp.dist.inv_abs_b_mean / b.r * f * jacobian *
                 overlap(sb.x, sb.z, mup, comp.charge_sign);
      if (!std::isfinite(w) || w < 0)
        throw std::runtime_error("Invalid differential optical depth");
      result[j] = {j, beta, w};
    }
    return result;
  }

  bool
  may_resonate() const {
    if (comp.current_fraction == 0)
      return false;
    if (!options.resonance_bounds)
      return true;
    double beta = std::max(std::abs(comp.dist.b_min), std::abs(comp.dist.b_max));
    if (!(beta < 1 && geodesic.rs >= 0 && bfield.B_pole > 0 && bfield.p > 0 &&
          comp.charge_mass_ratio > 0 && omega_inf > 0 &&
          std::isfinite(bfield.angular_max)))
      return true;
    auto [rmin, rmax] = geodesic.radius_bounds();
    if (!(rmin > geodesic.rs && rmax >= rmin))
      return true;
    double scale = 0.5 * bfield.B_pole * B_to_omega *
                   comp.charge_mass_ratio / omega_inf;
    double xmin = scale * bfield.angular_min *
                  std::pow(bfield.R_star / rmax, 2 + bfield.p) *
                  std::sqrt(1 - geodesic.rs / rmin);
    double xmax = scale * bfield.angular_max *
                  std::pow(bfield.R_star / rmin, 2 + bfield.p) *
                  std::sqrt(1 - geodesic.rs / rmax);
    if (!std::isfinite(xmin) || !std::isfinite(xmax))
      return true;
    // Resonance requires x = gamma*(1-beta*mu), hence 1/edge <= x <= edge.
    double edge = std::sqrt((1 + beta) / (1 - beta));
    return !(xmin > edge * (1 + 1e-10) || xmax < (1 - 1e-10) / edge);
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

  std::vector<double>
  integration_cuts(double left, double right) const {
    std::vector<double> cuts{left, std::midpoint(left, right), right};
    struct ScanPoint {
      double l, x, mu;
    };
    auto sample = [&](double l) {
      if (options.cache_tolerance > 0) {
        auto b = transport(l);
        return ScanPoint{l, comp.charge_mass_ratio * b.x_e, b.mu};
      }
      auto b = get_basis(l);
      return ScanPoint{l, comp.charge_mass_ratio * b.x_e, b.mu};
    };
    // All boundary functions share this grid. Cache only the geometry they
    // need; root-finding trial points still evaluate the actual trajectory.
    std::vector<ScanPoint> grid;
    grid.reserve(options.n_scan + 1);
    grid.push_back(sample(left));
    for (int j = 1; j <= options.n_scan; ++j)
      grid.push_back(
          sample(std::lerp(left, right, double(j) / options.n_scan)));

    auto scan = [&](auto const& value) {
      auto f = [&](double l) {
        auto point = sample(l);
        return value(point.x, point.mu);
      };
      double a = grid.front().l, fa = value(grid.front().x, grid.front().mu);
      for (size_t j = 1; j < grid.size(); ++j) {
        double b = grid[j].l, fb = value(grid[j].x, grid[j].mu);
        if (fa == 0)
          cuts.push_back(a);
        if (fb == 0)
          cuts.push_back(b);
        if ((fa < 0 && fb > 0) || (fa > 0 && fb < 0))
          cuts.push_back(ridders_root(f, a, b, fa, fb));
        a = b;
        fa = fb;
      }
    };

    scan([](double x, double mu) { return discriminant(x, mu); });
    for (double beta : knots) {
      double inv_gamma = std::sqrt((1 - beta) * (1 + beta));
      scan([&](double x, double mu) { return x * inv_gamma + beta * mu - 1; });
    }

    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    return cuts;
  }

  template <class PolarizationAt>
  bool
  integrate_optical_depth(double left, double right,
                          PolarizationAt const& pol_at) {
    auto cuts = integration_cuts(left, right);
    auto integrand = [&](double l) {
      if (options.cache_tolerance > 0)
        return total_rate(rates_at(transport(l), [&] { return pol_at(l); }));
      return total_rate(rates_at(get_basis(l), [&] { return pol_at(l); }));
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
      if (dtau > 0 && tau + dtau > 0) {
        // Find the point where the optical depth reaches zero.
        auto residual = [&](double l) { return tau + integral(l); };
        double l = ridders_root(residual, a, b, tau, tau + dtau, options.event_rtol * (b - a));
        Polarization incoming = pol_at(l);
        Rate event = select_rate(rates(get_basis(l), incoming));
        scattered_data = {geodesic.interp(l), event.beta};
        pol = incoming;
        return true;
      }
      tau += dtau;
    }
    return false;
  }

  bool
  resonance_escape(double l) const {
    auto y = geodesic.interp(l);
    if (std::cos(y[2]) <= 0)
      return false; // Infinity reconstruction is only for the outgoing branch.
    auto b = basis(y);
    double h = 1e-3 * geodesic.stepper.h_old;
    auto minus = get_basis(l - h);
    double dx = (b.x_e - minus.x_e) / (h);
    double dmu = (b.mu - minus.mu) / (h);
    double x = comp.charge_mass_ratio * b.x_e;
    double tx = comp.charge_mass_ratio * dx;
    return discriminant(x, b.mu) < 0 && tx < 0 &&
           std::abs(b.mu - x * dmu / tx) < 1;
  }

  void
  record_escape(double l, Geodesic::YVector y) {
    auto [r, psi, alpha] = y;
    double impact = r * std::sin(alpha) / std::sqrt(1 - geodesic.rs / r);
    // u = r/rho maps [r,infinity) to [1,0], with a regular integrand at 0.
    auto integrand = [&](double u) {
      double t = impact / r;
      return t / std::sqrt(1 - t * t * u * u * (1 - geodesic.rs * u / r));
    };
    double psi_inf = psi + quad.qags(integrand, 0, 1, 1e-11, 1e-11);
    double3 k = to_unit(std::cos(psi_inf) * e1 + std::sin(psi_inf) * e2);
    Stokes sky = pol.observer_stokes(k, n);
    escaped_data = {omega_inf, std::clamp(k.z, -1.0, 1.0), sky.x, sky.y, sky.z};
  }

  enum class EvolveResult { Absorbed, Escaped, Scattered };

  struct EscapedData {
    double omega_inf, muk;
    double Q, U, V;
  } escaped_data;

  struct ScatteredData {
    Geodesic::YVector r_psi_alpha;
    double beta;
  } scattered_data;

  EvolveResult
  evolve() {
    while (true) {
      bool absorbed = geodesic.step();
      cache.leaves.clear();
      bool possible_scattering = may_resonate();
      if (pol.is_adiabatic) {
        auto f_couple = [&](double l) {
          Basis basis = get_basis(l);
          return options.eta_couple * basis.kappa * basis.r - 1;
        };
        if (f_couple(geodesic.l_old()) < 0) {
          Basis b = get_basis(geodesic.l_old());
          pol.set_stokes(b.c, b.s);
        }
      }
      double left = geodesic.l_old(), right = geodesic.l_new();
      if (options.cache_tolerance > 0 && (possible_scattering || !pol.is_adiabatic)) {
        cache.tolerance = options.cache_tolerance;
        cache.build(left, right, [&](double l) { return sample_transport(l); });
      }
      if (pol.is_adiabatic) {
        if (possible_scattering &&
            integrate_optical_depth(left, right, [&](double) { return pol; }))
          return EvolveResult::Scattered;
      } else {
        StepperMagnus stepper;
        stepper.init(geodesic.l_old(), geodesic.l_new(),
                     0.1 * (geodesic.l_new() - geodesic.l_old()), pol.stokes);
        auto omega = [&](double l) {
          return options.cache_tolerance > 0 ? transport(l).omega : get_basis(l).omega;
        };
        while (stepper.do_step(omega)) {
          if (possible_scattering && integrate_optical_depth(
                  stepper.x_old, stepper.x_new,
                  [&](double l) { return Polarization(stepper.dense_out(l)); }))
            return EvolveResult::Scattered;

          pol.set_stokes(stepper.y_new);
          auto y = geodesic.interp(stepper.x_new);
          double change = (stepper.y_new + (-1.0) * stepper.y_old).length() *
                          y[0] / (stepper.x_new - stepper.x_old);
          if (change < options.escape_tolerance &&
                  resonance_escape(stepper.x_new) ||
              (y[0] > 1e4 * geodesic.R_star && std::cos(y[2]) > 0)) {
            record_escape(stepper.x_new, y);
            return EvolveResult::Escaped;
          }
          stepper.update_old();
        }
      }
      if (absorbed || geodesic.stepper.y_new[0] <= geodesic.R_star)
        return EvolveResult::Absorbed;
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
    init(InitializedPhoton{.n = normal,
                           .e1 = b.r_hat,
                           .e2 = cross(normal, b.r_hat),
                           .omega_inf = energy,
                           .r0 = b.r,
                           .psi0 = 0,
                           .alpha0 = std::atan2(cross_len, dot(k, b.r_hat)),
                           .pol = mode});
  }
};
