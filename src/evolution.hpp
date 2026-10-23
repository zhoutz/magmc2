#pragma once

#include "bfield.hpp"
#include "find_root.hpp"
#include "geodesic.hpp"
#include "initialize.hpp"
#include "magnus.hpp"
#include "quad.hpp"
#include "ran.hpp"

template <PhotonInitializer Init, typename Comp>
struct PhotonEvolution {
  struct Options {
    double eta_couple = 1e-3, escape_tolerance = 1e-3;
    double quad_atol = 1e-8, quad_rtol = 1e-8;
    int n_knots = 8, n_scan = 8;
  } options;
  BField const& bfield;
  Init const& initializer;
  Comp comp;
  Ran ran;
  Geodesic geodesic;
  Quad quad;

  std::vector<double> knots;

  double3 n, e1, e2;
  double omega_inf, tau;

  bool is_adiabatic;
  AdiabaticPolarization apol;
  Stokes stokes;

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

  Stokes
  eigenmode(Basis const& b) const {
    double sign = apol == AdiabaticPolarization::O ? 1 : -1;
    return {sign * b.c, sign * b.s, 0};
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
    tau = std::log(ran.U());
    is_adiabatic = true;
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
  rates(Basis const& b, Stokes S, bool adiabatic) const {
    Rates result{};
    double qb = adiabatic ? (apol == AdiabaticPolarization::O ? 1 : -1)
                          : b.c * S.x + b.s * S.y;
    double v = adiabatic ? 0 : S.z;

    if (comp.current_fraction == 0 || b.twist == 0)
      return result;
    double x = comp.charge_mass_ratio * b.x_e, D = discriminant(x, b.mu);
    // Exact coalescence is an integration boundary, not a sample.
    if (!(D > 0) || !(x > 0))
      return result;
    auto roots = resonance_roots(x, b.mu, D);
    for (int j = 0; j < 2; ++j) {
      double beta = roots[j];
      double f = comp.dist.f(comp.charge_sign * beta);
      if (!(beta > -1 && beta < 1) || f == 0)
        continue;
      double t = (1 - beta) * (1 + beta);
      double mup = std::clamp((b.mu - beta) / (1 - beta * b.mu), -1.0, 1.0);
      // Equivalent to (1-beta*mu)^2(1-beta^2)/abs(mu-beta),
      // but avoids cancellation at coalescing roots.
      double jacobian = x * x * t * std::sqrt(t) / std::sqrt(D);
      double w = comp.current_fraction * (bfield.p + 1) * pi * b.twist *
                 comp.dist.inv_abs_b_mean / b.r * f * jacobian *
                 overlap(qb, v, mup, comp.charge_sign);
      if (!std::isfinite(w) || w < 0)
        throw std::runtime_error("Invalid differential optical depth");
      result[j] = {j, beta, w};
    }
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

  std::vector<double>
  integration_cuts(double left, double right) const {
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
          cuts.push_back(ridders_root(f, a, b, fa, fb));
        a = b;
        fa = fb;
      }
    };

    scan([&](double l) {
      auto b = get_basis(l);
      return discriminant(comp.charge_mass_ratio * b.x_e, b.mu);
    });
    for (double beta : knots) {
      double inv_gamma = std::sqrt((1 - beta) * (1 + beta));
      scan([&](double l) {
        auto b = get_basis(l);
        return comp.charge_mass_ratio * b.x_e * inv_gamma + beta * b.mu - 1;
      });
    }

    std::sort(cuts.begin(), cuts.end());
    cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
    return cuts;
  }

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
      if (dtau > 0 && tau + dtau > 0) {
        // Find the point where the optical depth reaches zero.
        auto residual = [&](double l) { return tau + integral(l); };
        double l = ridders_root(residual, a, b, tau, tau + dtau);
        Stokes incoming = pol(l);
        Rate event = select_rate(rates(get_basis(l), incoming, is_adiabatic));
        scattered_data = {geodesic.interp(l), event.beta};
        stokes = incoming;
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
    Stokes sky = observer_stokes(stokes, k, n);
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
      if (is_adiabatic) {
        auto f_couple = [&](double l) {
          Basis basis = get_basis(l);
          return options.eta_couple * basis.kappa * basis.r - 1;
        };
        if (f_couple(geodesic.l_old()) < 0) {
          stokes = eigenmode(get_basis(geodesic.l_old()));
          is_adiabatic = false;
        }
      }
      double left = geodesic.l_old(), right = geodesic.l_new();
      if (is_adiabatic) {
        if (integrate_scattering(left, right,
                                 [&](double) { return Stokes{0, 0, 0}; }))
          return EvolveResult::Scattered;
      } else {
        StepperMagnus stepper;
        stepper.init(geodesic.l_old(), geodesic.l_new(),
                     0.1 * (geodesic.l_new() - geodesic.l_old()), stokes);
        auto omega = [&](double l) { return get_basis(l).omega; };
        while (stepper.do_step(omega)) {
          if (integrate_scattering(stepper.x_old, stepper.x_new, [&](double l) {
                return stepper.dense_out(l);
              }))
            return EvolveResult::Scattered;

          stokes = stepper.y_new;
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
