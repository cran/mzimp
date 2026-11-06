#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

#include <cmath>
#include <vector>
#include <algorithm>

// ============================================================
// Utilities
// ============================================================

inline double clamp_u_betaprime(double u, double eps_u) {
  if (u < eps_u) return eps_u;
  if (u > 1.0 - eps_u) return 1.0 - eps_u;
  return u;
}

inline double logit_betaprime(double v, double eps_u) {
  v = clamp_u_betaprime(v, eps_u);
  return std::log(v) - std::log1p(-v);
}

// ============================================================
// Parameter conversion
//
// X = beta * V/(1-V), V ~ Beta(a,a), a > 2.
//
// E(X)   = beta*a/(a-1)
// Var(X) = beta^2*a*(2a-1)/[(a-2)(a-1)^2]
//
// Required:
// E(X)=mu, Var(X)=phi*mu^var_power.
//
// cv2 = phi*mu^(var_power-2)
//     = (2a-1)/[a(a-2)]
//
// a = [cv2+1+sqrt(cv2^2+cv2+1)]/cv2 > 2
// beta = mu*(a-1)/a
// ============================================================

inline void betaprime_mu_phi_parameters(
    double mu,
    double phi,
    double var_power,
    double& cv2,
    double& a,
    double& beta,
    double& log_beta,
    double& da_dmu,
    double& da_dphi,
    double& dbeta_dmu,
    double& dbeta_dphi) {

  if (!R_finite(mu) || mu <= 0.0)
    stop("mu must be positive and finite.");
  if (!R_finite(phi) || phi <= 0.0)
    stop("phi must be positive and finite.");
  if (!R_finite(var_power))
    stop("var_power must be finite.");

  double log_cv2 =
    std::log(phi) +
    (var_power - 2.0) * std::log(mu);

  if (!R_finite(log_cv2) || log_cv2 > 700.0)
    stop("phi * mu^(var_power-2) is numerically too large.");

  cv2 = std::exp(log_cv2);

  if (!R_finite(cv2) || cv2 <= 0.0)
    stop("phi * mu^(var_power-2) must be positive and finite.");

  double s = std::sqrt(cv2 * cv2 + cv2 + 1.0);

  a = (cv2 + 1.0 + s) / cv2;

  if (!R_finite(a) || a <= 2.0)
    stop("The implied beta-prime shape a must be greater than 2.");

  beta = mu * (a - 1.0) / a;

  if (!R_finite(beta) || beta <= 0.0)
    stop("The implied beta-prime scale must be positive and finite.");

  log_beta = std::log(beta);

  double dcv2_dmu =
    cv2 * (var_power - 2.0) / mu;

  double dcv2_dphi =
    cv2 / phi;

  // From cv2(a)=(2a-1)/[a(a-2)]:
  double da_dcv2 =
    -a * a * (a - 2.0) * (a - 2.0) /
    (2.0 * (a * a - a + 1.0));

  da_dmu  = da_dcv2 * dcv2_dmu;
  da_dphi = da_dcv2 * dcv2_dphi;

  double dbeta_da = mu / (a * a);

  dbeta_dmu =
    (a - 1.0) / a +
    dbeta_da * da_dmu;

  dbeta_dphi =
    dbeta_da * da_dphi;
}

// [[Rcpp::export]]
NumericMatrix betaprime_mu_phi_parameters_cpp(
    NumericVector mu,
    double phi,
    double var_power = 2.0) {

  int n = mu.size();
  NumericMatrix out(n, 13);

  colnames(out) = CharacterVector::create(
    "mu", "variance", "cv2", "a", "beta", "log_beta",
    "da_dmu", "da_dphi", "dbeta_dmu", "dbeta_dphi",
    "mean_check", "variance_check", "tail_index"
  );

  for (int i = 0; i < n; ++i) {
    double cv2, a, beta, log_beta;
    double da_dmu, da_dphi, dbeta_dmu, dbeta_dphi;

    betaprime_mu_phi_parameters(
      mu[i], phi, var_power,
      cv2, a, beta, log_beta,
      da_dmu, da_dphi, dbeta_dmu, dbeta_dphi
    );

    double variance = phi * std::pow(mu[i], var_power);
    double mean_check = beta * a / (a - 1.0);
    double variance_check =
      beta * beta * a * (2.0 * a - 1.0) /
      ((a - 2.0) * (a - 1.0) * (a - 1.0));

    out(i, 0) = mu[i];
    out(i, 1) = variance;
    out(i, 2) = cv2;
    out(i, 3) = a;
    out(i, 4) = beta;
    out(i, 5) = log_beta;
    out(i, 6) = da_dmu;
    out(i, 7) = da_dphi;
    out(i, 8) = dbeta_dmu;
    out(i, 9) = dbeta_dphi;
    out(i,10) = mean_check;
    out(i,11) = variance_check;
    out(i,12) = a;
  }

  return out;
}

// ============================================================
// Density
// ============================================================

// [[Rcpp::export]]
NumericVector dbetaprime_mu_phi_cpp(
    NumericVector x,
    NumericVector mu,
    double phi,
    double var_power = 2.0,
    bool log_out = false) {

  int n = x.size();

  if (mu.size() == 1 && n > 1) mu = rep(mu, n);
  if (mu.size() != n)
    stop("length(mu) must equal length(x), or mu must be scalar.");

  NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    if (NumericVector::is_na(x[i])) {
      out[i] = NA_REAL;
      continue;
    }

    if (!R_finite(x[i])) {
      if (x[i] == R_PosInf) {
        out[i] = log_out ? R_NegInf : 0.0;
      } else {
        out[i] = NA_REAL;
      }
      continue;
    }

    if (x[i] <= 0.0) {
      out[i] = log_out ? R_NegInf : 0.0;
      continue;
    }

    double cv2, a, beta, log_beta;
    double da_dmu, da_dphi, dbeta_dmu, dbeta_dphi;

    betaprime_mu_phi_parameters(
      mu[i], phi, var_power,
      cv2, a, beta, log_beta,
      da_dmu, da_dphi, dbeta_dmu, dbeta_dphi
    );

    double log_density =
      (a - 1.0) * std::log(x[i]) -
      a * log_beta -
      R::lbeta(a, a) -
      2.0 * a * std::log1p(x[i] / beta);

    out[i] = log_out ? log_density : std::exp(log_density);
  }

  return out;
}

// ============================================================
// CDF
// ============================================================

// [[Rcpp::export]]
NumericVector pbetaprime_mu_phi_cpp(
    NumericVector q,
    NumericVector mu,
    double phi,
    double var_power = 2.0,
    bool lower_tail = true,
    bool log_p = false) {

  int n = q.size();

  if (mu.size() == 1 && n > 1) mu = rep(mu, n);
  if (mu.size() != n)
    stop("length(mu) must equal length(q), or mu must be scalar.");

  NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    if (NumericVector::is_na(q[i])) {
      out[i] = NA_REAL;
      continue;
    }

    if (q[i] <= 0.0) {
      if (lower_tail) out[i] = log_p ? R_NegInf : 0.0;
      else            out[i] = log_p ? 0.0 : 1.0;
      continue;
    }

    if (q[i] == R_PosInf) {
      if (lower_tail) out[i] = log_p ? 0.0 : 1.0;
      else            out[i] = log_p ? R_NegInf : 0.0;
      continue;
    }

    if (!R_finite(q[i])) {
      out[i] = NA_REAL;
      continue;
    }

    double cv2, a, beta, log_beta;
    double da_dmu, da_dphi, dbeta_dmu, dbeta_dphi;

    betaprime_mu_phi_parameters(
      mu[i], phi, var_power,
      cv2, a, beta, log_beta,
      da_dmu, da_dphi, dbeta_dmu, dbeta_dphi
    );

    double y = q[i] / (q[i] + beta);

    out[i] = R::pbeta(y, a, a, lower_tail, log_p);
  }

  return out;
}

// ============================================================
// Quantile
// ============================================================

// [[Rcpp::export]]
NumericVector qbetaprime_mu_phi_cpp(
    NumericVector p,
    NumericVector mu,
    double phi,
    double var_power = 2.0,
    bool lower_tail = true,
    bool log_p = false,
    double eps_u = 1e-15) {

  int n = p.size();

  if (mu.size() == 1 && n > 1) mu = rep(mu, n);
  if (mu.size() != n)
    stop("length(mu) must equal length(p), or mu must be scalar.");

  if (!R_finite(eps_u) || eps_u <= 0.0 || eps_u >= 0.5)
    stop("eps_u must lie in (0,0.5).");

  NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    if (NumericVector::is_na(p[i])) {
      out[i] = NA_REAL;
      continue;
    }

    double prob;

    if (log_p) {
      if (p[i] > 0.0)
        stop("When log_p=TRUE, all probabilities must be <= 0.");
      prob = std::exp(p[i]);
    } else {
      prob = p[i];
    }

    if (!lower_tail) prob = 1.0 - prob;

    if (prob < 0.0 || prob > 1.0)
      stop("Probabilities must lie in [0,1].");

    if (prob == 0.0) {
      out[i] = 0.0;
      continue;
    }

    if (prob == 1.0) {
      out[i] = R_PosInf;
      continue;
    }

    double cv2, a, beta, log_beta;
    double da_dmu, da_dphi, dbeta_dmu, dbeta_dphi;

    betaprime_mu_phi_parameters(
      mu[i], phi, var_power,
      cv2, a, beta, log_beta,
      da_dmu, da_dphi, dbeta_dmu, dbeta_dphi
    );

    double v = R::qbeta(prob, a, a, true, false);
    v = clamp_u_betaprime(v, eps_u);

    out[i] = std::exp(log_beta + logit_betaprime(v, eps_u));
  }

  return out;
}

// ============================================================
// Random generation
// ============================================================

// [[Rcpp::export]]
NumericVector rbetaprime_mu_phi_cpp(
    int n,
    NumericVector mu,
    double phi,
    double var_power = 2.0,
    double eps_u = 1e-15) {

  if (n <= 0) stop("n must be a positive integer.");

  if (mu.size() == 1 && n > 1) mu = rep(mu, n);
  if (mu.size() != n)
    stop("length(mu) must equal n, or mu must be scalar.");

  if (!R_finite(eps_u) || eps_u <= 0.0 || eps_u >= 0.5)
    stop("eps_u must lie in (0,0.5).");

  NumericVector out(n);

  for (int i = 0; i < n; ++i) {
    double cv2, a, beta, log_beta;
    double da_dmu, da_dphi, dbeta_dmu, dbeta_dphi;

    betaprime_mu_phi_parameters(
      mu[i], phi, var_power,
      cv2, a, beta, log_beta,
      da_dmu, da_dphi, dbeta_dmu, dbeta_dphi
    );

    double v = clamp_u_betaprime(R::rbeta(a, a), eps_u);

    out[i] = std::exp(log_beta + logit_betaprime(v, eps_u));
  }

  return out;
}

// ============================================================
// d/da logit{Q_Beta(u;a,a)}
//
// Five-point central finite difference.
// ============================================================

inline double dlogit_qbeta_da_betaprime(
    double u,
    double a,
    double eps_u,
    double fd_rel_step) {

  u = clamp_u_betaprime(u, eps_u);

  if (!R_finite(fd_rel_step) || fd_rel_step <= 0.0)
    stop("fd_rel_step must be positive and finite.");

  double h = fd_rel_step * std::max(1.0, std::fabs(a));
  h = std::min(h, 0.20 * a);

  if (!R_finite(h) || h <= 0.0)
    stop("Invalid finite-difference step.");

  double am2 = a - 2.0 * h;
  double am1 = a - h;
  double ap1 = a + h;
  double ap2 = a + 2.0 * h;

  if (am2 <= 0.0) {
    h = std::min(h, 0.40 * a);
    am1 = a - h;
    ap1 = a + h;

    double gm1 = logit_betaprime(
      R::qbeta(u, am1, am1, true, false), eps_u
    );
    double gp1 = logit_betaprime(
      R::qbeta(u, ap1, ap1, true, false), eps_u
    );

    return (gp1 - gm1) / (2.0 * h);
  }

  double gm2 = logit_betaprime(
    R::qbeta(u, am2, am2, true, false), eps_u
  );
  double gm1 = logit_betaprime(
    R::qbeta(u, am1, am1, true, false), eps_u
  );
  double gp1 = logit_betaprime(
    R::qbeta(u, ap1, ap1, true, false), eps_u
  );
  double gp2 = logit_betaprime(
    R::qbeta(u, ap2, ap2, true, false), eps_u
  );

  return (gm2 - 8.0 * gm1 + 8.0 * gp1 - gp2) / (12.0 * h);
}


// ============================================================
// Fixed 32-point Gauss-Legendre rule on [-1,1].
// Used internally for all beta-prime cell integrals.  The legacy
// lag_nodes/lag_weights arguments are retained only for API compatibility.
// ============================================================
inline void gauss_legendre_32_betaprime(
    std::vector<double>& nodes,
    std::vector<double>& weights) {
  static const double x[16] = {
    0.048307665687738316235, 0.144471961582796493485,
    0.239287362252137074545, 0.331868602282127649780,
    0.421351276130635345364, 0.506899908932229390024,
    0.587715757240762329041, 0.663044266930215200975,
    0.732182118740289680387, 0.794483795967942406963,
    0.849367613732569970134, 0.896321155766052123965,
    0.934906075937739689171, 0.964762255587506430774,
    0.985611511545268335400, 0.997263861849481563545
  };
  static const double w[16] = {
    0.096540088514727800567, 0.095638720079274859419,
    0.093844399080804565639, 0.091173878695763884713,
    0.087652093004403811143, 0.083311924226946755222,
    0.078193895787070306472, 0.072345794108848506225,
    0.065822222776361846838, 0.058684093478535547145,
    0.050998059262376176196, 0.042835898022226680657,
    0.034273862913021433103, 0.025392065309262059456,
    0.016274394730905670605, 0.007018610009470096600
  };
  nodes.resize(32); weights.resize(32);
  for (int j=0;j<16;++j) {
    nodes[15-j] = -x[j]; weights[15-j] = w[j];
    nodes[16+j] =  x[j]; weights[16+j] = w[j];
  }
}

// ============================================================
// Density-based mixed-Poisson cell probabilities
//
// Let T ~ Beta(a,a) and Lambda = beta * T/(1-T).
// For cell m, define probability boundaries p_m=m/d and
// t_m = qbeta(p_m; a,a). Then
//
// I_i^(m)(z)
// = integral_{t_{m-1}}^{t_m}
//     Pois(z; beta_i*t/(1-t)) * f_Beta(t;a_i,a_i) dt.
//
// Therefore qbeta is used only for the d-1 cell boundaries,
// not at every quadrature node. For d=1 no qbeta call is needed.
// Each finite cell [lower, upper] is integrated directly using a
// fixed 32-point Gauss-Legendre rule. If x_k,w_k are the standard
// nodes/weights on [-1,1], then
//
//   t_k = midpoint + half_width*x_k,
//   I ~= half_width * sum_k w_k h(t_k).
//
// No transformation to (0,infinity) is used.
// ============================================================

inline double qbeta_sym_betaprime(double prob, double a, double eps_u) {
  if (prob <= 0.0) return 0.0;
  if (prob >= 1.0) return 1.0;
  double q = R::qbeta(prob, a, a, true, false);
  return clamp_u_betaprime(q, eps_u);
}

// Five-point derivative of qbeta(p; a,a) with respect to a.
// This is evaluated only at cell boundaries, rather than at every node.
inline double dqbeta_sym_da_betaprime(
    double prob,
    double a,
    double eps_u,
    double fd_rel_step) {

  if (prob <= 0.0 || prob >= 1.0) return 0.0;
  if (!R_finite(fd_rel_step) || fd_rel_step <= 0.0)
    stop("fd_rel_step must be positive and finite.");

  double h = fd_rel_step * std::max(1.0, std::fabs(a));
  h = std::min(h, 0.20 * a);

  if (!R_finite(h) || h <= 0.0)
    stop("Invalid finite-difference step for qbeta boundary derivative.");

  if (a - 2.0 * h <= 0.0) {
    h = std::min(h, 0.40 * a);
    double qm1 = qbeta_sym_betaprime(prob, a - h, eps_u);
    double qp1 = qbeta_sym_betaprime(prob, a + h, eps_u);
    return (qp1 - qm1) / (2.0 * h);
  }

  double qm2 = qbeta_sym_betaprime(prob, a - 2.0 * h, eps_u);
  double qm1 = qbeta_sym_betaprime(prob, a - h, eps_u);
  double qp1 = qbeta_sym_betaprime(prob, a + h, eps_u);
  double qp2 = qbeta_sym_betaprime(prob, a + 2.0 * h, eps_u);

  return (qm2 - 8.0 * qm1 + 8.0 * qp1 - qp2) / (12.0 * h);
}

inline void betaprime_cell_boundaries(
    double a,
    int d,
    double eps_u,
    bool need_derivative,
    double fd_rel_step,
    std::vector<double>& bounds,
    std::vector<double>& dbounds_da) {

  bounds.assign(d + 1, 0.0);
  dbounds_da.assign(d + 1, 0.0);

  bounds[0] = 0.0;
  bounds[d] = 1.0;

  if (d == 1) return;

  for (int m = 1; m < d; ++m) {
    double prob = static_cast<double>(m) / static_cast<double>(d);
    bounds[m] = qbeta_sym_betaprime(prob, a, eps_u);

    if (need_derivative) {
      dbounds_da[m] = dqbeta_sym_da_betaprime(
        prob, a, eps_u, fd_rel_step
      );
    }
  }
}

inline double betaprime_log_integrand_t(
    double t,
    int z,
    double a,
    double beta,
    double log_beta,
    double lgamma_z1,
    double eps_u,
    double& lambda) {

  t = clamp_u_betaprime(t, eps_u);

  double log_t = std::log(t);
  double log_1mt = std::log1p(-t);
  double log_lambda = log_beta + log_t - log_1mt;

  if (!R_finite(log_lambda) || log_lambda > 700.0) {
    lambda = R_PosInf;
    return R_NegInf;
  }

  lambda = std::exp(log_lambda);

  double log_pois =
    ((z == 0) ? 0.0 : static_cast<double>(z) * log_lambda) -
    lambda - lgamma_z1;

  double log_beta_density =
    (a - 1.0) * (log_t + log_1mt) - R::lbeta(a, a);

  return log_pois + log_beta_density;
}

// [[Rcpp::export]]
NumericMatrix I_betaprime_mu_phi_cpp(
    IntegerVector z,
    NumericVector mu,
    double phi,
    int d = 1,
    NumericVector lag_nodes = NumericVector::create(),
    NumericVector lag_weights = NumericVector::create(),
    double var_power = 2.0,
    bool log_out = false,
    double eps_u = 1e-15) {

  int n = z.size();
  const int K = 32;

  if (mu.size() != n)
    stop("length(mu) must equal length(z).");
  if (d <= 0)
    stop("d must be a positive integer.");
  if (!R_finite(eps_u) || eps_u <= 0.0 || eps_u >= 0.5)
    stop("eps_u must lie in (0,0.5).");

  // Legacy lag_nodes/lag_weights arguments are intentionally ignored.
  std::vector<double> gl_nodes, gl_weights;
  gauss_legendre_32_betaprime(gl_nodes, gl_weights);

  NumericMatrix out(n, d);

  bool common_shape = std::fabs(var_power - 2.0) < 1e-14;
  std::vector<double> cached_bounds, dummy_deriv;

  if (common_shape && n > 0) {
    double cv2, a, beta, log_beta;
    double da_dmu, da_dphi, dbeta_dmu, dbeta_dphi;
    betaprime_mu_phi_parameters(
      mu[0], phi, var_power,
      cv2, a, beta, log_beta,
      da_dmu, da_dphi, dbeta_dmu, dbeta_dphi
    );
    betaprime_cell_boundaries(
      a, d, eps_u, false, 1e-5, cached_bounds, dummy_deriv
    );
  }

  for (int i = 0; i < n; ++i) {
    if (z[i] < 0)
      stop("z must contain nonnegative integers.");
    if (!R_finite(mu[i]) || mu[i] <= 0.0)
      stop("mu must contain positive finite values.");

    double cv2, a, beta, log_beta;
    double da_dmu, da_dphi, dbeta_dmu, dbeta_dphi;

    betaprime_mu_phi_parameters(
      mu[i], phi, var_power,
      cv2, a, beta, log_beta,
      da_dmu, da_dphi, dbeta_dmu, dbeta_dphi
    );

    std::vector<double> bounds, dummy;
    if (common_shape) {
      bounds = cached_bounds;
    } else {
      betaprime_cell_boundaries(
        a, d, eps_u, false, 1e-5, bounds, dummy
      );
    }

    double lgamma_z1 = R::lgammafn(static_cast<double>(z[i]) + 1.0);

    for (int m = 0; m < d; ++m) {
      double lower = bounds[m];
      double upper = bounds[m + 1];
      double width = upper - lower;

      if (!(width > 0.0) || !R_finite(width))
        stop("Invalid beta-prime cell width.");

      double half_width = 0.5 * width;
      double midpoint = 0.5 * (lower + upper);
      double log_half_width = std::log(half_width);
      double max_log_term = R_NegInf;
      std::vector<double> log_terms(K);

      for (int k = 0; k < K; ++k) {
        double t = midpoint + half_width * gl_nodes[k];
        double lambda;
        double log_h = betaprime_log_integrand_t(
          t, z[i], a, beta, log_beta, lgamma_z1, eps_u, lambda
        );

        double log_term = log_half_width + std::log(gl_weights[k]) + log_h;
        log_terms[k] = log_term;
        if (log_term > max_log_term) max_log_term = log_term;
      }

      double log_value = R_NegInf;

      if (R_finite(max_log_term)) {
        double scaled_sum = 0.0;
        for (int k = 0; k < K; ++k)
          scaled_sum += std::exp(log_terms[k] - max_log_term);

        if (scaled_sum > 0.0)
          log_value = max_log_term + std::log(scaled_sum);
      }

      out(i, m) = log_out ? log_value : std::exp(log_value);
    }
  }

  return out;
}

// ============================================================
// Analytic gradient of the density integral (32-point Gauss-Legendre).
//
// Only derivatives of the beta cell boundaries with respect to a
// are evaluated numerically. For d=1 these terms vanish exactly and
// no qbeta call is made.
// ============================================================

// [[Rcpp::export]]
List I_betaprime_mu_phi_grad_cpp(
    IntegerVector z,
    NumericVector mu,
    double phi,
    int d = 1,
    NumericVector lag_nodes = NumericVector::create(),
    NumericVector lag_weights = NumericVector::create(),
    double var_power = 2.0,
    double eps_u = 1e-15,
    double fd_rel_step = 1e-5) {

  int n = z.size();
  const int K = 32;

  if (mu.size() != n)
    stop("length(mu) must equal length(z).");
  if (d <= 0)
    stop("d must be a positive integer.");
  if (!R_finite(eps_u) || eps_u <= 0.0 || eps_u >= 0.5)
    stop("eps_u must lie in (0,0.5).");
  if (!R_finite(fd_rel_step) || fd_rel_step <= 0.0)
    stop("fd_rel_step must be positive and finite.");

  // Legacy lag_nodes/lag_weights arguments are intentionally ignored.
  std::vector<double> gl_nodes, gl_weights;
  gauss_legendre_32_betaprime(gl_nodes, gl_weights);

  NumericMatrix I_out(n, d);
  NumericMatrix d_mu_out(n, d);
  NumericMatrix d_phi_out(n, d);

  bool common_shape = std::fabs(var_power - 2.0) < 1e-14;
  std::vector<double> cached_bounds, cached_dbounds;

  if (common_shape && n > 0) {
    double cv2, a, beta, log_beta;
    double da_dmu, da_dphi, dbeta_dmu, dbeta_dphi;
    betaprime_mu_phi_parameters(
      mu[0], phi, var_power,
      cv2, a, beta, log_beta,
      da_dmu, da_dphi, dbeta_dmu, dbeta_dphi
    );

    betaprime_cell_boundaries(
      a, d, eps_u, true, fd_rel_step,
      cached_bounds, cached_dbounds
    );
  }

  for (int i = 0; i < n; ++i) {
    if (z[i] < 0)
      stop("z must contain nonnegative integers.");
    if (!R_finite(mu[i]) || mu[i] <= 0.0)
      stop("mu must contain positive finite values.");

    double cv2, a, beta, log_beta;
    double da_dmu, da_dphi, dbeta_dmu, dbeta_dphi;

    betaprime_mu_phi_parameters(
      mu[i], phi, var_power,
      cv2, a, beta, log_beta,
      da_dmu, da_dphi, dbeta_dmu, dbeta_dphi
    );

    double dlogbeta_dmu = dbeta_dmu / beta;
    double dlogbeta_dphi = dbeta_dphi / beta;

    std::vector<double> bounds, dbounds_da;
    if (common_shape) {
      bounds = cached_bounds;
      dbounds_da = cached_dbounds;
    } else {
      betaprime_cell_boundaries(
        a, d, eps_u, true, fd_rel_step,
        bounds, dbounds_da
      );
    }

    double lgamma_z1 = R::lgammafn(static_cast<double>(z[i]) + 1.0);
    double shape_normalizer_score =
      -2.0 * R::digamma(a) + 2.0 * R::digamma(2.0 * a);

    for (int m = 0; m < d; ++m) {
      double lower = bounds[m];
      double upper = bounds[m + 1];
      double width = upper - lower;

      double dlower_da = dbounds_da[m];
      double dupper_da = dbounds_da[m + 1];
      double dwidth_da = dupper_da - dlower_da;

      if (!(width > 0.0) || !R_finite(width))
        stop("Invalid beta-prime cell width.");

      double half_width = 0.5 * width;
      double midpoint = 0.5 * (lower + upper);
      double dmidpoint_da = 0.5 * (dlower_da + dupper_da);
      double dhalfwidth_da = 0.5 * dwidth_da;
      double log_half_width = std::log(half_width);
      double dlogwidth_da = dwidth_da / width;

      double max_log_term = R_NegInf;
      std::vector<double> log_terms(K);
      std::vector<double> score_mu(K);
      std::vector<double> score_phi(K);

      for (int k = 0; k < K; ++k) {
        double xk = gl_nodes[k];
        double t = midpoint + half_width * xk;
        t = clamp_u_betaprime(t, eps_u);

        double dt_da = dmidpoint_da + dhalfwidth_da * xk;

        double lambda;
        double log_h = betaprime_log_integrand_t(
          t, z[i], a, beta, log_beta, lgamma_z1, eps_u, lambda
        );

        double log_term = log_half_width + std::log(gl_weights[k]) + log_h;
        log_terms[k] = log_term;
        if (log_term > max_log_term) max_log_term = log_term;

        if (!R_finite(log_h) || !R_finite(lambda)) {
          score_mu[k] = 0.0;
          score_phi[k] = 0.0;
          continue;
        }

        double log_t = std::log(t);
        double log_1mt = std::log1p(-t);

        // Partial derivative of log h with respect to a, holding t fixed.
        double dlogh_da_fixed =
          log_t + log_1mt + shape_normalizer_score;

        // Partial derivative of log h with respect to t.
        double inv_t1mt = 1.0 / (t * (1.0 - t));
        double dlogh_dt =
          (static_cast<double>(z[i]) - lambda) * inv_t1mt +
          (a - 1.0) * (1.0 - 2.0 * t) * inv_t1mt;

        double common_a_score =
          dlogh_da_fixed + dlogwidth_da + dlogh_dt * dt_da;

        score_mu[k] =
          (static_cast<double>(z[i]) - lambda) * dlogbeta_dmu +
          common_a_score * da_dmu;

        score_phi[k] =
          (static_cast<double>(z[i]) - lambda) * dlogbeta_dphi +
          common_a_score * da_dphi;
      }

      double value = 0.0;
      double derivative_mu = 0.0;
      double derivative_phi = 0.0;

      if (R_finite(max_log_term)) {
        double scaled_value = 0.0;
        double scaled_dmu = 0.0;
        double scaled_dphi = 0.0;

        for (int k = 0; k < K; ++k) {
          double scaled_term = std::exp(log_terms[k] - max_log_term);
          scaled_value += scaled_term;
          scaled_dmu += scaled_term * score_mu[k];
          scaled_dphi += scaled_term * score_phi[k];
        }

        double scale_factor = std::exp(max_log_term);
        value = scale_factor * scaled_value;
        derivative_mu = scale_factor * scaled_dmu;
        derivative_phi = scale_factor * scaled_dphi;
      }

      I_out(i, m) = value;
      d_mu_out(i, m) = derivative_mu;
      d_phi_out(i, m) = derivative_phi;
    }
  }

  return List::create(
    _["I"] = I_out,
    _["value"] = I_out,
    _["d_mu"] = d_mu_out,
    _["d_phi"] = d_phi_out
  );
}
