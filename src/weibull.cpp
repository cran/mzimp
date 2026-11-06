#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

#include <cmath>
#include <vector>
#include <algorithm>

// ============================================================
// Utilities
// ============================================================

inline double clamp_u_weibull(double u, double eps_u) {
  if (u < eps_u) return eps_u;
  if (u > 1.0 - eps_u) return 1.0 - eps_u;
  return u;
}


// ============================================================
// Weibull moment-ratio function
//
// Let a = 1 / shape. Then
//
//   E(X)   = scale * Gamma(1 + a)
//   Var(X) = E(X)^2 *
//            {Gamma(1 + 2a) / Gamma(1 + a)^2 - 1}
//
// Define
//
//   g(a) = log Gamma(1 + 2a)
//          - 2 log Gamma(1 + a).
//
// We solve
//
//   g(a) = log(1 + target_cv2),
//
// where
//
//   target_cv2 = Var(X) / E(X)^2
//              = phi * mu^(var_power - 2).
// ============================================================

inline double weibull_log_moment_ratio(double a) {
  return R::lgammafn(1.0 + 2.0 * a)
    - 2.0 * R::lgammafn(1.0 + a);
}


// Derivative:
//
// g'(a) = 2 digamma(1 + 2a) - 2 digamma(1 + a)

inline double weibull_log_moment_ratio_deriv(double a) {
  return 2.0 * (
    R::digamma(1.0 + 2.0 * a)
      - R::digamma(1.0 + a)
  );
}


// ============================================================
// Solve for a = 1 / Weibull shape
//
// Uses safeguarded Newton iterations:
// - Newton step when it remains inside the current bracket;
// - bisection otherwise.
// ============================================================

inline double solve_weibull_a(double target_cv2,
                              double tol = 1e-12,
                              int maxit = 100) {

  if (!R_finite(target_cv2) || target_cv2 <= 0.0) {
    stop("The target squared coefficient of variation must be positive.");
  }

  const double target_log_ratio = std::log1p(target_cv2);

  // ----------------------------------------------------------
  // Construct a bracket [lower, upper]
  // ----------------------------------------------------------

  double lower = 0.0;
  double upper = 1.0;

  double g_upper = weibull_log_moment_ratio(upper);

  while (g_upper < target_log_ratio) {
    upper *= 2.0;

    if (!R_finite(upper) || upper > 1e6) {
      stop("Unable to bracket the Weibull shape parameter.");
    }

    g_upper = weibull_log_moment_ratio(upper);

    if (!R_finite(g_upper)) {
      stop("Non-finite value encountered while solving Weibull parameters.");
    }
  }

  // A reasonable initial value inside the bracket
  double a = 0.5 * (lower + upper);

  // ----------------------------------------------------------
  // Safeguarded Newton iterations
  // ----------------------------------------------------------

  for (int iter = 0; iter < maxit; ++iter) {

    double g = weibull_log_moment_ratio(a);
    double diff = g - target_log_ratio;

    if (std::abs(diff) <
        tol * (1.0 + std::abs(target_log_ratio))) {
      return a;
    }

    // Update the bracket
    if (diff < 0.0) {
      lower = a;
    } else {
      upper = a;
    }

    double gp = weibull_log_moment_ratio_deriv(a);

    double a_new;

    if (R_finite(gp) && gp > 0.0) {
      a_new = a - diff / gp;
    } else {
      a_new = 0.5 * (lower + upper);
    }

    // Safeguard the Newton step
    if (!R_finite(a_new) ||
        a_new <= lower ||
        a_new >= upper) {
      a_new = 0.5 * (lower + upper);
    }

    if (std::abs(a_new - a) <
        tol * (1.0 + std::abs(a))) {
      return a_new;
    }

    a = a_new;
  }

  return a;
}


// ============================================================
// Convert (mu, phi, var_power) to Weibull parameters
//
// Parameterization:
//
//   E(X)   = mu
//   Var(X) = phi * mu^var_power
//
// Conventional Weibull parameterization:
//
//   F(x) = 1 - exp{-(x / scale)^shape}
//
// We use
//
//   a     = 1 / shape
//   scale = mu / Gamma(1 + a).
// ============================================================

inline void weibull_mu_phi_parameters(double mu,
                                      double phi,
                                      int var_power,
                                      double& a,
                                      double& shape,
                                      double& scale,
                                      double& log_scale,
                                      double& da_dmu,
                                      double& da_dphi) {

  if (!R_finite(mu) || mu <= 0.0) {
    stop("mu must be positive and finite.");
  }

  if (!R_finite(phi) || phi <= 0.0) {
    stop("phi must be positive and finite.");
  }

  if (!(var_power == 1 ||
        var_power == 2 ||
        var_power == 3)) {
    stop("var_power must be one of 1, 2, or 3.");
  }

  // Squared coefficient of variation
  double target_cv2 =
    phi * std::pow(mu, static_cast<double>(var_power - 2));

  if (!R_finite(target_cv2) || target_cv2 <= 0.0) {
    stop("Invalid Weibull variance-to-mean ratio.");
  }

  a = solve_weibull_a(target_cv2);

  shape = 1.0 / a;

  double lgamma_1pa = R::lgammafn(1.0 + a);

  log_scale = std::log(mu) - lgamma_1pa;
  scale = std::exp(log_scale);

  // ----------------------------------------------------------
  // Implicit derivatives of a
  //
  // g(a) = log(1 + target_cv2)
  //
  // da/dtheta =
  //   [d target_cv2/dtheta / (1 + target_cv2)] / g'(a)
  // ----------------------------------------------------------

  double gp = weibull_log_moment_ratio_deriv(a);

  if (!R_finite(gp) || gp <= 0.0) {
    stop("Invalid derivative while solving Weibull parameters.");
  }

  double dlogtarget_dmu = 0.0;

  if (var_power != 2) {
    dlogtarget_dmu =
      static_cast<double>(var_power - 2) *
      target_cv2 /
      (mu * (1.0 + target_cv2));
  }

  double dlogtarget_dphi =
    target_cv2 /
    (phi * (1.0 + target_cv2));

  da_dmu = dlogtarget_dmu / gp;
  da_dphi = dlogtarget_dphi / gp;
}


// ============================================================
// Optional exported parameter-conversion function
// ============================================================

// [[Rcpp::export]]
NumericMatrix weibull_mu_phi_parameters_cpp(
    NumericVector mu,
    double phi,
    int var_power = 2) {

  int n = mu.size();

  NumericMatrix out(n, 7);

  colnames(out) = CharacterVector::create(
    "mu",
    "variance",
    "cv2",
    "a",
    "shape",
    "scale",
    "log_scale"
  );

  for (int i = 0; i < n; ++i) {

    double a;
    double shape;
    double scale;
    double log_scale;
    double da_dmu;
    double da_dphi;

    weibull_mu_phi_parameters(
      mu[i],
      phi,
      var_power,
      a,
      shape,
      scale,
      log_scale,
      da_dmu,
      da_dphi
    );

    double variance =
      phi * std::pow(mu[i], static_cast<double>(var_power));

    double cv2 =
      phi * std::pow(mu[i], static_cast<double>(var_power - 2));

    out(i, 0) = mu[i];
    out(i, 1) = variance;
    out(i, 2) = cv2;
    out(i, 3) = a;
    out(i, 4) = shape;
    out(i, 5) = scale;
    out(i, 6) = log_scale;
  }

  return out;
}


// ============================================================
// I_weibull_mu_phi_cpp
//
// Computes
//
// I_i^(m)(z)
// = integral_{(m-1)/d}^{m/d}
//     exp{-Q_i(u)} Q_i(u)^z / Gamma(z+1) du.
//
// Parameterization:
//
//   E(Lambda_i)   = mu_i
//   Var(Lambda_i) = phi * mu_i^var_power.
//
// Weibull quantile:
//
//   Q_i(u)
//   = scale_i [-log(1-u)]^a_i,
//
// where
//
//   a_i     = 1 / shape_i,
//   scale_i = mu_i / Gamma(1 + a_i).
//
// Each checkerboard cell is transformed using
//
//   u = lower + width * (1 - exp(-s)),
//   du = width * exp(-s) ds,
//
// and evaluated by Gauss-Laguerre quadrature.
// ============================================================

// [[Rcpp::export]]
NumericMatrix I_weibull_mu_phi_cpp(
    IntegerVector z,
    NumericVector mu,
    double phi,
    int d = 1,
    int var_power = 2,
    NumericVector lag_nodes = NumericVector::create(),
    NumericVector lag_weights = NumericVector::create(),
    bool log_out = false,
    double eps_u = 1e-15) {

  int n = z.size();
  int K = lag_nodes.size();

  if (mu.size() != n) {
    stop("length(mu) must equal length(z).");
  }

  if (!R_finite(phi) || phi <= 0.0) {
    stop("phi must be positive.");
  }

  if (d <= 0) {
    stop("d must be a positive integer.");
  }

  if (!(var_power == 1 ||
        var_power == 2 ||
        var_power == 3)) {
    stop("var_power must be one of 1, 2, or 3.");
  }

  if (K == 0 || lag_weights.size() != K) {
    stop("lag_nodes and lag_weights must be supplied and have equal lengths.");
  }

  if (!R_finite(eps_u) || eps_u <= 0.0 || eps_u >= 0.5) {
    stop("eps_u must lie in (0, 0.5).");
  }

  for (int k = 0; k < K; ++k) {
    if (!R_finite(lag_nodes[k]) || lag_nodes[k] < 0.0) {
      stop("lag_nodes must be finite and nonnegative.");
    }

    if (!R_finite(lag_weights[k]) || lag_weights[k] <= 0.0) {
      stop("lag_weights must be positive and finite.");
    }
  }

  NumericMatrix out(n, d);

  // ----------------------------------------------------------
  // Precompute cell-specific values:
  //
  // t = -log(1-u)
  // ----------------------------------------------------------

  std::vector< std::vector<double> > log_t_grid(
    d,
    std::vector<double>(K)
  );

  std::vector<double> log_effective_weight(K);

  for (int k = 0; k < K; ++k) {
    log_effective_weight[k] = std::log(lag_weights[k]);
  }

  for (int m = 0; m < d; ++m) {

    double lower =
      static_cast<double>(m) /
      static_cast<double>(d);

    double upper =
      static_cast<double>(m + 1) /
      static_cast<double>(d);

    double width = upper - lower;

    for (int k = 0; k < K; ++k) {

      double s = lag_nodes[k];
      double one_minus_u;

      if (upper < 1.0) {
        one_minus_u =
          (1.0 - upper) + width * std::exp(-s);
      } else {
        one_minus_u = width * std::exp(-s);
      }

      one_minus_u =
        std::max(one_minus_u, eps_u);

      double t = -std::log(one_minus_u);

      // For the first cell and very small quadrature nodes,
      // t can approach zero.
      t = std::max(t, 1e-300);

      log_t_grid[m][k] = std::log(t);
    }
  }

  // ----------------------------------------------------------
  // Main loop
  // ----------------------------------------------------------

  for (int i = 0; i < n; ++i) {

    if (z[i] < 0) {
      stop("z must contain nonnegative integers.");
    }

    if (!R_finite(mu[i]) || mu[i] <= 0.0) {
      stop("mu must contain positive finite values.");
    }

    double a;
    double shape;
    double scale;
    double log_scale;
    double da_dmu;
    double da_dphi;

    weibull_mu_phi_parameters(
      mu[i],
      phi,
      var_power,
      a,
      shape,
      scale,
      log_scale,
      da_dmu,
      da_dphi
    );

    double zi = static_cast<double>(z[i]);
    double lgamma_z1 = R::lgammafn(zi + 1.0);

    for (int m = 0; m < d; ++m) {

      double log_width = -std::log(static_cast<double>(d));

      double max_log_term = R_NegInf;
      std::vector<double> log_terms(K);

      for (int k = 0; k < K; ++k) {

        double logQ =
          log_scale +
          a * log_t_grid[m][k];

        double log_integrand = R_NegInf;

        if (logQ <= 700.0) {

          double Q = std::exp(logQ);

          double zlogQ =
            (z[i] == 0) ? 0.0 : zi * logQ;

          log_integrand =
            zlogQ -
            Q -
            lgamma_z1;
        }

        double log_term =
          log_width +
          log_effective_weight[k] +
          log_integrand;

        log_terms[k] = log_term;

        if (log_term > max_log_term) {
          max_log_term = log_term;
        }
      }

      double log_value = R_NegInf;

      if (R_finite(max_log_term)) {

        double scaled_sum = 0.0;

        for (int k = 0; k < K; ++k) {
          scaled_sum +=
            std::exp(log_terms[k] - max_log_term);
        }

        if (scaled_sum > 0.0) {
          log_value =
            max_log_term +
            std::log(scaled_sum);
        }
      }

      out(i, m) =
        log_out ? log_value : std::exp(log_value);
    }
  }

  return out;
}


// ============================================================
// I_weibull_mu_phi_grad_cpp
//
// Returns:
//   I       : cell probabilities/integrals
//   d_mu    : derivative with respect to mu_i
//   d_phi   : derivative with respect to phi
//
// The derivative is analytic after solving for a_i.
//
// Since
//
//   log Q = log(mu)
//           - log Gamma(1+a)
//           + a log{-log(1-u)},
//
// we have
//
//   d log Q / d theta
//   = 1(theta=mu)/mu
//     + da/dtheta *
//       [log{-log(1-u)} - digamma(1+a)].
//
// For the Poisson kernel
//
//   h(Q) = exp(-Q) Q^z / Gamma(z+1),
//
//   d h / d theta
//   = h(Q) (z-Q) d log Q/d theta.
// ============================================================

// [[Rcpp::export]]
List I_weibull_mu_phi_grad_cpp(
    IntegerVector z,
    NumericVector mu,
    double phi,
    int d = 1,
    int var_power = 2,
    NumericVector lag_nodes = NumericVector::create(),
    NumericVector lag_weights = NumericVector::create(),
    double eps_u = 1e-15) {

  int n = z.size();
  int K = lag_nodes.size();

  if (mu.size() != n) {
    stop("length(mu) must equal length(z).");
  }

  if (!R_finite(phi) || phi <= 0.0) {
    stop("phi must be positive.");
  }

  if (d <= 0) {
    stop("d must be a positive integer.");
  }

  if (!(var_power == 1 ||
        var_power == 2 ||
        var_power == 3)) {
    stop("var_power must be one of 1, 2, or 3.");
  }

  if (K == 0 || lag_weights.size() != K) {
    stop("lag_nodes and lag_weights must be supplied and have equal lengths.");
  }

  if (!R_finite(eps_u) || eps_u <= 0.0 || eps_u >= 0.5) {
    stop("eps_u must lie in (0, 0.5).");
  }

  for (int k = 0; k < K; ++k) {
    if (!R_finite(lag_nodes[k]) || lag_nodes[k] < 0.0) {
      stop("lag_nodes must be finite and nonnegative.");
    }

    if (!R_finite(lag_weights[k]) || lag_weights[k] <= 0.0) {
      stop("lag_weights must be positive and finite.");
    }
  }

  NumericMatrix prob(n, d);
  NumericMatrix d_mu_out(n, d);
  NumericMatrix d_phi_out(n, d);

  // ----------------------------------------------------------
  // Precompute transformed quadrature grid
  // ----------------------------------------------------------

  std::vector< std::vector<double> > log_t_grid(
    d,
    std::vector<double>(K)
  );

  std::vector<double> log_lag_weight(K);

  for (int k = 0; k < K; ++k) {
    log_lag_weight[k] = std::log(lag_weights[k]);
  }

  for (int m = 0; m < d; ++m) {

    double lower =
      static_cast<double>(m) /
      static_cast<double>(d);

    double upper =
      static_cast<double>(m + 1) /
      static_cast<double>(d);

    double width = upper - lower;

    for (int k = 0; k < K; ++k) {

      double s = lag_nodes[k];
      double one_minus_u;

      if (upper < 1.0) {
        one_minus_u =
          (1.0 - upper) +
          width * std::exp(-s);
      } else {
        one_minus_u =
          width * std::exp(-s);
      }

      one_minus_u =
        std::max(one_minus_u, eps_u);

      double t = -std::log(one_minus_u);
      t = std::max(t, 1e-300);

      log_t_grid[m][k] = std::log(t);
    }
  }

  // ----------------------------------------------------------
  // Main loop
  // ----------------------------------------------------------

  for (int i = 0; i < n; ++i) {

    if (z[i] < 0) {
      stop("z must contain nonnegative integers.");
    }

    if (!R_finite(mu[i]) || mu[i] <= 0.0) {
      stop("mu must contain positive finite values.");
    }

    double mui = mu[i];
    double zi = static_cast<double>(z[i]);
    double lgamma_z1 = R::lgammafn(zi + 1.0);

    double a;
    double shape;
    double scale;
    double log_scale;
    double da_dmu;
    double da_dphi;

    weibull_mu_phi_parameters(
      mui,
      phi,
      var_power,
      a,
      shape,
      scale,
      log_scale,
      da_dmu,
      da_dphi
    );

    double digamma_1pa = R::digamma(1.0 + a);

    for (int m = 0; m < d; ++m) {

      double log_width =
        -std::log(static_cast<double>(d));

      double max_log_term = R_NegInf;

      std::vector<double> log_terms(K);
      std::vector<double> score_mu(K);
      std::vector<double> score_phi(K);

      for (int k = 0; k < K; ++k) {

        double log_t = log_t_grid[m][k];

        double logQ =
          log_scale +
          a * log_t;

        double log_integrand = R_NegInf;
        double Q = R_PosInf;

        if (logQ <= 700.0) {

          Q = std::exp(logQ);

          double zlogQ =
            (z[i] == 0) ? 0.0 : zi * logQ;

          log_integrand =
            zlogQ -
            Q -
            lgamma_z1;
        }

        double log_term =
          log_width +
          log_lag_weight[k] +
          log_integrand;

        log_terms[k] = log_term;

        if (log_term > max_log_term) {
          max_log_term = log_term;
        }

        double common_a_term =
          log_t -
          digamma_1pa;

        double dlogQ_dmu =
          1.0 / mui +
          da_dmu * common_a_term;

        double dlogQ_dphi =
          da_dphi * common_a_term;

        if (R_finite(Q)) {
          score_mu[k] =
            (zi - Q) * dlogQ_dmu;

          score_phi[k] =
            (zi - Q) * dlogQ_dphi;
        } else {
          score_mu[k] = 0.0;
          score_phi[k] = 0.0;
        }
      }

      double value = 0.0;
      double derivative_mu = 0.0;
      double derivative_phi = 0.0;

      if (R_finite(max_log_term)) {

        double scaled_value = 0.0;
        double scaled_dmu = 0.0;
        double scaled_dphi = 0.0;

        for (int k = 0; k < K; ++k) {

          double scaled_term =
            std::exp(log_terms[k] - max_log_term);

          scaled_value += scaled_term;
          scaled_dmu += scaled_term * score_mu[k];
          scaled_dphi += scaled_term * score_phi[k];
        }

        double scale_factor = std::exp(max_log_term);

        value = scale_factor * scaled_value;
        derivative_mu = scale_factor * scaled_dmu;
        derivative_phi = scale_factor * scaled_dphi;
      }

      prob(i, m) = value;
      d_mu_out(i, m) = derivative_mu;
      d_phi_out(i, m) = derivative_phi;
    }
  }

  return List::create(
    _["prob"] = prob,
    _["d_mu"] = d_mu_out,
    _["d_phi"] = d_phi_out
  );
}
