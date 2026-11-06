#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

#include <cmath>
#include <vector>
#include <algorithm>
#include <limits>

// ============================================================
// Utilities
// ============================================================

inline double clamp_u_loglogistic(double u, double eps_u) {
  if (u < eps_u) return eps_u;
  if (u > 1.0 - eps_u) return 1.0 - eps_u;
  return u;
}


// Stable evaluation of log(sin(x)) for x in (0, pi).
inline double log_sin_positive(double x) {

  double sx = std::sin(x);

  if (!R_finite(sx) || sx <= 0.0) {
    stop("Invalid argument encountered in log_sin_positive().");
  }

  return std::log(sx);
}


// ============================================================
// Log-logistic moment-ratio function
//
// Let
//
//   a = 1 / kappa,
//   0 < a < 1/2.
//
// For a log-logistic random variable,
// provided r * a < 1,
//
//   E(X^r)
//   = beta^r * (r * pi * a) / sin(r * pi * a).
//
// Therefore,
//
//   E(X) = beta * A1(a),
//
// where
//
//   A1(a) = pi * a / sin(pi * a),
//
// and
//
//   E(X^2) = beta^2 * A2(a),
//
// where
//
//   A2(a) = 2 * pi * a / sin(2 * pi * a).
//
// Hence,
//
//   1 + CV^2
//   = A2(a) / A1(a)^2.
//
// Define
//
//   g(a)
//   = log A2(a) - 2 log A1(a).
//
// We solve
//
//   g(a) = log(1 + target_cv2),
//
// where
//
//   target_cv2
//   = Var(X) / E(X)^2
//   = phi * mu^(var_power - 2).
// ============================================================

inline double loglogistic_log_A1(double a) {

  if (!R_finite(a) || a <= 0.0 || a >= 0.5) {
    stop("a must lie in (0, 0.5).");
  }

  double x = M_PI * a;

  return std::log(M_PI) +
    std::log(a) -
    log_sin_positive(x);
}


inline double loglogistic_log_A2(double a) {

  if (!R_finite(a) || a <= 0.0 || a >= 0.5) {
    stop("a must lie in (0, 0.5).");
  }

  double x = 2.0 * M_PI * a;

  return std::log(2.0 * M_PI) +
    std::log(a) -
    log_sin_positive(x);
}


inline double loglogistic_log_moment_ratio(double a) {

  return loglogistic_log_A2(a) -
    2.0 * loglogistic_log_A1(a);
}


// ============================================================
// Derivative of the log moment ratio
//
// Since
//
//   log A1(a)
//   = log(pi) + log(a) - log sin(pi a),
//
// we have
//
//   d log A1(a) / da
//   = 1/a - pi cot(pi a).
//
// Similarly,
//
//   d log A2(a) / da
//   = 1/a - 2 pi cot(2 pi a).
//
// Therefore,
//
//   g'(a)
//   = -1/a
//     - 2 pi cot(2 pi a)
//     + 2 pi cot(pi a).
// ============================================================

inline double loglogistic_log_moment_ratio_deriv(double a) {

  if (!R_finite(a) || a <= 0.0 || a >= 0.5) {
    stop("a must lie in (0, 0.5).");
  }

  double x1 = M_PI * a;
  double x2 = 2.0 * M_PI * a;

  double sin_x1 = std::sin(x1);
  double sin_x2 = std::sin(x2);

  if (!R_finite(sin_x1) || sin_x1 <= 0.0 ||
      !R_finite(sin_x2) || sin_x2 <= 0.0) {
    stop("Invalid trigonometric value in moment-ratio derivative.");
  }

  double cot_x1 = std::cos(x1) / sin_x1;
  double cot_x2 = std::cos(x2) / sin_x2;

  return
    -1.0 / a
    - 2.0 * M_PI * cot_x2
    + 2.0 * M_PI * cot_x1;
}


// ============================================================
// Solve for a = 1 / kappa
//
// The admissible interval is
//
//   0 < a < 1/2,
//
// because the log-logistic variance exists only for kappa > 2.
//
// As a -> 0,
//   CV^2 -> 0.
//
// As a -> 1/2,
//   CV^2 -> infinity.
//
// Thus, every positive target_cv2 has a unique solution.
//
// A safeguarded Newton algorithm is used:
//
// - Newton's method when the proposed step remains inside the
//   current bracket;
// - bisection otherwise.
// ============================================================

inline double solve_loglogistic_a(
    double target_cv2,
    double tol = 1e-12,
    int maxit = 100) {

  if (!R_finite(target_cv2) || target_cv2 <= 0.0) {
    stop("The target squared coefficient of variation must be positive.");
  }

  if (!R_finite(tol) || tol <= 0.0) {
    stop("tol must be positive and finite.");
  }

  if (maxit <= 0) {
    stop("maxit must be a positive integer.");
  }

  const double target_log_ratio = std::log1p(target_cv2);

  // Keep the bracket strictly inside (0, 1/2).
  double lower = 1e-10;
  double upper = 0.5 - 1e-10;

  double g_lower =
    loglogistic_log_moment_ratio(lower) -
    target_log_ratio;

  double g_upper =
    loglogistic_log_moment_ratio(upper) -
    target_log_ratio;

  // Numerically, the value near a = 0 can be extremely close
  // to zero. A negative or zero lower residual is expected.
  if (!R_finite(g_lower)) {
    stop("Non-finite lower bracket value while solving log-logistic parameters.");
  }

  if (!R_finite(g_upper) || g_upper <= 0.0) {
    stop("Unable to bracket the log-logistic shape parameter.");
  }

  // A moderate starting value usually gives stable Newton steps.
  double a = 0.25;

  if (a <= lower || a >= upper) {
    a = 0.5 * (lower + upper);
  }

  for (int iter = 0; iter < maxit; ++iter) {

    double g =
      loglogistic_log_moment_ratio(a) -
      target_log_ratio;

    if (!R_finite(g)) {
      a = 0.5 * (lower + upper);
      continue;
    }

    if (std::abs(g) <
        tol * (1.0 + std::abs(target_log_ratio))) {
      return a;
    }

    // Update the bracket.
    if (g < 0.0) {
      lower = a;
    } else {
      upper = a;
    }

    double gp =
      loglogistic_log_moment_ratio_deriv(a);

    double a_new;

    if (R_finite(gp) && gp > 0.0) {
      a_new = a - g / gp;
    } else {
      a_new = 0.5 * (lower + upper);
    }

    // Safeguard the Newton step.
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
// Convert (mu, phi, var_power) to log-logistic parameters
//
// Parameterization:
//
//   E(X)   = mu,
//   Var(X) = phi * mu^var_power.
//
// Define
//
//   target_cv2
//   = phi * mu^(var_power - 2).
//
// Let a = 1/kappa solve
//
//   log A2(a) - 2 log A1(a)
//   = log(1 + target_cv2).
//
// Then
//
//   kappa = 1/a,
//
// and, because
//
//   E(X) = beta * A1(a),
//
// the scale parameter is
//
//   beta
//   = mu / A1(a)
//   = mu * sin(pi a) / (pi a).
//
// The derivatives da/dmu and da/dphi are obtained by
// implicit differentiation.
// ============================================================

inline void loglogistic_mu_phi_parameters(
    double mu,
    double phi,
    int var_power,
    double& a,
    double& kappa,
    double& beta,
    double& log_beta,
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

  // Squared coefficient of variation.
  double target_cv2 =
    phi *
    std::pow(
      mu,
      static_cast<double>(var_power - 2)
    );

  if (!R_finite(target_cv2) || target_cv2 <= 0.0) {
    stop("Invalid log-logistic variance-to-mean ratio.");
  }

  a = solve_loglogistic_a(target_cv2);

  if (!R_finite(a) || a <= 0.0 || a >= 0.5) {
    stop("Invalid solution for the log-logistic shape parameter.");
  }

  kappa = 1.0 / a;

  // log beta
  // = log mu + log sin(pi a) - log(pi a).
  double pi_a = M_PI * a;

  log_beta =
    std::log(mu) +
    log_sin_positive(pi_a) -
    std::log(M_PI) -
    std::log(a);

  beta = std::exp(log_beta);

  if (!R_finite(beta) || beta <= 0.0) {
    stop("Invalid log-logistic scale parameter.");
  }

  // ----------------------------------------------------------
  // Implicit derivatives of a
  //
  // Let
  //
  //   g(a) = log(1 + target_cv2).
  //
  // Then
  //
  //   da/dtheta
  //   =
  //   {
  //     d target_cv2/dtheta
  //     /
  //     (1 + target_cv2)
  //   }
  //   /
  //   g'(a).
  // ----------------------------------------------------------

  double gp =
    loglogistic_log_moment_ratio_deriv(a);

  if (!R_finite(gp) || gp <= 0.0) {
    stop("Invalid derivative while solving log-logistic parameters.");
  }

  double rhs_deriv_mu = 0.0;

  if (var_power != 2) {

    rhs_deriv_mu =
      static_cast<double>(var_power - 2) *
      target_cv2 /
      (
        mu *
        (1.0 + target_cv2)
      );
  }

  double rhs_deriv_phi =
    target_cv2 /
    (
      phi *
      (1.0 + target_cv2)
    );

  da_dmu = rhs_deriv_mu / gp;
  da_dphi = rhs_deriv_phi / gp;

  if (!R_finite(da_dmu) ||
      !R_finite(da_dphi)) {
    stop("Non-finite implicit derivatives for log-logistic parameters.");
  }
}


// ============================================================
// Optional exported parameter-conversion function
// ============================================================

// [[Rcpp::export]]
NumericMatrix loglogistic_mu_phi_parameters_cpp(
    NumericVector mu,
    double phi,
    int var_power = 2) {

  int n = mu.size();

  if (n == 0) {
    stop("mu must contain at least one value.");
  }

  if (!R_finite(phi) || phi <= 0.0) {
    stop("phi must be positive and finite.");
  }

  if (!(var_power == 1 ||
        var_power == 2 ||
        var_power == 3)) {
    stop("var_power must be one of 1, 2, or 3.");
  }

  NumericMatrix out(n, 8);

  colnames(out) = CharacterVector::create(
    "mu",
    "variance",
    "cv2",
    "a",
    "kappa",
    "beta",
    "log_beta",
    "mean_check"
  );

  for (int i = 0; i < n; ++i) {

    double mui = mu[i];

    double a;
    double kappa;
    double beta;
    double log_beta;
    double da_dmu;
    double da_dphi;

    loglogistic_mu_phi_parameters(
      mui,
      phi,
      var_power,
      a,
      kappa,
      beta,
      log_beta,
      da_dmu,
      da_dphi
    );

    double variance =
      phi *
      std::pow(
        mui,
        static_cast<double>(var_power)
      );

    double cv2 =
      phi *
      std::pow(
        mui,
        static_cast<double>(var_power - 2)
      );

    double log_A1 =
      loglogistic_log_A1(a);

    double mean_check =
      std::exp(log_beta + log_A1);

    out(i, 0) = mui;
    out(i, 1) = variance;
    out(i, 2) = cv2;
    out(i, 3) = a;
    out(i, 4) = kappa;
    out(i, 5) = beta;
    out(i, 6) = log_beta;
    out(i, 7) = mean_check;
  }

  return out;
}

// ============================================================
// I_loglogistic_mu_phi_cpp
//
// Computes
//
//   I_i^(m)(z)
//   = integral_{(m-1)/d}^{m/d}
//       exp{-Q_i(u)} Q_i(u)^z / Gamma(z+1) du.
//
// Parameterization:
//
//   E(Lambda_i)   = mu_i,
//   Var(Lambda_i) = phi * mu_i^var_power.
//
// Log-logistic quantile:
//
//   Q_i(u)
//   = beta_i * {u / (1-u)}^a_i,
//
// where
//
//   a_i     = 1 / kappa_i,
//   beta_i  = mu_i * sin(pi a_i) / (pi a_i).
//
// Each checkerboard cell is transformed using
//
//   u = lower + width * (1 - exp(-s)),
//   du = width * exp(-s) ds,
//
// and evaluated by Gauss-Laguerre quadrature.
//
// The supplied Laguerre nodes and weights approximate
//
//   integral_0^infinity exp(-s) g(s) ds.
// ============================================================

// [[Rcpp::export]]
NumericMatrix I_loglogistic_mu_phi_cpp(
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

  // ----------------------------------------------------------
  // Input checks
  // ----------------------------------------------------------

  if (n == 0) {
    stop("z and mu must contain at least one observation.");
  }

  if (mu.size() != n) {
    stop("length(mu) must equal length(z).");
  }

  if (!R_finite(phi) || phi <= 0.0) {
    stop("phi must be positive and finite.");
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
    stop(
      "lag_nodes and lag_weights must be supplied "
      "and have equal positive lengths."
    );
  }

  if (!R_finite(eps_u) ||
      eps_u <= 0.0 ||
      eps_u >= 0.5) {
    stop("eps_u must lie in (0, 0.5).");
  }

  for (int k = 0; k < K; ++k) {

    if (!R_finite(lag_nodes[k]) ||
        lag_nodes[k] < 0.0) {
      stop("lag_nodes must be finite and nonnegative.");
    }

    if (!R_finite(lag_weights[k]) ||
        lag_weights[k] <= 0.0) {
      stop("lag_weights must be positive and finite.");
    }
  }

  // ----------------------------------------------------------
  // Output
  // ----------------------------------------------------------

  NumericMatrix out(n, d);

  // ----------------------------------------------------------
  // Precompute cell-specific transformed quantities
  //
  // For each cell and Laguerre node, compute
  //
  //   logit(u) = log(u) - log(1-u).
  //
  // Under
  //
  //   u = lower + width * (1-exp(-s)),
  //
  // we have
  //
  //   1-u = (1-upper) + width exp(-s),
  //
  // for upper < 1, and
  //
  //   1-u = width exp(-s),
  //
  // for the final cell.
  // ----------------------------------------------------------

  std::vector< std::vector<double> > logit_u_grid(
    d,
    std::vector<double>(K)
  );

  std::vector<double> log_lag_weight(K);

  for (int k = 0; k < K; ++k) {
    log_lag_weight[k] =
      std::log(lag_weights[k]);
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

      double exp_minus_s = std::exp(-s);

      double one_minus_u;

      if (upper < 1.0) {

        one_minus_u =
          (1.0 - upper) +
          width * exp_minus_s;

      } else {

        one_minus_u =
          width * exp_minus_s;
      }

      one_minus_u =
        std::min(
          1.0 - eps_u,
          std::max(one_minus_u, eps_u)
        );

      double u = 1.0 - one_minus_u;

      u = clamp_u_loglogistic(u, eps_u);

      one_minus_u = 1.0 - u;

      double logit_u =
        std::log(u) -
        std::log(one_minus_u);

      logit_u_grid[m][k] = logit_u;
    }
  }

  // ----------------------------------------------------------
  // Main loop
  // ----------------------------------------------------------

  const double log_width =
    -std::log(static_cast<double>(d));

  for (int i = 0; i < n; ++i) {

    if (z[i] < 0) {
      stop("z must contain nonnegative integers.");
    }

    if (!R_finite(mu[i]) || mu[i] <= 0.0) {
      stop("mu must contain positive finite values.");
    }

    double a;
    double kappa;
    double beta;
    double log_beta;
    double da_dmu;
    double da_dphi;

    loglogistic_mu_phi_parameters(
      mu[i],
      phi,
      var_power,
      a,
      kappa,
      beta,
      log_beta,
      da_dmu,
      da_dphi
    );

    double zi =
      static_cast<double>(z[i]);

    double lgamma_z1 =
      R::lgammafn(zi + 1.0);

    for (int m = 0; m < d; ++m) {

      double max_log_term = R_NegInf;

      std::vector<double> log_terms(K);

      for (int k = 0; k < K; ++k) {

        // ----------------------------------------------------
        // Log-logistic quantile
        //
        // log Q
        // = log beta + a log{u/(1-u)}.
        // ----------------------------------------------------

        double logQ =
          log_beta +
          a * logit_u_grid[m][k];

        double log_integrand = R_NegInf;

        // exp(logQ) overflows near logQ = 709.
        // A slightly conservative threshold is used.
        if (R_finite(logQ) && logQ <= 700.0) {

          double Q = std::exp(logQ);

          double zlogQ =
            (z[i] == 0) ?
            0.0 :
            zi * logQ;

          log_integrand =
            zlogQ -
            Q -
            lgamma_z1;
        }

        // The cell width appears outside the Laguerre weight:
        //
        // I_m
        // = width *
        //   integral exp(-s) h{Q(u_m(s))} ds.
        double log_term =
          log_width +
          log_lag_weight[k] +
          log_integrand;

        log_terms[k] = log_term;

        if (log_term > max_log_term) {
          max_log_term = log_term;
        }
      }

      // ------------------------------------------------------
      // Stable log-sum-exp
      // ------------------------------------------------------

      double log_value = R_NegInf;

      if (R_finite(max_log_term)) {

        double scaled_sum = 0.0;

        for (int k = 0; k < K; ++k) {

          if (R_finite(log_terms[k])) {

            scaled_sum +=
              std::exp(
                log_terms[k] -
                max_log_term
              );
          }
        }

        if (scaled_sum > 0.0 &&
            R_finite(scaled_sum)) {

          log_value =
            max_log_term +
            std::log(scaled_sum);
        }
      }

      out(i, m) =
        log_out ?
        log_value :
        std::exp(log_value);
    }
  }

  return out;
}

// ============================================================
// I_loglogistic_mu_phi_grad_cpp
//
// Returns:
//
//   I       : checkerboard cell probabilities,
//   value   : alias of I,
//   d_mu    : derivative with respect to mu_i,
//   d_phi   : derivative with respect to phi.
//
// After solving for a_i = 1/kappa_i, the derivatives are
// analytic.
//
// Since
//
//   log beta
//   = log mu
//     + log sin(pi a)
//     - log(pi a),
//
// we have
//
//   d log beta / d theta
//   = 1(theta=mu)/mu
//     + da/dtheta *
//       {pi cot(pi a) - 1/a}.
//
// Furthermore,
//
//   log Q(u)
//   = log beta
//     + a log{u/(1-u)},
//
// and therefore
//
//   d log Q / d theta
//   = 1(theta=mu)/mu
//     + da/dtheta *
//       [
//         pi cot(pi a)
//         - 1/a
//         + log{u/(1-u)}
//       ].
//
// For the Poisson kernel
//
//   h(Q) = exp(-Q) Q^z / Gamma(z+1),
//
//   d h / d theta
//   = h(Q) (z-Q) d log Q/d theta.
// ============================================================

// [[Rcpp::export]]
List I_loglogistic_mu_phi_grad_cpp(
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

  // ----------------------------------------------------------
  // Input checks
  // ----------------------------------------------------------

  if (n == 0) {
    stop("z and mu must contain at least one observation.");
  }

  if (mu.size() != n) {
    stop("length(mu) must equal length(z).");
  }

  if (!R_finite(phi) || phi <= 0.0) {
    stop("phi must be positive and finite.");
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
    stop(
      "lag_nodes and lag_weights must be supplied "
      "and have equal positive lengths."
    );
  }

  if (!R_finite(eps_u) ||
      eps_u <= 0.0 ||
      eps_u >= 0.5) {
    stop("eps_u must lie in (0, 0.5).");
  }

  for (int k = 0; k < K; ++k) {

    if (!R_finite(lag_nodes[k]) ||
        lag_nodes[k] < 0.0) {
      stop("lag_nodes must be finite and nonnegative.");
    }

    if (!R_finite(lag_weights[k]) ||
        lag_weights[k] <= 0.0) {
      stop("lag_weights must be positive and finite.");
    }
  }

  // ----------------------------------------------------------
  // Outputs
  // ----------------------------------------------------------

  NumericMatrix prob(n, d);
  NumericMatrix d_mu_out(n, d);
  NumericMatrix d_phi_out(n, d);

  // ----------------------------------------------------------
  // Precompute transformed quadrature grid
  //
  // Store
  //
  //   logit(u) = log{u/(1-u)}
  //
  // for each cell and each Laguerre node.
  // ----------------------------------------------------------

  std::vector< std::vector<double> > logit_u_grid(
    d,
    std::vector<double>(K)
  );

  std::vector<double> log_lag_weight(K);

  for (int k = 0; k < K; ++k) {
    log_lag_weight[k] =
      std::log(lag_weights[k]);
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

      double exp_minus_s = std::exp(-s);

      double one_minus_u;

      if (upper < 1.0) {

        one_minus_u =
          (1.0 - upper) +
          width * exp_minus_s;

      } else {

        one_minus_u =
          width * exp_minus_s;
      }

      one_minus_u =
        std::min(
          1.0 - eps_u,
          std::max(one_minus_u, eps_u)
        );

      double u = 1.0 - one_minus_u;

      u = clamp_u_loglogistic(u, eps_u);

      one_minus_u = 1.0 - u;

      logit_u_grid[m][k] =
        std::log(u) -
        std::log(one_minus_u);
    }
  }

  // ----------------------------------------------------------
  // Main loop
  // ----------------------------------------------------------

  const double log_width =
    -std::log(static_cast<double>(d));

  for (int i = 0; i < n; ++i) {

    if (z[i] < 0) {
      stop("z must contain nonnegative integers.");
    }

    if (!R_finite(mu[i]) || mu[i] <= 0.0) {
      stop("mu must contain positive finite values.");
    }

    double mui = mu[i];

    double zi =
      static_cast<double>(z[i]);

    double lgamma_z1 =
      R::lgammafn(zi + 1.0);

    double a;
    double kappa;
    double beta;
    double log_beta;
    double da_dmu;
    double da_dphi;

    loglogistic_mu_phi_parameters(
      mui,
      phi,
      var_power,
      a,
      kappa,
      beta,
      log_beta,
      da_dmu,
      da_dphi
    );

    // --------------------------------------------------------
    // Derivative contribution from log beta:
    //
    // d log beta / da
    // = pi cot(pi a) - 1/a.
    // --------------------------------------------------------

    double pi_a = M_PI * a;

    double sin_pi_a = std::sin(pi_a);

    if (!R_finite(sin_pi_a) ||
        sin_pi_a <= 0.0) {
      stop(
        "Invalid trigonometric value while computing "
        "log-logistic gradients."
      );
    }

    double cot_pi_a =
      std::cos(pi_a) /
      sin_pi_a;

    double dlogbeta_da =
      M_PI * cot_pi_a -
      1.0 / a;

    for (int m = 0; m < d; ++m) {

      double max_log_term = R_NegInf;

      std::vector<double> log_terms(K);
      std::vector<double> score_mu(K);
      std::vector<double> score_phi(K);

      for (int k = 0; k < K; ++k) {

        double logit_u =
          logit_u_grid[m][k];

        // ----------------------------------------------------
        // Quantile
        // ----------------------------------------------------

        double logQ =
          log_beta +
          a * logit_u;

        double log_integrand = R_NegInf;
        double Q = R_PosInf;

        if (R_finite(logQ) && logQ <= 700.0) {

          Q = std::exp(logQ);

          double zlogQ =
            (z[i] == 0) ?
            0.0 :
            zi * logQ;

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

        // ----------------------------------------------------
        // Analytic derivatives of log Q
        //
        // common_a_term
        // = d log Q / da, holding mu fixed except through a.
        // ----------------------------------------------------

        double common_a_term =
          dlogbeta_da +
          logit_u;

        double dlogQ_dmu =
          1.0 / mui +
          da_dmu * common_a_term;

        double dlogQ_dphi =
          da_dphi * common_a_term;

        if (R_finite(Q)) {

          score_mu[k] =
            (zi - Q) *
            dlogQ_dmu;

          score_phi[k] =
            (zi - Q) *
            dlogQ_dphi;

        } else {

          // The corresponding Poisson kernel is numerically zero.
          score_mu[k] = 0.0;
          score_phi[k] = 0.0;
        }
      }

      // ------------------------------------------------------
      // Simultaneous stable computation of
      //
      //   I,
      //   dI/dmu,
      //   dI/dphi.
      // ------------------------------------------------------

      double value = 0.0;
      double derivative_mu = 0.0;
      double derivative_phi = 0.0;

      if (R_finite(max_log_term)) {

        double scaled_value = 0.0;
        double scaled_dmu = 0.0;
        double scaled_dphi = 0.0;

        for (int k = 0; k < K; ++k) {

          if (!R_finite(log_terms[k])) {
            continue;
          }

          double scaled_term =
            std::exp(
              log_terms[k] -
              max_log_term
            );

          scaled_value += scaled_term;

          scaled_dmu +=
            scaled_term *
            score_mu[k];

          scaled_dphi +=
            scaled_term *
            score_phi[k];
        }

        double scale_factor =
          std::exp(max_log_term);

        value =
          scale_factor *
          scaled_value;

        derivative_mu =
          scale_factor *
          scaled_dmu;

        derivative_phi =
          scale_factor *
          scaled_dphi;
      }

      prob(i, m) = value;

      d_mu_out(i, m) =
        derivative_mu;

      d_phi_out(i, m) =
        derivative_phi;
    }
  }

  return List::create(
    _["prob"] = prob,
    _["d_mu"] = d_mu_out,
    _["d_phi"] = d_phi_out
  );
}
