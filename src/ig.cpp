#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

inline double clamp_prob(double p) {
  if (p < 1e-14) return 1e-14;
  if (p > 1.0 - 1e-14) return 1.0 - 1e-14;
  return p;
}

inline double ig_lambda(double mu, double phi, int var_power) {
  // Var(X) = mu^3 / lambda = phi * mu^var_power
  // lambda = mu^(3 - var_power) / phi
  return std::pow(mu, 3.0 - var_power) / phi;
}

inline double ig_cdf_one(double x, double mu, double lambda) {
  if (x <= 0.0) return 0.0;
  if (!R_finite(x)) return 1.0;
  
  double sx = std::sqrt(lambda / x);
  double a1 = sx * (x / mu - 1.0);
  double a2 = -sx * (x / mu + 1.0);
  
  double term1 = R::pnorm(a1, 0.0, 1.0, 1, 0);
  double log_term2 = 2.0 * lambda / mu +
    R::pnorm(a2, 0.0, 1.0, 1, 1);
  
  double term2 = (log_term2 > 700.0) ? 0.0 : std::exp(log_term2);
  
  double out = term1 + term2;
  if (out < 0.0) out = 0.0;
  if (out > 1.0) out = 1.0;
  return out;
}

inline double ig_quantile_one(double p, double mu, double lambda) {
  p = clamp_prob(p);
  
  double low = 0.0;
  double high = mu;
  
  while (ig_cdf_one(high, mu, lambda) < p && high < 1e100) {
    high *= 2.0;
  }
  
  for (int it = 0; it < 100; ++it) {
    double mid = 0.5 * (low + high);
    double Fmid = ig_cdf_one(mid, mu, lambda);
    
    if (Fmid < p) low = mid;
    else high = mid;
  }
  
  return 0.5 * (low + high);
}

inline double logsumexp_vec(const std::vector<double>& x) {
  double m = R_NegInf;
  for (double v : x) if (v > m) m = v;
  if (!R_finite(m)) return m;
  
  double s = 0.0;
  for (double v : x) s += std::exp(v - m);
  
  return m + std::log(s);
}

struct LogBesselKHalfResult {
  double log_value;
  double dlog_dx;
};

inline LogBesselKHalfResult log_besselK_half_integer_with_deriv_one(
    double x,
    int z) {
  // K_{z-1/2}(x) = K_{|z-1/2|}(x). For z = 0, this is K_{1/2}(x).
  int m = (z == 0) ? 0 : z - 1;
  
  double log_poly = 0.0;
  double mean_k = 0.0;
  
  if (m > 0) {
    std::vector<double> terms(m + 1);
    
    for (int k = 0; k <= m; ++k) {
      terms[k] =
        R::lgammafn(m + k + 1.0) -
        R::lgammafn(k + 1.0) -
        R::lgammafn(m - k + 1.0) -
        static_cast<double>(k) * std::log(2.0 * x);
    }
    
    log_poly = logsumexp_vec(terms);
    
    // P'(x)/P(x) = -E_w(k)/x, where w_k is proportional
    // to the kth term in the half-integer Bessel polynomial.
    double max_term = *std::max_element(terms.begin(), terms.end());
    double denom = 0.0;
    double numer = 0.0;
    
    for (int k = 0; k <= m; ++k) {
      double wk = std::exp(terms[k] - max_term);
      denom += wk;
      numer += static_cast<double>(k) * wk;
    }
    
    mean_k = numer / denom;
  }
  
  LogBesselKHalfResult out;
  out.log_value =
    0.5 * (std::log(M_PI) - std::log(2.0 * x)) -
    x + log_poly;
  out.dlog_dx = -1.0 - (0.5 + mean_k) / x;
  
  return out;
}

inline double log_besselK_half_integer_one(double x, int z) {
  return log_besselK_half_integer_with_deriv_one(x, z).log_value;
}

struct IGClosedResult {
  double log_value;
  double dlog_dmu;
  double dlog_dphi;
};

inline IGClosedResult ig_I_d1_closed_log_grad_one(int z,
                                                   double mu,
                                                   double phi,
                                                   int var_power) {
  double lambda = ig_lambda(mu, phi, var_power);
  
  double nu = static_cast<double>(z) - 0.5;
  double beta = 1.0 + lambda / (2.0 * mu * mu);
  double gamma = lambda / 2.0;
  double x0 = 2.0 * std::sqrt(beta * gamma);
  
  LogBesselKHalfResult bessel =
    log_besselK_half_integer_with_deriv_one(x0, z);
  
  double logI =
    0.5 * (std::log(lambda) - std::log(2.0 * M_PI)) -
    R::lgammafn(static_cast<double>(z) + 1.0) +
    lambda / mu +
    std::log(2.0) +
    0.5 * nu * (std::log(gamma) - std::log(beta)) +
    bessel.log_value;
  
  // Generic derivatives for Var(Lambda) = phi * mu^var_power,
  // lambda = mu^(3-var_power) / phi.
  double dloglambda_dmu = (3.0 - var_power) / mu;
  double dloglambda_dphi = -1.0 / phi;
  
  double t = lambda / (2.0 * mu * mu); // beta - 1
  
  double dbeta_dmu =
    t * (dloglambda_dmu - 2.0 / mu);
  double dbeta_dphi =
    t * dloglambda_dphi;
  
  double dlogbeta_dmu = dbeta_dmu / beta;
  double dlogbeta_dphi = dbeta_dphi / beta;
  
  double dloggamma_dmu = dloglambda_dmu;
  double dloggamma_dphi = dloglambda_dphi;
  
  double dlambda_over_mu_dmu =
    (lambda / mu) * (dloglambda_dmu - 1.0 / mu);
  double dlambda_over_mu_dphi =
    (lambda / mu) * dloglambda_dphi;
  
  double dx0_dmu =
    0.5 * x0 * (dlogbeta_dmu + dloggamma_dmu);
  double dx0_dphi =
    0.5 * x0 * (dlogbeta_dphi + dloggamma_dphi);
  
  double dlogI_dmu =
    0.5 * dloglambda_dmu +
    dlambda_over_mu_dmu +
    0.5 * nu * (dloggamma_dmu - dlogbeta_dmu) +
    bessel.dlog_dx * dx0_dmu;
  
  double dlogI_dphi =
    0.5 * dloglambda_dphi +
    dlambda_over_mu_dphi +
    0.5 * nu * (dloggamma_dphi - dlogbeta_dphi) +
    bessel.dlog_dx * dx0_dphi;
  
  // The exact probability is at most one. Clipping is only a safeguard
  // against tiny positive round-off errors in logI.
  if (logI > 0.0) {
    logI = 0.0;
    dlogI_dmu = 0.0;
    dlogI_dphi = 0.0;
  }
  
  IGClosedResult out;
  out.log_value = logI;
  out.dlog_dmu = dlogI_dmu;
  out.dlog_dphi = dlogI_dphi;
  return out;
}

inline double ig_I_d1_closed_log_one(int z,
                                     double mu,
                                     double phi,
                                     int var_power) {
  return ig_I_d1_closed_log_grad_one(z, mu, phi, var_power).log_value;
}

inline double ig_log_density_one(double x,
                                 double mu,
                                 double phi,
                                 int var_power,
                                 double eps_x = 1e-300) {
  if (x < eps_x) x = eps_x;
  
  double lambda = ig_lambda(mu, phi, var_power);
  
  return 0.5 * (std::log(lambda) - std::log(2.0 * M_PI)) -
    1.5 * std::log(x) -
    lambda * (x - mu) * (x - mu) / (2.0 * mu * mu * x);
}

inline double ig_log_h_one(double x,
                           int z,
                           double mu,
                           double phi,
                           int var_power,
                           double eps_x = 1e-300) {
  if (x < eps_x) x = eps_x;
  
  double zlogx = (z == 0) ? 0.0 : static_cast<double>(z) * std::log(x);
  
  return zlogx -
    x -
    R::lgammafn(static_cast<double>(z) + 1.0) +
    ig_log_density_one(x, mu, phi, var_power, eps_x);
}

inline double ig_score_mu_logf_one(double x,
                                   double mu,
                                   double phi,
                                   int var_power,
                                   double eps_x = 1e-300) {
  if (x < eps_x) x = eps_x;
  
  double lambda = ig_lambda(mu, phi, var_power);
  double dlambda_dmu = (3.0 - var_power) * lambda / mu;
  
  double D = (x - mu) * (x - mu) / (2.0 * mu * mu * x);
  double dD_dmu = -x / (mu * mu * mu) + 1.0 / (mu * mu);
  
  return 0.5 * dlambda_dmu / lambda -
    dlambda_dmu * D -
    lambda * dD_dmu;
}

inline double ig_score_phi_logf_one(double x,
                                    double mu,
                                    double phi,
                                    int var_power,
                                    double eps_x = 1e-300) {
  if (x < eps_x) x = eps_x;
  
  double lambda = ig_lambda(mu, phi, var_power);
  double dlambda_dphi = -lambda / phi;
  
  double D = (x - mu) * (x - mu) / (2.0 * mu * mu * x);
  
  return 0.5 * dlambda_dphi / lambda -
    dlambda_dphi * D;
}

struct IGCDFDerivResult {
  double d_mu;
  double d_phi;
};

inline double safe_exp_one(double logx) {
  if (logx <= std::log(std::numeric_limits<double>::min())) return 0.0;
  if (logx >= std::log(std::numeric_limits<double>::max())) {
    return std::numeric_limits<double>::infinity();
  }
  return std::exp(logx);
}

inline IGCDFDerivResult ig_cdf_deriv_one(double x,
                                        double mu,
                                        double phi,
                                        int var_power) {
  IGCDFDerivResult out;
  out.d_mu = 0.0;
  out.d_phi = 0.0;

  if (x <= 0.0 || !R_finite(x)) return out;

  double lambda = ig_lambda(mu, phi, var_power);
  double s = std::sqrt(lambda / x);
  double r = x / mu;

  double a1 = s * (r - 1.0);
  double a2 = -s * (r + 1.0);
  double c = 2.0 * lambda / mu;

  double dloglambda_dmu = (3.0 - var_power) / mu;
  double dloglambda_dphi = -1.0 / phi;

  double ds_dmu = 0.5 * s * dloglambda_dmu;
  double ds_dphi = 0.5 * s * dloglambda_dphi;

  double dr_dmu = -r / mu;

  double da1_dmu = ds_dmu * (r - 1.0) + s * dr_dmu;
  double da2_dmu = -ds_dmu * (r + 1.0) - s * dr_dmu;

  double da1_dphi = ds_dphi * (r - 1.0);
  double da2_dphi = -ds_dphi * (r + 1.0);

  double dc_dmu = c * (dloglambda_dmu - 1.0 / mu);
  double dc_dphi = c * dloglambda_dphi;

  double normal_density_a1 = R::dnorm(a1, 0.0, 1.0, 0);

  // exp(c) Phi(a2), evaluated on the log scale.
  double log_Phi_a2 = R::pnorm(a2, 0.0, 1.0, 1, 1);
  double log_term2 = c + log_Phi_a2;
  double term2 = safe_exp_one(log_term2);

  // phi(a2) / Phi(a2), evaluated as a log inverse Mills ratio.
  double log_normal_density_a2 = R::dnorm(a2, 0.0, 1.0, 1);
  double mills_a2 = safe_exp_one(log_normal_density_a2 - log_Phi_a2);

  out.d_mu =
    normal_density_a1 * da1_dmu +
    term2 * (dc_dmu + mills_a2 * da2_dmu);

  out.d_phi =
    normal_density_a1 * da1_dphi +
    term2 * (dc_dphi + mills_a2 * da2_dphi);

  return out;
}

struct IGQuantileDerivResult {
  double d_mu;
  double d_phi;
};

inline IGQuantileDerivResult ig_quantile_deriv_one(double q,
                                                   double mu,
                                                   double phi,
                                                   int var_power,
                                                   double eps_density = 1e-300) {
  IGQuantileDerivResult out;
  out.d_mu = 0.0;
  out.d_phi = 0.0;

  if (q <= 0.0 || !R_finite(q)) return out;

  IGCDFDerivResult dF =
    ig_cdf_deriv_one(q, mu, phi, var_power);

  double log_density =
    ig_log_density_one(q, mu, phi, var_power, eps_density);
  double density =
    std::exp(std::max(log_density, std::log(eps_density)));

  // F(Q(p; theta); theta) = p implies
  // dQ/dtheta = -F_theta(Q; theta) / f(Q; theta).
  out.d_mu = -dF.d_mu / density;
  out.d_phi = -dF.d_phi / density;

  return out;
}


inline void ig_finite_interval_one(int z,
                                   double mu,
                                   double phi,
                                   int var_power,
                                   double L,
                                   double U,
                                   double dL_dmu,
                                   double dU_dmu,
                                   double dL_dphi,
                                   double dU_dphi,
                                   NumericVector leg_nodes,
                                   NumericVector leg_weights,
                                   double& value,
                                   double& d_mu,
                                   double& d_phi) {
  int K = leg_nodes.size();
  
  double half = 0.5 * (U - L);
  double mid  = 0.5 * (U + L);
  
  double val = 0.0;
  double dm = 0.0;
  double dp = 0.0;
  
  for (int k = 0; k < K; ++k) {
    double x = mid + half * leg_nodes[k];
    double w = half * leg_weights[k];
    
    double logh = ig_log_h_one(x, z, mu, phi, var_power);
    double hval = std::exp(logh);
    
    val += w * hval;
    dm  += w * hval * ig_score_mu_logf_one(x, mu, phi, var_power);
    dp  += w * hval * ig_score_phi_logf_one(x, mu, phi, var_power);
  }
  
  double hU = std::exp(ig_log_h_one(U, z, mu, phi, var_power));
  double hL = (L <= 0.0) ? 0.0 :
    std::exp(ig_log_h_one(L, z, mu, phi, var_power));
  
  value = val;
  d_mu  = dm + hU * dU_dmu - hL * dL_dmu;
  d_phi = dp + hU * dU_dphi - hL * dL_dphi;
}

// ============================================================
// I_ig_mu_phi_cpp
// ============================================================

// [[Rcpp::export]]
NumericMatrix I_ig_mu_phi_cpp(IntegerVector z,
                              NumericVector mu,
                              double phi,
                              int d = 1,
                              int var_power = 3,
                              NumericVector leg_nodes = NumericVector::create(),
                              NumericVector leg_weights = NumericVector::create(),
                              bool log_out = false,
                              double rel_step_mu = 1e-6,
                              double rel_step_phi = 1e-6,
                              double eps_tail = 1e-300) {
  
  int n = z.size();
  if (mu.size() != n) stop("length(mu) must equal length(z).");
  if (!R_finite(phi) || phi <= 0.0) stop("phi must be positive.");
  if (d <= 0) stop("d must be positive.");
  if (!(var_power == 1 || var_power == 2 || var_power == 3)) {
    stop("var_power must be one of 1, 2, or 3.");
  }
  if (d > 1 && leg_nodes.size() == 0) {
    stop("leg_nodes and leg_weights must be supplied when d > 1.");
  }
  
  NumericMatrix out(n, d);
  
  for (int i = 0; i < n; ++i) {
    if (z[i] < 0) stop("z must be nonnegative.");
    if (!R_finite(mu[i]) || mu[i] <= 0.0) stop("mu must be positive.");
    
    double mui = mu[i];
    double I_closed = std::exp(ig_I_d1_closed_log_one(z[i], mui, phi, var_power));
    
    if (d == 1) {
      out(i, 0) = log_out ? std::log(std::max(I_closed, eps_tail)) : I_closed;
      continue;
    }
    
    double sum_prev = 0.0;
    
    for (int m = 0; m < d - 1; ++m) {
      double pL = static_cast<double>(m) / static_cast<double>(d);
      double pU = static_cast<double>(m + 1) / static_cast<double>(d);
      
      double lambda = ig_lambda(mui, phi, var_power);
      
      double L = (m == 0) ? 0.0 : ig_quantile_one(pL, mui, lambda);
      double U = ig_quantile_one(pU, mui, lambda);
      
      double tmp, dm, dp;
      ig_finite_interval_one(
        z[i], mui, phi, var_power,
        L, U,
        0.0, 0.0, 0.0, 0.0,
        leg_nodes, leg_weights,
        tmp, dm, dp
      );
      
      out(i, m) = tmp;
      sum_prev += tmp;
    }
    
    double last = I_closed - sum_prev;
    if (last < eps_tail) last = eps_tail;
    out(i, d - 1) = last;
    
    if (log_out) {
      for (int m = 0; m < d; ++m) {
        out(i, m) = std::log(std::max(out(i, m), eps_tail));
      }
    }
  }
  
  return out;
}

// ============================================================
// I_ig_mu_phi_grad_cpp
// ============================================================

// [[Rcpp::export]]
List I_ig_mu_phi_grad_cpp(IntegerVector z,
                          NumericVector mu,
                          double phi,
                          int d = 1,
                          int var_power = 3,
                          NumericVector leg_nodes = NumericVector::create(),
                          NumericVector leg_weights = NumericVector::create(),
                          double rel_step_mu = 1e-6,
                          double rel_step_phi = 1e-6,
                          double eps_tail = 1e-300) {
  
  // Retained for backward compatibility with the R interface.
  // Quantile derivatives are now computed analytically.
  (void) rel_step_mu;
  (void) rel_step_phi;

  int n = z.size();
  if (mu.size() != n) stop("length(mu) must equal length(z).");
  if (!R_finite(phi) || phi <= 0.0) stop("phi must be positive.");
  if (d <= 0) stop("d must be positive.");
  if (!(var_power == 1 || var_power == 2 || var_power == 3)) {
    stop("var_power must be one of 1, 2, or 3.");
  }
  if (d > 1 && leg_nodes.size() == 0) {
    stop("leg_nodes and leg_weights must be supplied when d > 1.");
  }
  
  NumericMatrix prob(n, d);
  NumericMatrix d_mu(n, d);
  NumericMatrix d_phi(n, d);
  
  for (int i = 0; i < n; ++i) {
    if (z[i] < 0) stop("z must be nonnegative.");
    if (!R_finite(mu[i]) || mu[i] <= 0.0) stop("mu must be positive.");
    
    double mui = mu[i];
    
    IGClosedResult closed =
      ig_I_d1_closed_log_grad_one(z[i], mui, phi, var_power);
    
    double I_closed =
      std::exp(std::max(closed.log_value, std::log(eps_tail)));
    double d_mu_closed = I_closed * closed.dlog_dmu;
    double d_phi_closed = I_closed * closed.dlog_dphi;
    
    if (d == 1) {
      prob(i, 0) = I_closed;
      d_mu(i, 0) = d_mu_closed;
      d_phi(i, 0) = d_phi_closed;
      continue;
    }
    
    double sum_val = 0.0;
    double sum_dm = 0.0;
    double sum_dp = 0.0;
    
    for (int m = 0; m < d - 1; ++m) {
      double pL = static_cast<double>(m) / static_cast<double>(d);
      double pU = static_cast<double>(m + 1) / static_cast<double>(d);
      
      double lambda = ig_lambda(mui, phi, var_power);
      
      double L = (m == 0) ? 0.0 : ig_quantile_one(pL, mui, lambda);
      double U = ig_quantile_one(pU, mui, lambda);
      
      IGQuantileDerivResult dL;
      dL.d_mu = 0.0;
      dL.d_phi = 0.0;

      if (m > 0) {
        dL = ig_quantile_deriv_one(
          L, mui, phi, var_power, eps_tail
        );
      }

      IGQuantileDerivResult dU =
        ig_quantile_deriv_one(
          U, mui, phi, var_power, eps_tail
        );

      double dL_dmu = dL.d_mu;
      double dU_dmu = dU.d_mu;
      double dL_dphi = dL.d_phi;
      double dU_dphi = dU.d_phi;

      double val_m, dm_m, dp_m;
      
      ig_finite_interval_one(
        z[i], mui, phi, var_power,
        L, U,
        dL_dmu, dU_dmu,
        dL_dphi, dU_dphi,
        leg_nodes, leg_weights,
        val_m, dm_m, dp_m
      );
      
      prob(i, m) = val_m;
      d_mu(i, m) = dm_m;
      d_phi(i, m) = dp_m;
      
      sum_val += val_m;
      sum_dm += dm_m;
      sum_dp += dp_m;
    }
    
    prob(i, d - 1) = std::max(I_closed - sum_val, eps_tail);
    d_mu(i, d - 1) = d_mu_closed - sum_dm;
    d_phi(i, d - 1) = d_phi_closed - sum_dp;
  }
  
  return List::create(
    _["prob"] = prob,
    _["d_mu"] = d_mu,
    _["d_phi"] = d_phi
  );
}
