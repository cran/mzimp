#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

inline double clamp_u(double u, double eps_u) {
  if (u < eps_u) return eps_u;
  if (u > 1.0 - eps_u) return 1.0 - eps_u;
  return u;
}

inline double safe_exp_700(double x) {
  if (x > 700.0) return R_PosInf;
  return std::exp(x);
}

inline double log_integrand_lnorm(int z, double logQ) {
  if (logQ > 700.0) return R_NegInf;
  
  double Q = std::exp(logQ);
  double zlogQ = (z == 0) ? 0.0 : static_cast<double>(z) * logQ;
  
  double out = zlogQ - Q - R::lgammafn(static_cast<double>(z) + 1.0);
  if (!R_finite(out)) return R_NegInf;
  return out;
}

inline double row_logsumexp_weighted(const std::vector<double>& logf,
                                     NumericVector w) {
  int K = logf.size();
  double maxv = R_NegInf;
  
  for (int k = 0; k < K; ++k) {
    double val = logf[k] + std::log(w[k]);
    if (val > maxv) maxv = val;
  }
  
  if (!R_finite(maxv)) return R_NegInf;
  
  double s = 0.0;
  for (int k = 0; k < K; ++k) {
    s += std::exp(logf[k] + std::log(w[k]) - maxv);
  }
  
  return maxv + std::log(s);
}


// ============================================================
// I_lnorm_mu_phi_cpp
// Optimized version
// ============================================================

// [[Rcpp::export]]
NumericMatrix I_lnorm_mu_phi_cpp(IntegerVector z,
                                 NumericVector mu,
                                 double phi,
                                 int d = 1,
                                 int var_power = 2,
                                 NumericVector gh_nodes = NumericVector::create(),
                                 NumericVector gh_weights = NumericVector::create(),
                                 bool log_out = false,
                                 double eps_u = 1e-16) {
  
  int n = z.size();
  int K = gh_nodes.size();
  
  if (mu.size() != n) stop("length(mu) must equal length(z).");
  if (K == 0 || gh_weights.size() != K) {
    stop("gh_nodes and gh_weights must be supplied and have the same length.");
  }
  if (!R_finite(phi) || phi <= 0.0) stop("phi must be positive.");
  if (d <= 0) stop("d must be a positive integer.");
  if (!(var_power == 1 || var_power == 2 || var_power == 3)) {
    stop("var_power must be one of 1, 2, or 3.");
  }
  
  NumericMatrix out(n, d);
  
  // ------------------------------------------------------------
  // Precompute qnorm grid and log weights
  // ------------------------------------------------------------
  
  std::vector<double> base_u(K);
  std::vector<double> log_w(K);
  
  for (int k = 0; k < K; ++k) {
    base_u[k] = clamp_u(
      R::pnorm(std::sqrt(2.0) * gh_nodes[k], 0.0, 1.0, 1, 0),
      eps_u
    );
    log_w[k] = std::log(gh_weights[k]);
  }
  
  std::vector< std::vector<double> > qn_grid(
    d, std::vector<double>(K)
  );
  
  for (int m = 0; m < d; ++m) {
    double a = static_cast<double>(m) / static_cast<double>(d);
    double b = static_cast<double>(m + 1) / static_cast<double>(d);
    double width = b - a;
    
    for (int k = 0; k < K; ++k) {
      double u = clamp_u(a + width * base_u[k], eps_u);
      qn_grid[m][k] = R::qnorm(u, 0.0, 1.0, 1, 0);
    }
  }
  
  // ------------------------------------------------------------
  // Main loop
  // ------------------------------------------------------------
  
  for (int i = 0; i < n; ++i) {
    
    if (z[i] < 0) stop("z must contain nonnegative integers.");
    if (!R_finite(mu[i]) || mu[i] <= 0.0) {
      stop("mu must contain positive finite values.");
    }
    
    double mui = mu[i];
    double logmu = std::log(mui);
    
    double B = phi * std::pow(mui, static_cast<double>(var_power - 2));
    double sigma2 = std::log1p(B);
    double sigma = std::sqrt(sigma2);
    
    double lgam_z1 = R::lgammafn(static_cast<double>(z[i]) + 1.0);
    
    for (int m = 0; m < d; ++m) {
      
      double width = 1.0 / static_cast<double>(d);
      double log_width = std::log(width);
      
      double maxv = R_NegInf;
      std::vector<double> logf(K);
      
      for (int k = 0; k < K; ++k) {
        
        double qn = qn_grid[m][k];
        double logQ = logmu - 0.5 * sigma2 + sigma * qn;
        
        double val;
        if (logQ > 700.0) {
          val = R_NegInf;
        } else {
          double Q = std::exp(logQ);
          double zlogQ = (z[i] == 0) ? 0.0 : static_cast<double>(z[i]) * logQ;
          val = zlogQ - Q - lgam_z1;
          if (!R_finite(val)) val = R_NegInf;
        }
        
        logf[k] = val + log_w[k];
        if (logf[k] > maxv) maxv = logf[k];
      }
      
      double log_sum = R_NegInf;
      
      if (R_finite(maxv)) {
        double s = 0.0;
        for (int k = 0; k < K; ++k) {
          s += std::exp(logf[k] - maxv);
        }
        log_sum = maxv + std::log(s);
      }
      
      double log_val =
        log_width -
        0.5 * std::log(M_PI) +
        log_sum;
      
      out(i, m) = log_out ? log_val : std::exp(log_val);
    }
  }
  
  return out;
}

// ============================================================
// I_lnorm_mu_phi_grad_cpp
// Optimized version
// ============================================================

// [[Rcpp::export]]
List I_lnorm_mu_phi_grad_cpp(IntegerVector z,
                             NumericVector mu,
                             double phi,
                             int d = 1,
                             int var_power = 2,
                             NumericVector gh_nodes = NumericVector::create(),
                             NumericVector gh_weights = NumericVector::create(),
                             double eps_u = 1e-16) {
  
  int n = z.size();
  int K = gh_nodes.size();
  
  if (mu.size() != n) stop("length(mu) must equal length(z).");
  if (K == 0 || gh_weights.size() != K) {
    stop("gh_nodes and gh_weights must be supplied and have the same length.");
  }
  if (!R_finite(phi) || phi <= 0.0) stop("phi must be positive.");
  if (d <= 0) stop("d must be a positive integer.");
  if (!(var_power == 1 || var_power == 2 || var_power == 3)) {
    stop("var_power must be one of 1, 2, or 3.");
  }
  
  NumericMatrix prob(n, d);
  NumericMatrix d_mu_out(n, d);
  NumericMatrix d_phi_out(n, d);
  
  // ------------------------------------------------------------
  // Precompute qnorm grid and effective weights
  // ------------------------------------------------------------
  
  std::vector<double> base_u(K);
  
  for (int k = 0; k < K; ++k) {
    base_u[k] = clamp_u(
      R::pnorm(std::sqrt(2.0) * gh_nodes[k], 0.0, 1.0, 1, 0),
      eps_u
    );
  }
  
  std::vector< std::vector<double> > qn_grid(
    d, std::vector<double>(K)
  );
  
  std::vector< std::vector<double> > w_eff_grid(
    d, std::vector<double>(K)
  );
  
  for (int m = 0; m < d; ++m) {
    double a = static_cast<double>(m) / static_cast<double>(d);
    double b = static_cast<double>(m + 1) / static_cast<double>(d);
    double width = b - a;
    
    for (int k = 0; k < K; ++k) {
      double u = clamp_u(a + width * base_u[k], eps_u);
      qn_grid[m][k] = R::qnorm(u, 0.0, 1.0, 1, 0);
      w_eff_grid[m][k] = gh_weights[k] * width / std::sqrt(M_PI);
    }
  }
  
  // ------------------------------------------------------------
  // Main loop
  // ------------------------------------------------------------
  
  for (int i = 0; i < n; ++i) {
    
    if (z[i] < 0) stop("z must contain nonnegative integers.");
    if (!R_finite(mu[i]) || mu[i] <= 0.0) {
      stop("mu must contain positive finite values.");
    }
    
    double mui = mu[i];
    double logmu = std::log(mui);
    double zi = static_cast<double>(z[i]);
    double lgam_z1 = R::lgammafn(zi + 1.0);
    
    double B = phi * std::pow(mui, static_cast<double>(var_power - 2));
    double sigma2 = std::log1p(B);
    double sigma = std::sqrt(sigma2);
    
    double d_sigma2_dmu = 0.0;
    if (var_power != 2) {
      d_sigma2_dmu =
        phi * static_cast<double>(var_power - 2) *
        std::pow(mui, static_cast<double>(var_power - 3)) /
        (1.0 + B);
    }
    
    double d_sigma2_dphi =
      std::pow(mui, static_cast<double>(var_power - 2)) /
      (1.0 + B);
    
    double d_sigma_dmu = 0.0;
    double d_sigma_dphi = 0.0;
    
    if (sigma > 0.0) {
      d_sigma_dmu  = 0.5 * d_sigma2_dmu  / sigma;
      d_sigma_dphi = 0.5 * d_sigma2_dphi / sigma;
    }
    
    for (int m = 0; m < d; ++m) {
      
      double I_sum = 0.0;
      double dmu_sum = 0.0;
      double dphi_sum = 0.0;
      
      for (int k = 0; k < K; ++k) {
        
        double qn = qn_grid[m][k];
        double logQ = logmu - 0.5 * sigma2 + sigma * qn;
        
        double f = 0.0;
        double Q = 0.0;
        
        if (logQ <= 700.0) {
          Q = std::exp(logQ);
          
          double zlogQ = (z[i] == 0) ? 0.0 : zi * logQ;
          double logf = zlogQ - Q - lgam_z1;
          
          if (R_finite(logf)) {
            f = std::exp(logf);
          }
        }
        
        double A = zi - Q;
        
        double dlogQ_dmu =
          1.0 / mui -
          0.5 * d_sigma2_dmu +
          d_sigma_dmu * qn;
        
        double dlogQ_dphi =
          -0.5 * d_sigma2_dphi +
          d_sigma_dphi * qn;
        
        double w_eff = w_eff_grid[m][k];
        
        I_sum    += w_eff * f;
        dmu_sum  += w_eff * f * A * dlogQ_dmu;
        dphi_sum += w_eff * f * A * dlogQ_dphi;
      }
      
      prob(i, m) = I_sum;
      d_mu_out(i, m) = dmu_sum;
      d_phi_out(i, m) = dphi_sum;
    }
  }
  
  return List::create(
    _["prob"] = prob,
    _["d_mu"] = d_mu_out,
    _["d_phi"] = d_phi_out
  );
}
