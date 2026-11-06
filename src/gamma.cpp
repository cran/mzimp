#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]


// ============================================================
// log(exp(a) - exp(b)), assuming a >= b
// ============================================================

inline double logspace_sub_cpp(double log_a, double log_b) {
  if (log_b == R_NegInf) return log_a;
  if (log_b >= log_a) return R_NegInf;
  return log_a + std::log1p(-std::exp(log_b - log_a));
}





// ------------------------------------------------------------
// stable log( P(yL < Gamma(s,1) <= yU) )
// ------------------------------------------------------------

inline double stable_pdiff_cpp(double yL,
                               double yU,
                               double shape)
{
    if (yL <= 0.0 && R_finite(yU))
        return R::pgamma(yU, shape, 1.0, 1, 1);

    if (!R_finite(yU))
        return (yL <= 0.0)
            ? 0.0
            : R::pgamma(yL, shape, 1.0, 0, 1);

    double logFU = R::pgamma(yU, shape, 1.0, 1, 1);
    double logFL = R::pgamma(yL, shape, 1.0, 1, 1);

    if (logFL >= logFU)
        return R_NegInf;

    return logFU + std::log1p(-std::exp(logFL - logFU));
}




// ============================================================
// lower_gamma_scaled_cpp
//
// Computes:
//   gamma(s, x) / Gamma(z + 1)
//
// where gamma(s,x) is the lower incomplete gamma.
// In the R code this corresponds to:
//   lower_gamma_scaled(s, x, lgam_z1)
// ============================================================

inline double lower_gamma_scaled_cpp(double s,
                                     double x,
                                     double lgam_z1) {
  
  if (x <= 0.0) {
    return 0.0;
  }
  
  if (!R_finite(x)) {
    return std::exp(R::lgammafn(s) - lgam_z1);
  }
  
  double logP = R::pgamma(x, s, 1.0, 1, 1);
  
  if (!R_finite(logP)) {
    return 0.0;
  }
  
  return std::exp(R::lgammafn(s) + logP - lgam_z1);
}


// ============================================================
// lower_gamma_scaled_dshape_cpp
//
// Central finite difference derivative with respect to shape s:
//   d/ds { gamma(s, x) / Gamma(z + 1) }
// ============================================================

inline double lower_gamma_scaled_dshape_cpp(double s,
                                            double x,
                                            double lgam_z1,
                                            double eps_shape) {
  
  double up = lower_gamma_scaled_cpp(
    s + eps_shape,
    x,
    lgam_z1
  );
  
  double dn = lower_gamma_scaled_cpp(
    s - eps_shape,
    x,
    lgam_z1
  );
  
  return (up - dn) / (2.0 * eps_shape);
}



// ============================================================
// h_scaled_cpp
//
// Computes:
//   exp(-x) * x^z / Gamma(z + 1)
//
// In the R code this corresponds to:
//   h_scaled(x, z, lgam_z1)
// ============================================================

inline double h_scaled_cpp(double x,
                           int z,
                           double lgam_z1) {
  
  if (!R_finite(x)) {
    return 0.0;
  }
  
  if (x < 0.0) {
    return R_NaN;
  }
  
  if (x == 0.0) {
    if (z == 0) {
      return std::exp(-lgam_z1);
    } else {
      return 0.0;
    }
  }
  
  return std::exp(
    -x + static_cast<double>(z) * std::log(x) - lgam_z1
  );
}


// ============================================================
// raw lower incomplete gamma:
//   gamma(shape, x)
// used in d_shape boundary terms
// ============================================================

inline double lower_gamma_raw_cpp(double shape,
                                  double x) {
  
  if (x <= 0.0) {
    return 0.0;
  }
  
  if (!R_finite(x)) {
    return R::gammafn(shape);
  }
  
  double logP = R::pgamma(x, shape, 1.0, 1, 1);
  
  if (!R_finite(logP)) {
    return 0.0;
  }
  
  return std::exp(R::lgammafn(shape) + logP);
}


// ============================================================
// derivative of raw lower incomplete gamma wrt shape
// by central finite difference
// ============================================================

inline double lower_gamma_raw_dshape_cpp(double shape,
                                         double x,
                                         double eps_shape) {
  
  if (!R_finite(x)) {
    return R::gammafn(shape) * R::digamma(shape);
  }
  
  if (x <= 0.0) {
    return 0.0;
  }
  
  double up = lower_gamma_raw_cpp(shape + eps_shape, x);
  double dn = lower_gamma_raw_cpp(shape - eps_shape, x);
  
  return (up - dn) / (2.0 * eps_shape);
}



// ============================================================
// I_gamma_mu_phi_cpp
// E(Lambda) = mu
// Var(Lambda) = phi * mu^var_power
// ============================================================

// [[Rcpp::export]]
NumericMatrix I_gamma_mu_phi_cpp(IntegerVector z,
                                  NumericVector mu,
                                  double phi,
                                  int d,
                                  int var_power = 2,
                                  bool log_out = false) {
  
  int n = z.size();
  
  if (mu.size() != n) {
    stop("length(mu) must equal length(z).");
  }
  
  if (!R_finite(phi) || phi <= 0.0) {
    stop("phi must be a positive scalar.");
  }
  
  if (d <= 0) {
    stop("d must be a positive integer.");
  }
  
  if (!(var_power == 1 || var_power == 2 || var_power == 3)) {
    stop("var_power must be one of 1, 2, or 3.");
  }
  
  NumericMatrix out(n, d);
  
  for (int i = 0; i < n; ++i) {
    
    if (z[i] < 0) {
      stop("z must contain nonnegative integers.");
    }
    
    if (!R_finite(mu[i]) || mu[i] <= 0.0) {
      stop("mu must contain positive finite values.");
    }
    
    double mui = mu[i];
    double zi = static_cast<double>(z[i]);
    
    double shape = std::pow(mui, 2.0 - var_power) / phi;
    double rate  = std::pow(mui, 1.0 - var_power) / phi;
    
    if (!R_finite(shape) || shape <= 0.0 ||
        !R_finite(rate)  || rate  <= 0.0) {
      stop("Invalid Gamma shape/rate implied by mu, phi, and var_power.");
    }
    
    // --------------------------------------------------------
    // d = 1: Gamma mixed Poisson = Negative Binomial
    // --------------------------------------------------------
    
    if (d == 1) {
      double logp = R::dnbinom_mu(zi, shape, mui, 1);
      out(i, 0) = log_out ? logp : std::exp(logp);
      continue;
    }
    
    // --------------------------------------------------------
    // d >= 2
    // --------------------------------------------------------
    
    double s = zi + shape;
    
    double log_const =
      shape * std::log(rate)
      - R::lgammafn(shape)
      + R::lgammafn(s)
      - R::lgammafn(zi + 1.0)
      - s * std::log(rate + 1.0);
    
    std::vector<double> x(d + 1);
    std::vector<double> y(d + 1);
    
    x[0] = 0.0;
    y[0] = 0.0;
    
    for (int m = 1; m < d; ++m) {
      double p = static_cast<double>(m) / static_cast<double>(d);
      x[m] = R::qgamma(p, shape, 1.0 / rate, 1, 0);
      y[m] = x[m] * (rate + 1.0);
    }
    
    x[d] = R_PosInf;
    y[d] = R_PosInf;
    
    for (int m = 0; m < d; ++m) {
      
      double log_prob_diff = stable_pdiff_cpp(y[m], y[m + 1], s);
      double log_val = log_const + log_prob_diff;
      
      out(i, m) = log_out ? log_val : std::exp(log_val);
    }
  }
  
  return out;
}




// ============================================================
// I_gamma_mu_phi_grad_cpp
// Gradient wrt mu and phi
// E(Lambda)   = mu
// Var(Lambda) = phi * mu^var_power
// ============================================================

// [[Rcpp::export]]
List I_gamma_mu_phi_grad_cpp(IntegerVector z,
                             NumericVector mu,
                             double phi,
                             int d,
                             int var_power = 2,
                             double eps = 1e-6) {
  
  int n = z.size();
  
  if (mu.size() != n) {
    stop("length(mu) must equal length(z).");
  }
  
  if (!R_finite(phi) || phi <= 0.0) {
    stop("phi must be a positive scalar.");
  }
  
  if (d <= 0) {
    stop("d must be a positive integer.");
  }
  
  if (!(var_power == 1 || var_power == 2 || var_power == 3)) {
    stop("var_power must be one of 1, 2, or 3.");
  }
  
  if (!R_finite(eps) || eps <= 0.0) {
    stop("eps must be a positive scalar.");
  }
  
  NumericMatrix prob(n, d);
  NumericMatrix d_mu(n, d);
  NumericMatrix d_phi(n, d);
  
  // ----------------------------------------------------------
  // Special case d = 1
  // ----------------------------------------------------------
  
  if (d == 1) {
    
    for (int i = 0; i < n; ++i) {
      
      if (z[i] < 0) {
        stop("z must contain nonnegative integers.");
      }
      
      if (!R_finite(mu[i]) || mu[i] <= 0.0) {
        stop("mu must contain positive finite values.");
      }
      
      double mui = mu[i];
      double zi = static_cast<double>(z[i]);
      
      double shape_i = std::pow(mui, 2.0 - var_power) / phi;
      
      if (!R_finite(shape_i) || shape_i <= 0.0) {
        stop("Invalid Gamma shape implied by mu, phi, and var_power.");
      }
      
      double log_p = R::dnbinom_mu(zi, shape_i, mui, 1);
      double prob_i = std::exp(log_p);
      
      double dlog_dmu_fixed =
        zi / mui - (zi + shape_i) / (mui + shape_i);
      
      double d_mu_fixed_i = prob_i * dlog_dmu_fixed;
      
      double dlog_dshape =
        R::digamma(zi + shape_i) -
        R::digamma(shape_i) +
        std::log(shape_i / (shape_i + mui)) +
        1.0 -
        (shape_i + zi) / (shape_i + mui);
      
      double d_shape_i = prob_i * dlog_dshape;
      
      double d_mu_i =
        d_mu_fixed_i +
        d_shape_i * (2.0 - var_power) * shape_i / mui;
      
      double d_phi_i =
        d_shape_i * (-shape_i / phi);
      
      prob(i, 0) = prob_i;
      d_mu(i, 0) = d_mu_i;
      d_phi(i, 0) = d_phi_i;
    }
    
    return List::create(
      _["prob"] = prob,
      _["d_mu"] = d_mu,
      _["d_phi"] = d_phi
    );
  }
  
  // ----------------------------------------------------------
  // General case d >= 2
  // ----------------------------------------------------------
  
  for (int i = 0; i < n; ++i) {
    
    if (z[i] < 0) {
      stop("z must contain nonnegative integers.");
    }
    
    if (!R_finite(mu[i]) || mu[i] <= 0.0) {
      stop("mu must contain positive finite values.");
    }
    
    int zi_int = z[i];
    double zi = static_cast<double>(z[i]);
    double mui = mu[i];
    
    double shape_i = std::pow(mui, 2.0 - var_power) / phi;
    double r_i = shape_i / mui;
    double s_i = zi + shape_i;
    double A_i = r_i + 1.0;
    
    if (!R_finite(shape_i) || shape_i <= 0.0 ||
        !R_finite(r_i) || r_i <= 0.0) {
      stop("Invalid Gamma shape/rate implied by mu, phi, and var_power.");
    }
    
    double lgam_z1_i = R::lgammafn(zi + 1.0);
    double gamma_shape_i = R::gammafn(shape_i);
    double digamma_shape_i = R::digamma(shape_i);
    
    double eps_shape_i = std::min(eps, 0.5 * shape_i);
    
    if (!R_finite(eps_shape_i) || eps_shape_i <= 0.0) {
      stop("eps_shape contains non-positive values.");
    }
    
    for (int mm = 0; mm < d; ++mm) {
      
      double pL = static_cast<double>(mm) / static_cast<double>(d);
      double pU = static_cast<double>(mm + 1) / static_cast<double>(d);
      
      double L;
      double U;
      
      if (pL == 0.0) {
        L = 0.0;
      } else {
        L = R::qgamma(pL, shape_i, 1.0 / r_i, 1, 0);
      }
      
      if (pU == 1.0) {
        U = R_PosInf;
      } else {
        U = R::qgamma(pU, shape_i, 1.0 / r_i, 1, 0);
      }
      
      double AL = A_i * L;
      double AU = A_i * U;
      
      double rL = r_i * L;
      double rU = r_i * U;
      
      // ------------------------------
      // Lower bound
      // ------------------------------
      
      double gam_s_L =
        lower_gamma_scaled_cpp(s_i, AL, lgam_z1_i);
      
      double gam_sp1_L =
        lower_gamma_scaled_cpp(s_i + 1.0, AL, lgam_z1_i);
      
      double gam1_s_L =
        lower_gamma_scaled_dshape_cpp(
          s_i,
          AL,
          lgam_z1_i,
          eps_shape_i
        );
      
      double gam_shape_L_raw =
        lower_gamma_raw_cpp(shape_i, rL);
      
      double gam1_shape_L_raw =
        lower_gamma_raw_dshape_cpp(shape_i, rL, eps_shape_i);
      
      double hL =
        h_scaled_cpp(L, zi_int, lgam_z1_i);
      
      double h_dF_dshape_L = 0.0;
      double h_dF_dmu_L = 0.0;
      
      if (R_finite(L) && L > 0.0) {
        
        h_dF_dshape_L =
          hL / gamma_shape_i *
          (
            -digamma_shape_i * gam_shape_L_raw +
              gam1_shape_L_raw
          ) +
          std::exp(
            shape_i * std::log(r_i) +
              (shape_i + zi) * std::log(L) -
              A_i * L -
              lgam_z1_i -
              R::lgammafn(shape_i) -
              std::log(shape_i)
          );
        
        h_dF_dmu_L =
          -std::exp(
            shape_i * std::log(r_i) +
              (zi + shape_i) * std::log(L) -
              A_i * L -
              lgam_z1_i -
              R::lgammafn(shape_i) -
              std::log(mui)
          );
      }
      
      // ------------------------------
      // Upper bound
      // ------------------------------
      
      double gam_s_U =
        lower_gamma_scaled_cpp(s_i, AU, lgam_z1_i);
      
      double gam_sp1_U =
        lower_gamma_scaled_cpp(s_i + 1.0, AU, lgam_z1_i);
      
      double gam1_s_U =
        lower_gamma_scaled_dshape_cpp(
          s_i,
          AU,
          lgam_z1_i,
          eps_shape_i
        );
      
      double gam_shape_U_raw =
        lower_gamma_raw_cpp(shape_i, rU);
      
      double gam1_shape_U_raw =
        lower_gamma_raw_dshape_cpp(shape_i, rU, eps_shape_i);
      
      double hU =
        h_scaled_cpp(U, zi_int, lgam_z1_i);
      
      double h_dF_dshape_U = 0.0;
      double h_dF_dmu_U = 0.0;
      
      if (R_finite(U) && U > 0.0) {
        
        h_dF_dshape_U =
          hU / gamma_shape_i *
          (
            -digamma_shape_i * gam_shape_U_raw +
              gam1_shape_U_raw
          ) +
          std::exp(
            shape_i * std::log(r_i) +
              (shape_i + zi) * std::log(U) -
              A_i * U -
              lgam_z1_i -
              R::lgammafn(shape_i) -
              std::log(shape_i)
          );
        
        h_dF_dmu_U =
          -std::exp(
            shape_i * std::log(r_i) +
              (zi + shape_i) * std::log(U) -
              A_i * U -
              lgam_z1_i -
              R::lgammafn(shape_i) -
              std::log(mui)
          );
      }
      
      // ------------------------------
      // Interior terms
      // ------------------------------
      
      double Delta_s = gam_s_U - gam_s_L;
      double Delta_sp1 = gam_sp1_U - gam_sp1_L;
      double Delta1_s = gam1_s_U - gam1_s_L;

      double prob_i =
        std::pow(r_i, shape_i) /
        (gamma_shape_i * std::pow(A_i, s_i)) *
        Delta_s;
      
      double boundary_mu =
        -h_dF_dmu_U + h_dF_dmu_L;
      
      double integral_mu =
        (shape_i / (mui * mui)) *
        (std::pow(r_i, shape_i) / gamma_shape_i) *
        (
          Delta_sp1 / std::pow(A_i, s_i + 1.0) -
            mui * Delta_s / std::pow(A_i, s_i)
        );
      
      double d_mu_fixed_i =
        boundary_mu + integral_mu;
      
      double boundary_shape =
        -h_dF_dshape_U + h_dF_dshape_L;
      
      double A1 =
        (
          std::log(r_i) +
            1.0 -
            digamma_shape_i -
            std::log(A_i)
        ) * Delta_s;
      
      double A2 =
        Delta_sp1 / (mui * A_i);
      
      double A3 =
        Delta1_s;
      
      double integral_shape =
        std::pow(r_i, shape_i) /
        (gamma_shape_i * std::pow(A_i, s_i)) *
        (A1 - A2 + A3);
      
      double d_shape_i =
        boundary_shape + integral_shape;
      
      double d_mu_i =
        d_mu_fixed_i +
        d_shape_i * (2.0 - var_power) * shape_i / mui;
      
      double d_phi_i =
        d_shape_i * (-shape_i / phi);
      
      prob(i, mm) = prob_i;
      d_mu(i, mm) = d_mu_i;
      d_phi(i, mm) = d_phi_i;
    }
  }
  
  return List::create(
    _["prob"] = prob,
    _["d_mu"] = d_mu,
    _["d_phi"] = d_phi
  );
}
