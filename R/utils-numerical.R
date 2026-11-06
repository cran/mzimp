# ============================================================
# Stable helper functions for gamma function
# ============================================================
logspace_sub <- function(logx, logy, tol = 1e-12) {
  out <- rep(NA_real_, length(logx))
  
  idx1 <- is.finite(logx) & is.infinite(logy) & (logy < 0)
  out[idx1] <- logx[idx1]
  
  idx2 <- is.finite(logx) & is.finite(logy)
  if (any(idx2)) {
    lx <- logx[idx2]
    ly <- logy[idx2]
    
    bad <- lx < ly - tol
    if (any(bad)) {
      stop("logspace_sub received logx < logy.")
    }
    
    eq <- abs(lx - ly) <= tol
    tmp <- rep(NA_real_, length(lx))
    tmp[eq] <- -Inf
    
    gt <- !eq
    tmp[gt] <- lx[gt] + log1p(-exp(ly[gt] - lx[gt]))
    
    out[idx2] <- tmp
  }
  
  idx3 <- is.infinite(logx) & (logx < 0) &
    is.infinite(logy) & (logy < 0)
  out[idx3] <- -Inf
  
  out
}


stable_pdiff <- function(yL, yU, shape_vec) {
  if (length(yL) != length(yU) || length(yL) != length(shape_vec)) {
    stop("yL, yU, and shape_vec must have the same length.")
  }
  
  if (any(yU < yL, na.rm = TRUE)) {
    stop("Found yU < yL in stable_pdiff.")
  }
  
  pU <- pgamma(yU, shape = shape_vec, rate = 1)
  use_lower <- (pU <= 0.5)
  
  out_log <- rep(NA_real_, length(shape_vec))
  
  if (any(use_lower)) {
    idx <- which(use_lower)
    logFU <- pgamma(yU[idx], shape = shape_vec[idx], rate = 1, log.p = TRUE)
    logFL <- pgamma(yL[idx], shape = shape_vec[idx], rate = 1, log.p = TRUE)
    out_log[idx] <- logspace_sub(logFU, logFL)
  }
  
  if (any(!use_lower)) {
    idx <- which(!use_lower)
    logSU <- pgamma(yU[idx], shape = shape_vec[idx], rate = 1,
                    lower.tail = FALSE, log.p = TRUE)
    logSL <- pgamma(yL[idx], shape = shape_vec[idx], rate = 1,
                    lower.tail = FALSE, log.p = TRUE)
    out_log[idx] <- logspace_sub(logSL, logSU)
  }
  
  out_log
}


lower_gamma_scaled <- function(shape_vec, x, lgam_denom) {
  out <- numeric(length(shape_vec))
  
  idx_inf <- is.infinite(x)
  if (any(idx_inf)) {
    out[idx_inf] <- exp(lgamma(shape_vec[idx_inf]) - lgam_denom[idx_inf])
  }
  
  idx_pos <- (!idx_inf) & (x > 0)
  if (any(idx_pos)) {
    logP <- pgamma(x[idx_pos], shape = shape_vec[idx_pos],
                   rate = 1, log.p = TRUE)
    out[idx_pos] <- exp(lgamma(shape_vec[idx_pos]) -
                          lgam_denom[idx_pos] + logP)
  }
  
  out
}


lower_gamma_scaled_dshape <- function(shape_vec, x, lgam_denom, eps_shape) {
  out <- numeric(length(shape_vec))
  
  idx_inf <- is.infinite(x)
  if (any(idx_inf)) {
    out[idx_inf] <- exp(lgamma(shape_vec[idx_inf]) - lgam_denom[idx_inf]) *
      digamma(shape_vec[idx_inf])
  }
  
  idx_pos <- (!idx_inf) & (x > 0)
  if (any(idx_pos)) {
    up <- lower_gamma_scaled(
      shape_vec[idx_pos] + eps_shape,
      x[idx_pos],
      lgam_denom[idx_pos]
    )
    
    dn <- lower_gamma_scaled(
      shape_vec[idx_pos] - eps_shape,
      x[idx_pos],
      lgam_denom[idx_pos]
    )
    
    out[idx_pos] <- (up - dn) / (2 * eps_shape)
  }
  
  out
}

h_scaled <- function(x, z, lgam_denom) {
  out <- numeric(length(z))
  idx <- (!is.infinite(x)) & (x > 0)
  if (any(idx)) {
    out[idx] <- exp(-x[idx] + z[idx] * log(x[idx]) - lgam_denom[idx])
  }  
  out
}


######################
# FUNCTION FOR ỔN ĐỊNH SỐ
########################

safe_plogis <- function(x, lower_eta = -30, upper_eta = 30) {
  x <- pmin(pmax(as.numeric(x), lower_eta), upper_eta)
  plogis(x)
}

# Chặn eta để exp(eta) không overflow / underflow
# Có thể chỉnh lower/upper nếu muốn
safe_mu_from_Xbe <- function(X, be, lower_eta = -30, upper_eta = 30) {
  eta <- as.numeric(as.matrix(X) %*% be)
  eta <- pmin(pmax(eta, lower_eta), upper_eta)
  exp(eta)
}


############################
# ============================================================
# Helpers
# ============================================================

lnorm_gh_rule <- function(n_gh) {
  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop("Package 'statmod' is required.")
  }
  statmod::gauss.quad(as.integer(n_gh), kind = "hermite")
}

lnorm_row_logsumexp <- function(A) {
  rmax <- apply(A, 1L, max)
  rmax + log(rowSums(exp(A - rmax)))
}

lnorm_weighted_logsum <- function(logf, w) {
  A <- sweep(logf, 2L, log(w), "+")
  lnorm_row_logsumexp(A)
}

lnorm_weighted_rowsum <- function(K, w) {
  rowSums(sweep(K, 2L, w, "*"))
}

lnorm_log_integrand <- function(u, z, mu, phi,
                                var_power = 2,
                                eps_u = 1e-16) {
  
  u <- pmin(pmax(u, eps_u), 1 - eps_u)
  
  sigma2 <- log1p(phi * mu^(var_power - 2))
  sigma  <- sqrt(sigma2)
  
  qn <- qnorm(u)
  
  logQ <- outer(log(mu) - 0.5 * sigma2, rep(1, length(u))) +
    sigma %o% qn
  
  z_mat <- matrix(z, nrow = length(z), ncol = length(u))
  lgam_mat <- matrix(lgamma(z + 1), nrow = length(z), ncol = length(u))
  
  zlogQ <- z_mat * logQ
  zlogQ[z_mat == 0L] <- 0
  
  logf <- zlogQ - exp(pmin(logQ, 700)) - lgam_mat
  logf[logQ > 700] <- -Inf
  logf[!is.finite(logf)] <- -Inf
  
  logf
}

lnorm_eval_kernels <- function(u, z, mu, phi,
                               var_power = 2,
                               eps_u = 1e-16) {
  
  u <- pmin(pmax(u, eps_u), 1 - eps_u)
  
  n <- length(z)
  qn <- qnorm(u)
  
  sigma2 <- log1p(phi * mu^(var_power - 2))
  sigma  <- sqrt(sigma2)
  
  logQ <- outer(log(mu) - 0.5 * sigma2, rep(1, length(u))) +
    sigma %o% qn
  
  z_mat <- matrix(z, nrow = n, ncol = length(u))
  lgam_mat <- matrix(lgamma(z + 1), nrow = n, ncol = length(u))
  
  zlogQ <- z_mat * logQ
  zlogQ[z_mat == 0L] <- 0
  
  logf <- zlogQ - exp(pmin(logQ, 700)) - lgam_mat
  logf[logQ > 700] <- -Inf
  logf[!is.finite(logf)] <- -Inf
  
  f <- exp(logf)
  Q <- exp(pmin(logQ, 700))
  A <- z_mat - Q
  
  B <- phi * mu^(var_power - 2)
  d_sigma2_dmu  <- phi * (var_power - 2) * mu^(var_power - 3) / (1 + B)
  d_sigma2_dphi <- mu^(var_power - 2) / (1 + B)
  
  d_sigma_dmu  <- 0.5 * d_sigma2_dmu / sigma
  d_sigma_dphi <- 0.5 * d_sigma2_dphi / sigma
  
  dlogQ_dmu <- outer(
    1 / mu - 0.5 * d_sigma2_dmu,
    rep(1, length(u))
  ) + d_sigma_dmu %o% qn
  
  dlogQ_dphi <- outer(
    -0.5 * d_sigma2_dphi,
    rep(1, length(u))
  ) + d_sigma_dphi %o% qn
  
  dmu  <- f * A * dlogQ_dmu
  dphi <- f * A * dlogQ_dphi
  
  list(
    value = f,
    dmu = dmu,
    dphi = dphi
  )
}


# ============================================================
# Inverse Gaussian:
#   E(G)   = mu
#   Var(G) = mu^3 * phi
#
# statmod convention:
#   shape = lambda = 1 / phi
# ============================================================

# ============================================================
# 1. Helpers
# ============================================================

ig_safe_mu_from_Xbe <- function(X, be, lower_eta = -30, upper_eta = 30) {
  eta <- as.numeric(as.matrix(X) %*% be)
  eta <- pmin(pmax(eta, lower_eta), upper_eta)
  exp(eta)
}

ig_safe_plogis <- function(x, lower_eta = -30, upper_eta = 30) {
  x <- pmin(pmax(as.numeric(x), lower_eta), upper_eta)
  plogis(x)
}

ig_row_logsumexp <- function(log_mat) {
  if (!is.matrix(log_mat)) log_mat <- as.matrix(log_mat)
  m <- apply(log_mat, 1L, max)
  m[!is.finite(m)] <- 0
  s <- rowSums(exp(log_mat - m))
  out <- m + log(s)
  out[s <= 0] <- -Inf
  out
}

ig_row_logsumexp_signed <- function(log_abs_mat, sign_mat) {
  if (!is.matrix(log_abs_mat)) log_abs_mat <- as.matrix(log_abs_mat)
  if (!is.matrix(sign_mat))    sign_mat    <- as.matrix(sign_mat)
  
  m <- apply(log_abs_mat, 1L, max)
  m[!is.finite(m)] <- 0
  
  scaled <- sign_mat * exp(log_abs_mat - m)
  s <- rowSums(scaled)
  
  out_sign <- sign(s)
  out_log_abs <- m + log(abs(s))
  out_log_abs[s == 0] <- -Inf
  out_sign[s == 0] <- 0
  
  list(sign = out_sign, log_abs = out_log_abs)
}

ig_logdiffexp_vec <- function(a, b) {
  out <- a + log1p(-exp(b - a))
  out[is.infinite(b) & b < 0] <- a[is.infinite(b) & b < 0]
  out
}

ig_gauss_legendre <- function(n_nodes) {
  gl <- statmod::gauss.quad(n_nodes, kind = "legendre")
  list(nodes = gl$nodes, weights = gl$weights)
}

ig_gauss_laguerre <- function(n_nodes) {
  gl <- statmod::gauss.quad(n_nodes, kind = "laguerre")
  list(nodes = gl$nodes, weights = gl$weights)
}




log_besselK_half_integer <- function(x, z) {
  logsumexp <- function(a) {
    m <- max(a)
    if (!is.finite(m)) return(m)
    m + log(sum(exp(a - m)))
  }
  x <- as.numeric(x)
  z <- as.integer(round(z))
  
  if (length(x) != length(z)) stop("length(x) must equal length(z).")
  if (any(x <= 0) || any(!is.finite(x))) stop("x must be positive.")
  if (any(z < 0) || any(!is.finite(z))) stop("z must be nonnegative integers.")
  
  out <- numeric(length(x))
  
  for (i in seq_along(x)) {
    xi <- x[i]
    zi <- z[i]
    
    # K_{z - 1/2}(x) = K_{|z - 1/2|}(x)
    # For z = 0: K_{-1/2}(x) = K_{1/2}(x), so m = 0
    # For z >= 1: z - 1/2 = (z - 1) + 1/2, so m = z - 1
    m <- if (zi == 0L) 0L else zi - 1L
    
    if (m == 0L) {
      log_poly <- 0
    } else {
      k <- 0:m
      
      log_terms <-
        lgamma(m + k + 1) -
        lgamma(k + 1) -
        lgamma(m - k + 1) -
        k * log(2 * xi)
      
      log_poly <- logsumexp(log_terms)
    }
    
    out[i] <-
      0.5 * (log(pi) - log(2 * xi)) -
      xi +
      log_poly
  }
  
  out
}










# ============================================================
# 2. IG quantile / density / score
# ============================================================

ig_qinvgauss_mu_phi <- function(u, mu, phi) {
  if (any(!is.finite(u)) || any(u <= 0) || any(u >= 1)) {
    stop("u must lie strictly inside (0,1).")
  }
  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop("mu must be positive.")
  }
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi) || phi <= 0) {
    stop("phi must be a positive scalar.")
  }
  
  statmod::qinvgauss(u, mean = mu, shape = 1 / phi)
}

ig_log_dinvgauss_mu_phi <- function(x, mu, phi, eps_x = 1e-300) {
  x <- pmax(x, eps_x)
  
  if (is.matrix(x)) {
    mu_mat <- matrix(mu, nrow = nrow(x), ncol = ncol(x))
  } else {
    mu_mat <- mu
  }
  
  -0.5 * log(2 * pi * phi) -
    1.5 * log(x) -
    (x - mu_mat)^2 / (2 * phi * mu_mat^2 * x)
}

ig_score_mu_logf <- function(x, mu, phi, eps_x = 1e-300) {
  x <- pmax(x, eps_x)
  
  if (is.matrix(x)) {
    mu_mat <- matrix(mu, nrow = nrow(x), ncol = ncol(x))
  } else {
    mu_mat <- mu
  }
  
  (x - mu_mat) / (phi * mu_mat^3)
}

ig_score_phi_logf <- function(x, mu, phi, eps_x = 1e-300) {
  x <- pmax(x, eps_x)
  
  if (is.matrix(x)) {
    mu_mat <- matrix(mu, nrow = nrow(x), ncol = ncol(x))
  } else {
    mu_mat <- mu
  }
  
  -1 / (2 * phi) +
    (x - mu_mat)^2 / (2 * phi^2 * mu_mat^2 * x)
}

# ============================================================
# 3. Numerical derivatives of quantile bounds
#    vectorized over observations and probabilities
# ============================================================

ig_qbounds_mu_phi_grad <- function(probs, mu, phi,
                                   rel_step_mu = 1e-6,
                                   rel_step_phi = 1e-6) {
  probs <- as.numeric(probs)
  mu    <- as.numeric(mu)
  
  if (any(!is.finite(probs)) || any(probs <= 0) || any(probs >= 1)) {
    stop("probs must lie strictly inside (0,1).")
  }
  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop("mu must be positive.")
  }
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi) || phi <= 0) {
    stop("phi must be a positive scalar.")
  }
  
  h_mu <- pmax(rel_step_mu * mu, 1e-8)
  mu_plus  <- mu + h_mu
  mu_minus <- pmax(mu - h_mu, 1e-12)
  
  Q <- outer(mu, probs, function(m, p) {
    statmod::qinvgauss(p, mean = m, shape = 1 / phi)
  })
  
  Q_plus_mu <- outer(mu_plus, probs, function(m, p) {
    statmod::qinvgauss(p, mean = m, shape = 1 / phi)
  })
  
  Q_minus_mu <- outer(mu_minus, probs, function(m, p) {
    statmod::qinvgauss(p, mean = m, shape = 1 / phi)
  })
  
  dQ_dmu <- (Q_plus_mu - Q_minus_mu) / (mu_plus - mu_minus)
  
  h_phi <- max(rel_step_phi * phi, 1e-8)
  phi_plus  <- phi + h_phi
  phi_minus <- max(phi - h_phi, 1e-10)
  
  Q_plus_phi <- outer(mu, probs, function(m, p) {
    statmod::qinvgauss(p, mean = m, shape = 1 / phi_plus)
  })
  
  Q_minus_phi <- outer(mu, probs, function(m, p) {
    statmod::qinvgauss(p, mean = m, shape = 1 / phi_minus)
  })
  
  dQ_dphi <- (Q_plus_phi - Q_minus_phi) / (phi_plus - phi_minus)
  
  list(
    Q = Q,
    dQ_dmu = dQ_dmu,
    dQ_dphi = dQ_dphi
  )
}

# ============================================================
# 4. Closed form for I_ig(z; 0,1)
# ============================================================

ig_I_d1_closed_log <- function(z, mu, phi) {
  z  <- as.integer(round(z))
  mu <- as.numeric(mu)
  
  if (length(mu) != length(z)) stop("length(mu) must equal length(z).")
  if (any(mu <= 0) || any(!is.finite(mu))) stop("mu must be positive.")
  if (!is.numeric(phi) || length(phi) != 1L || !is.finite(phi) || phi <= 0) {
    stop("phi must be a positive scalar.")
  }
  
  nu <- z - 0.5
  r  <- mu^2 / (1 + 2 * phi * mu^2)
  x0 <- sqrt(1 + 2 * phi * mu^2) / (phi * mu)
  
  logK <- log_besselK_half_integer(x = x0, z = z)
  
  logI <-
    0.5 * log(2 / (pi * phi)) -
    lgamma(z + 1) +
    1 / (phi * mu) +
    0.5 * nu * log(r) +
    logK
  
  pmin(logI, 0)
}

ig_I_d1_closed <- function(z, X, be, phi, return_log = FALSE) {
  X <- as.matrix(X)
  
  if (nrow(X) != length(z)) stop("nrow(X) must equal length(z).")
  if (length(be) != ncol(X)) stop("length(be) must equal ncol(X).")
  if (!is.numeric(z) || any(!is.finite(z)) ||
      any(z < 0) || any(abs(z - round(z)) > 1e-12)) {
    stop("z must be nonnegative integers.")
  }
  
  mu <- ig_safe_mu_from_Xbe(X, be)
  out_log <- ig_I_d1_closed_log(z = z, mu = mu, phi = phi)
  
  if (return_log) out_log else exp(out_log)
}

# ============================================================
# 5. Core log h(x)
#    h(x) = exp(-x) * x^z / Gamma(z+1) * f_IG(x)
# ============================================================

ig_log_h_matrix <- function(x_mat, z, mu, phi, eps_x = 1e-300) {
  x_mat <- pmax(x_mat, eps_x)
  
  logf <- ig_log_dinvgauss_mu_phi(x_mat, mu, phi, eps_x = eps_x)
  
  sweep(log(x_mat), 1L, z, "*") -
    x_mat -
    lgamma(z + 1) +
    logf
}

# ============================================================
# 6. Finite interval [L,U] by Gauss-Legendre
#    with Leibniz rule
# ============================================================

ig_I_finite_legendre_value_grad <- function(z, mu, phi,
                                            L, U,
                                            dL_dmu = NULL,
                                            dU_dmu = NULL,
                                            dL_dphi = NULL,
                                            dU_dphi = NULL,
                                            n_nodes = 32L,
                                            eps_x = 1e-300) {
  z <- as.integer(round(z))
  n <- length(z)
  
  if (length(mu) != n || length(L) != n || length(U) != n) {
    stop("z, mu, L, U must have same length.")
  }
  
  if (is.null(dL_dmu))   dL_dmu   <- rep(0, n)
  if (is.null(dU_dmu))   dU_dmu   <- rep(0, n)
  if (is.null(dL_dphi))  dL_dphi  <- rep(0, n)
  if (is.null(dU_dphi))  dU_dphi  <- rep(0, n)
  
  gl <- ig_gauss_legendre(n_nodes)
  xi <- gl$nodes
  wi <- gl$weights
  
  half <- (U - L) / 2
  mid  <- (U + L) / 2
  
  x_mat <- outer(half, xi, "*") + matrix(mid, nrow = n, ncol = n_nodes)
  log_h <- ig_log_h_matrix(x_mat, z, mu, phi, eps_x = eps_x)
  
  log_dxw <- log(outer(half, wi))
  log_terms <- log_h + log_dxw
  
  value_log <- ig_row_logsumexp(log_terms)
  value <- exp(value_log)
  
  score_mu <- ig_score_mu_logf(x_mat, mu, phi, eps_x = eps_x)
  log_abs_mu <- log_terms + log(abs(score_mu))
  log_abs_mu[!is.finite(log_abs_mu)] <- -Inf
  sign_mu <- sign(score_mu)
  tmp_mu <- ig_row_logsumexp_signed(log_abs_mu, sign_mu)
  int_dmu <- tmp_mu$sign * exp(tmp_mu$log_abs)
  
  score_phi <- ig_score_phi_logf(x_mat, mu, phi, eps_x = eps_x)
  log_abs_phi <- log_terms + log(abs(score_phi))
  log_abs_phi[!is.finite(log_abs_phi)] <- -Inf
  sign_phi <- sign(score_phi)
  tmp_phi <- ig_row_logsumexp_signed(log_abs_phi, sign_phi)
  int_dphi <- tmp_phi$sign * exp(tmp_phi$log_abs)
  
  log_hU <- as.numeric(ig_log_h_matrix(matrix(U, nrow = n, ncol = 1L),
                                       z, mu, phi, eps_x = eps_x))
  hU <- exp(log_hU)
  
  need_L_mu   <- abs(dL_dmu)  > 0
  need_L_phi  <- abs(dL_dphi) > 0
  need_L_any  <- need_L_mu | need_L_phi
  
  hL <- rep(0, n)
  if (any(need_L_any)) {
    log_hL <- as.numeric(ig_log_h_matrix(matrix(pmax(L, eps_x), nrow = n, ncol = 1L),
                                         z, mu, phi, eps_x = eps_x))
    hL[need_L_any] <- exp(log_hL[need_L_any])
  }
  
  d_mu  <- int_dmu  + hU * dU_dmu  - hL * dL_dmu
  d_phi <- int_dphi + hU * dU_dphi - hL * dL_dphi
  
  list(
    value = value,
    d_mu  = d_mu,
    d_phi = d_phi
  )
}

# ============================================================
# 7. Tail interval [L, +inf) by shifted Gauss-Laguerre
#    with Leibniz rule
# ============================================================

ig_I_tail_laguerre_value_grad <- function(z, mu, phi,
                                          L,
                                          dL_dmu = NULL,
                                          dL_dphi = NULL,
                                          n_nodes = 32L,
                                          eps_x = 1e-300) {
  z <- as.integer(round(z))
  n <- length(z)
  
  if (length(mu) != n || length(L) != n) {
    stop("z, mu, L must have same length.")
  }
  
  if (is.null(dL_dmu))  dL_dmu  <- rep(0, n)
  if (is.null(dL_dphi)) dL_dphi <- rep(0, n)
  
  gl <- ig_gauss_laguerre(n_nodes)
  t  <- gl$nodes
  w  <- gl$weights
  
  x_mat <- outer(L, t, "+")
  
  logf <- ig_log_dinvgauss_mu_phi(x_mat, mu, phi, eps_x = eps_x)
  log_core <- sweep(log(pmax(x_mat, eps_x)), 1L, z, "*") -
    lgamma(z + 1) +
    logf
  
  log_terms <- log_core +
    matrix(log(w), nrow = n, ncol = n_nodes, byrow = TRUE) -
    matrix(L, nrow = n, ncol = n_nodes)
  
  value_log <- ig_row_logsumexp(log_terms)
  value <- exp(value_log)
  
  score_mu <- ig_score_mu_logf(x_mat, mu, phi, eps_x = eps_x)
  log_abs_mu <- log_terms + log(abs(score_mu))
  log_abs_mu[!is.finite(log_abs_mu)] <- -Inf
  sign_mu <- sign(score_mu)
  tmp_mu <- ig_row_logsumexp_signed(log_abs_mu, sign_mu)
  int_dmu <- tmp_mu$sign * exp(tmp_mu$log_abs)
  
  score_phi <- ig_score_phi_logf(x_mat, mu, phi, eps_x = eps_x)
  log_abs_phi <- log_terms + log(abs(score_phi))
  log_abs_phi[!is.finite(log_abs_phi)] <- -Inf
  sign_phi <- sign(score_phi)
  tmp_phi <- ig_row_logsumexp_signed(log_abs_phi, sign_phi)
  int_dphi <- tmp_phi$sign * exp(tmp_phi$log_abs)
  
  need_L_mu   <- abs(dL_dmu)  > 0
  need_L_phi  <- abs(dL_dphi) > 0
  need_L_any  <- need_L_mu | need_L_phi
  
  hL <- rep(0, n)
  if (any(need_L_any)) {
    log_hL <- as.numeric(ig_log_h_matrix(matrix(pmax(L, eps_x), nrow = n, ncol = 1L),
                                         z, mu, phi, eps_x = eps_x))
    hL[need_L_any] <- exp(log_hL[need_L_any])
  }
  
  d_mu  <- int_dmu  - hL * dL_dmu
  d_phi <- int_dphi - hL * dL_dphi
  
  list(
    value = value,
    d_mu  = d_mu,
    d_phi = d_phi
  )
}


