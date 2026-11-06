# ============================================================
# I_gamma
# Gamma mixed Poisson interval probabilities
# ============================================================

I_gamma <- function(z, X, be, phi, d,
                    var_power = 2,
                    log_out = FALSE) {
  
  X <- as.matrix(X)
  be <- as.numeric(be)
  
  if (nrow(X) != length(z)) {
    stop("nrow(X) must equal length(z).")
  }
  
  if (length(be) != ncol(X)) {
    stop("length(be) must equal ncol(X).")
  }
  
  if (!is.numeric(z) || any(!is.finite(z)) ||
      any(z < 0) || any(abs(z - round(z)) > 1e-12)) {
    stop("z must be a numeric vector of nonnegative integers.")
  }
  z <- as.integer(round(z))
  
  if (!is.numeric(phi) || length(phi) != 1 ||
      !is.finite(phi) || phi <= 0) {
    stop("phi must be a positive scalar.")
  }
  
  if (!is.numeric(d) || length(d) != 1 || !is.finite(d) ||
      d <= 0 || abs(d - round(d)) > 1e-12) {
    stop("d must be one positive integer.")
  }
  d <- as.integer(round(d))
  
  if (!is.numeric(var_power) ||
      length(var_power) != 1 ||
      !is.finite(var_power) ||
      abs(var_power - round(var_power)) > 1e-12 ||
      !(as.integer(round(var_power)) %in% c(1L, 2L, 3L))) {
    stop("var_power must be one of 1, 2, or 3.")
  }
  var_power <- as.integer(round(var_power))
  
  eta <- as.vector(X %*% be)
  mu <- exp(eta)
  
  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop("mu contains non-finite or non-positive values.")
  }
  
  out <- I_gamma_mu_phi_cpp(
    z = z,
    mu = mu,
    phi = phi,
    d = d,
    var_power = var_power,
    log_out = log_out
  )
  
  rownames(out) <- paste0("obs_", seq_along(z))
  colnames(out) <- paste0("I_", seq_len(d))
  
  out
}

# ============================================================
# I_gamma_grad_R
# Gradient with respect to beta and phi
# ============================================================

I_gamma_grad <- function(z, X, be, phi, d,
                           var_power = 2,
                           eps = 1e-6) {
  
  X <- as.matrix(X)
  be <- as.numeric(be)
  
  if (nrow(X) != length(z)) {
    stop("nrow(X) must equal length(z).")
  }
  
  if (length(be) != ncol(X)) {
    stop("length(be) must equal ncol(X).")
  }
  
  if (!is.numeric(z) || any(!is.finite(z)) ||
      any(z < 0) || any(abs(z - round(z)) > 1e-12)) {
    stop("z must be a numeric vector of nonnegative integers.")
  }
  z <- as.integer(round(z))
  
  if (!is.numeric(phi) || length(phi) != 1 || !is.finite(phi) || phi <= 0) {
    stop("phi must be a positive scalar.")
  }
  
  if (!is.numeric(d) || length(d) != 1 || !is.finite(d) ||
      d <= 0 || abs(d - round(d)) > 1e-12) {
    stop("d must be one positive integer.")
  }
  d <- as.integer(round(d))
  
  if (!is.numeric(var_power) ||
      length(var_power) != 1 ||
      !is.finite(var_power) ||
      abs(var_power - round(var_power)) > 1e-12 ||
      !(as.integer(round(var_power)) %in% c(1L, 2L, 3L))) {
    stop("var_power must be one of 1, 2, or 3.")
  }
  var_power <- as.integer(round(var_power))
  
  if (!is.numeric(eps) || length(eps) != 1 || !is.finite(eps) || eps <= 0) {
    stop("eps must be a positive scalar.")
  }
  
  eta <- as.vector(X %*% be)
  mu <- exp(eta)
  
  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop("mu contains non-finite or non-positive values.")
  }
  
  n <- length(z)
  p <- ncol(X)
  g <- I_gamma_mu_phi_grad_cpp(
    z = z,
    mu = mu,
    phi = phi,
    d = d,
    var_power = var_power,
    eps = eps
  )
  prob <- g$prob
  d_mu <- g$d_mu
  d_phi <- g$d_phi
  
  d_be <- array(0, dim = c(n, p, d))
  
  for (m in seq_len(d)) {
    d_be[, , m] <- X * as.numeric(mu * d_mu[, m])
  }
  
  rownames(prob) <- paste0("obs_", seq_len(n))
  rownames(d_mu) <- paste0("obs_", seq_len(n))
  rownames(d_phi) <- paste0("obs_", seq_len(n))
  
  colnames(prob) <- paste0("I_", seq_len(d))
  colnames(d_mu) <- paste0("I_", seq_len(d))
  colnames(d_phi) <- paste0("I_", seq_len(d))
  
  dimnames(d_be) <- list(
    paste0("obs_", seq_len(n)),
    colnames(X),
    paste0("I_", seq_len(d))
  )
  
  d_be_total <- apply(d_be, 2, sum)
  d_be_by_interval <- apply(d_be, c(2, 3), sum)
  
  list(
    prob = prob,
    d_mu = d_mu,
    d_phi = d_phi,
    d_be = d_be,
    d_be_total = d_be_total,
    d_be_by_interval = d_be_by_interval
  )
}


# ============================================================
# ZI_gamma_R
# Zero-inflated Gamma mixed Poisson
# ============================================================

ZI_gamma <- function(y, X, be, Z, al, phi,
                       var_power = 2) {
  
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  be <- as.numeric(be)
  al <- as.numeric(al)
  
  if (length(y) != nrow(X)) {
    stop("length(y) must equal nrow(X).")
  }
  
  if (length(y) != nrow(Z)) {
    stop("length(y) must equal nrow(Z).")
  }
  
  if (length(be) != ncol(X)) {
    stop("length(be) must equal ncol(X).")
  }
  
  if (length(al) != ncol(Z)) {
    stop("length(al) must equal ncol(Z).")
  }
  
  if (!is.numeric(var_power) ||
      length(var_power) != 1 ||
      !is.finite(var_power) ||
      abs(var_power - round(var_power)) > 1e-12 ||
      !(as.integer(round(var_power)) %in% c(1L, 2L, 3L))) {
    stop("var_power must be one of 1, 2, or 3.")
  }
  var_power <- as.integer(round(var_power))
  
  eta_pr <- as.numeric(Z %*% al)
  pr <- plogis(eta_pr)
  
  base <- I_gamma(
    z = y,
    X = X,
    be = be,
    phi = phi,
    d = 1,
    var_power = var_power,
    log_out = FALSE
  )
  
  base <- as.numeric(base)
  
  out <- (1 - pr) * as.numeric(y == 0) + pr * base
  
  out
}




# ============================================================
# ZI_gamma_grad_R
# Gradient of zero-inflated Gamma mixed Poisson
# ============================================================

# ============================================================
# ZI_gamma_grad
# Gradient of zero-inflated Gamma mixed Poisson
# ============================================================

ZI_gamma_grad <- function(y, X, be, Z, al, phi,
                          var_power = 2,
                          eps = 1e-6) {
  
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  be <- as.numeric(be)
  al <- as.numeric(al)
  
  if (length(y) != nrow(X)) {
    stop("length(y) must equal nrow(X).")
  }
  
  if (length(y) != nrow(Z)) {
    stop("length(y) must equal nrow(Z).")
  }
  
  if (length(be) != ncol(X)) {
    stop("length(be) must equal ncol(X).")
  }
  
  if (length(al) != ncol(Z)) {
    stop("length(al) must equal ncol(Z).")
  }
  
  if (!is.numeric(var_power) ||
      length(var_power) != 1 ||
      !is.finite(var_power) ||
      abs(var_power - round(var_power)) > 1e-12 ||
      !(as.integer(round(var_power)) %in% c(1L, 2L, 3L))) {
    stop("var_power must be one of 1, 2, or 3.")
  }
  
  var_power <- as.integer(round(var_power))
  
  # Probability of belonging to the count component
  eta_pr <- as.numeric(Z %*% al)
  pr <- plogis(eta_pr)
  
  # Gradient of the Gamma-mixed Poisson component
  gr <- I_gamma_grad(
    z = y,
    X = X,
    be = be,
    phi = phi,
    d = 1,
    var_power = var_power,
    eps = eps
  )
  
  # Base probability is already computed by I_gamma_grad()
  base <- as.numeric(gr$prob[, 1L])
  
  # ----------------------------------------------------------
  # Gradient with respect to beta
  #
  # dP/dbeta = pr * df/dbeta
  # ----------------------------------------------------------
  
  d_be <- pr * gr$d_be[, , 1L]
  
  # ----------------------------------------------------------
  # Gradient with respect to phi
  #
  # dP/dphi = pr * df/dphi
  # ----------------------------------------------------------
  
  d_phi <- pr * as.numeric(gr$d_phi[, 1L])
  
  # ----------------------------------------------------------
  # Gradient with respect to alpha
  #
  # P = (1-pr) * 1(y=0) + pr * base
  #
  # dP/dpr = base - 1(y=0)
  # dpr/dalpha = pr * (1-pr) * Z
  # ----------------------------------------------------------
  
  d_pr <- base - as.numeric(y == 0)
  
  d_al <- Z * as.numeric(
    d_pr * pr * (1 - pr)
  )
  
  # Full zero-inflated probability
  prob <- (1 - pr) * as.numeric(y == 0) + pr * base
  
  # Quantities required by marginal estimation
  list(
    prob = prob,
    d_al = d_al,
    d_be = d_be,
    d_phi = d_phi
  )
}