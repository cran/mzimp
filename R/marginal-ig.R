# ============================================================
# I_ig
# Wrapper for IG mixed Poisson interval probabilities
# ============================================================

I_ig <- function(z, X, be, phi, d = 1L,
                 var_power = 3,
                 n_legendre = 32L,
                 log_out = FALSE,
                 rel_step_mu = 1e-6,
                 rel_step_phi = 1e-6,
                 eps_tail = 1e-300
                  ) {
  
  X <- as.matrix(X)
  be <- as.numeric(be)
  
  if (nrow(X) != length(z)) {
    stop("nrow(X) must equal length(z).")
  }
  if (length(be) != ncol(X)) {
    stop("length(be) must equal ncol(X).")
  }
  
  z <- as.integer(round(z))
  d <- as.integer(round(d))
  
  mu <- ig_safe_mu_from_Xbe(X, be)
  
  gl <- statmod::gauss.quad(n_legendre, kind = "legendre")

  out <- I_ig_mu_phi_cpp(
      z = z,
      mu = mu,
      phi = phi,
      d = d,
      var_power = var_power,
      leg_nodes = gl$nodes,
      leg_weights = gl$weights,
      log_out = log_out,
      rel_step_mu = rel_step_mu,
      rel_step_phi = rel_step_phi,
      eps_tail = eps_tail
    )

  
  rownames(out) <- paste0("obs_", seq_along(z))
  colnames(out) <- paste0("I_", seq_len(d))
  
  out
}

# ============================================================
# I_ig_grad
# Wrapper for IG mixed Poisson gradient
# ============================================================

I_ig_grad <- function(z, X, be, phi, d = 1L,
                      var_power = 3,
                      n_legendre = 32L,
                      rel_step_mu = 1e-6,
                      rel_step_phi = 1e-6,
                      eps_tail = 1e-300
                  ) {

  X <- as.matrix(X)
  be <- as.numeric(be)

  if (nrow(X) != length(z)) {
    stop("nrow(X) must equal length(z).")
  }
  if (length(be) != ncol(X)) {
    stop("length(be) must equal ncol(X).")
  }

  z <- as.integer(round(z))
  d <- as.integer(round(d))

  mu <- ig_safe_mu_from_Xbe(X, be)
  n  <- length(z)
  p  <- ncol(X)

  gl <- statmod::gauss.quad(n_legendre, kind = "legendre")

  gr <- I_ig_mu_phi_grad_cpp(
      z = z,
      mu = mu,
      phi = phi,
      d = d,
      var_power = var_power,
      leg_nodes = gl$nodes,
      leg_weights = gl$weights,
      rel_step_mu = rel_step_mu,
      rel_step_phi = rel_step_phi,
      eps_tail = eps_tail
    )

  prob  <- gr$prob
  d_mu  <- gr$d_mu
  d_phi <- gr$d_phi


  d_be <- array(NA_real_, dim = c(n, p, d))

  for (m in seq_len(d)) {
    d_be[, , m] <- X * as.numeric(mu * d_mu[, m])
  }

  rownames(prob) <- rownames(d_mu) <- rownames(d_phi) <-
    paste0("obs_", seq_len(n))

  colnames(prob) <- colnames(d_mu) <- colnames(d_phi) <-
    paste0("I_", seq_len(d))

  dimnames(d_be) <- list(
    paste0("obs_", seq_len(n)),
    colnames(X),
    paste0("I_", seq_len(d))
  )

  list(
    prob = prob,
    d_mu = d_mu,
    d_phi = d_phi,
    d_be = d_be
  )
}


# ============================================================
# ZI_ig
# Zero-inflated IG mixed Poisson
# ============================================================

ZI_ig <- function(y, X, be, Z, al, phi,
                  var_power = 3,
                  n_legendre = 32L,
                  rel_step_mu = 1e-6,
                  rel_step_phi = 1e-6,
                  eps_tail = 1e-300) {
  
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  be <- as.numeric(be)
  al <- as.numeric(al)
  
  if (length(y) != nrow(X)) stop("length(y) must equal nrow(X).")
  if (length(y) != nrow(Z)) stop("length(y) must equal nrow(Z).")
  if (length(be) != ncol(X)) stop("length(be) must equal ncol(X).")
  if (length(al) != ncol(Z)) stop("length(al) must equal ncol(Z).")
  
  y <- as.integer(round(y))
  
  eta_pr <- as.numeric(Z %*% al)
  pr <- plogis(eta_pr)
  
  base <- I_ig(
    z = y,
    X = X,
    be = be,
    phi = phi,
    d = 1L,
    var_power = var_power,
    n_legendre = n_legendre,
    log_out = FALSE,
    rel_step_mu = rel_step_mu,
    rel_step_phi = rel_step_phi,
    eps_tail = eps_tail
  )
  
  base <- as.numeric(base)
  
  out <- (1 - pr) * as.numeric(y == 0) + pr * base
  pmax(out, eps_tail)
}

# ============================================================
# ZI_ig
# Zero-inflated IG mixed Poisson
# ============================================================

ZI_ig <- function(y, X, be, Z, al, phi,
                  var_power = 3,
                  n_legendre = 32L,
                  rel_step_mu = 1e-6,
                  rel_step_phi = 1e-6,
                  eps_tail = 1e-300
                  ) {
  
  X <- as.matrix(X)
  Z <- as.matrix(Z)
  be <- as.numeric(be)
  al <- as.numeric(al)
  
  if (length(y) != nrow(X)) stop("length(y) must equal nrow(X).")
  if (length(y) != nrow(Z)) stop("length(y) must equal nrow(Z).")
  if (length(be) != ncol(X)) stop("length(be) must equal ncol(X).")
  if (length(al) != ncol(Z)) stop("length(al) must equal ncol(Z).")
  
  y <- as.integer(round(y))
  
  eta_pr <- as.numeric(Z %*% al)
  pr <- plogis(eta_pr)
  
  base <- I_ig(
    z = y,
    X = X,
    be = be,
    phi = phi,
    d = 1L,
    var_power = var_power,
    n_legendre = n_legendre,
    log_out = FALSE,
    rel_step_mu = rel_step_mu,
    rel_step_phi = rel_step_phi,
    eps_tail = eps_tail
  )
  
  base <- as.numeric(base)
  
  out <- (1 - pr) * as.numeric(y == 0) + pr * base
  pmax(out, eps_tail)
}

# ============================================================
# ZI_ig_grad
# Gradient of zero-inflated IG mixed Poisson
# ============================================================

ZI_ig_grad <- function(y, X, be, Z, al, phi,
                       var_power = 3,
                       n_legendre = 32L,
                       rel_step_mu = 1e-6,
                       rel_step_phi = 1e-6,
                       eps_tail = 1e-300) {

  X <- as.matrix(X)
  Z <- as.matrix(Z)
  be <- as.numeric(be)
  al <- as.numeric(al)

  if (length(y) != nrow(X)) stop("length(y) must equal nrow(X).")
  if (length(y) != nrow(Z)) stop("length(y) must equal nrow(Z).")
  if (length(be) != ncol(X)) stop("length(be) must equal ncol(X).")
  if (length(al) != ncol(Z)) stop("length(al) must equal ncol(Z).")

  y <- as.integer(round(y))

  eta_pr <- as.numeric(Z %*% al)
  pr <- plogis(eta_pr)

  # I_ig_grad() returns both probability and gradient.
  # Therefore I_ig() does not need to be evaluated separately.
  gr <- I_ig_grad(
    z = y,
    X = X,
    be = be,
    phi = phi,
    d = 1L,
    var_power = var_power,
    n_legendre = n_legendre,
    rel_step_mu = rel_step_mu,
    rel_step_phi = rel_step_phi,
    eps_tail = eps_tail
  )

  base <- as.numeric(gr$prob[, 1L])
  base <- pmax(base, eps_tail)

  db <- gr$d_be
  if (length(dim(db)) == 3L) {
    db <- db[, , 1L]
  }

  base_d_phi <- as.numeric(gr$d_phi[, 1L])

  prob <- (1 - pr) * as.numeric(y == 0) + pr * base
  prob <- pmax(prob, eps_tail)

  d_be  <- db * as.numeric(pr)
  d_phi <- pr * base_d_phi

  d_pr <- base - as.numeric(y == 0)
  d_al <- Z * as.numeric(d_pr * pr * (1 - pr))

  list(
    prob = prob,
    d_al = d_al,
    d_be = d_be,
    d_phi = d_phi
  )
}
