
# ============================================================
# I_lnorm
# ============================================================

I_lnorm <- function(z, X, be, phi, d = 1,
                    var_power = 2,
                    n_gh = 20,
                    gh = NULL,
                    log_out = FALSE
                    ) {
  
  X <- as.matrix(X)
  z <- as.integer(z)
  
  mu <- safe_mu_from_Xbe(X, be)
  
  if (is.null(gh)) {
    gh <- statmod::gauss.quad(n_gh, kind = "hermite")
  }
  
  I_lnorm_mu_phi_cpp(
      z = z,
      mu = mu,
      phi = phi,
      d = d,
      var_power = var_power,
      gh_nodes = gh$nodes,
      gh_weights = gh$weights,
      log_out = log_out
    )
  
}




# ============================================================
# I_lnorm_grad
# ============================================================


I_lnorm_grad <- function(z, X, be, phi, d = 1L,
                         var_power = 2,
                         n_gh = 20,
                         gh = NULL,
                         eps_u = 1e-16
                         ) {

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

  if (!is.numeric(eps_u) || length(eps_u) != 1 ||
      !is.finite(eps_u) || eps_u <= 0) {
    stop("eps_u must be a positive scalar.")
  }
  eta <- as.vector(X %*% be)
  mu <- exp(eta)

  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop("mu contains non-finite or non-positive values.")
  }

  n <- length(z)
  p <- ncol(X)

  if (is.null(gh)) {
    gh <- statmod::gauss.quad(n_gh, kind = "hermite")
  }

  if (is.null(gh$nodes) || is.null(gh$weights)) {
    stop("gh must be a list with components 'nodes' and 'weights'.")
  }

  g <- I_lnorm_mu_phi_grad_cpp(
    z = z,
    mu = mu,
    phi = phi,
    d = d,
    var_power = var_power,
    gh_nodes = gh$nodes,
    gh_weights = gh$weights,
    eps_u = eps_u
  )

  prob <- g$prob
  d_mu <- g$d_mu
  d_phi <- g$d_phi

  d_be <- array(0, dim = c(n, p, d))

  for (m in seq_len(d)) {
    d_be[, , m] <- X * as.numeric(mu * d_mu[, m])
  }

  obs_names <- paste0("obs_", seq_len(n))
  I_names <- paste0("I_", seq_len(d))

  rownames(prob) <- obs_names
  colnames(prob) <- I_names

  rownames(d_mu) <- obs_names
  colnames(d_mu) <- I_names

  rownames(d_phi) <- obs_names
  colnames(d_phi) <- I_names

  dimnames(d_be) <- list(
    obs_names,
    colnames(X),
    I_names
  )

  list(
    prob = prob,
    d_mu = d_mu,
    d_phi = d_phi,
    d_be = d_be
  )
}


# ============================================================
# ZI_lnorm
# ============================================================

ZI_lnorm <- function(y, X, Z, be, phi, al,
                     var_power = 2,
                     n_gh = 20,
                     gh = NULL,
                     eps_prob = 1e-300) {
  
  pr <- safe_plogis(as.numeric(as.matrix(Z) %*% al))
  
  I0 <- I_lnorm(
    z = y, X = X, be = be, phi = phi, d = 1,
    var_power = var_power,
    n_gh = n_gh,
    gh = gh
  )
  
  I0 <- as.numeric(I0)
  
  out <- pr * I0
  idx0 <- (y == 0)
  out[idx0] <- (1 - pr[idx0]) + pr[idx0] * I0[idx0]
  
  pmax(out, eps_prob)
}

# ============================================================
# ZI_lnorm_grad
# ============================================================

ZI_lnorm_grad <- function(y, X, Z, be, phi, al,
                          var_power = 2,
                          n_gh = 20,
                          gh = NULL) {

  X <- as.matrix(X)
  Z <- as.matrix(Z)

  pr <- safe_plogis(as.numeric(Z %*% al))
  dpr <- pr * (1 - pr)

  # I_lnorm_grad() now returns both probability and gradient.
  # Avoid a separate evaluation of I_lnorm().
  gr <- I_lnorm_grad(
    z = y, X = X, be = be, phi = phi, d = 1,
    var_power = var_power,
    n_gh = n_gh,
    gh = gh
  )

  I0 <- as.numeric(gr$prob[, 1L])
  db <- gr$d_be[, , 1L]
  dp <- as.numeric(gr$d_phi[, 1L])

  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)

  d_be <- matrix(0, n, p)
  d_phi <- numeric(n)
  d_al <- matrix(0, n, q)

  idx0 <- (y == 0)
  idxp <- !idx0

  ## y > 0: P = pr * I
  d_be[idxp, ] <- pr[idxp] * db[idxp, , drop = FALSE]
  d_phi[idxp] <- pr[idxp] * dp[idxp]
  d_al[idxp, ] <- Z[idxp, , drop = FALSE] *
    as.numeric(dpr[idxp] * I0[idxp])

  ## y = 0: P = 1 - pr + pr * I = 1 - pr*(1-I)
  d_be[idx0, ] <- pr[idx0] * db[idx0, , drop = FALSE]
  d_phi[idx0] <- pr[idx0] * dp[idx0]
  d_al[idx0, ] <- Z[idx0, , drop = FALSE] *
    as.numeric(dpr[idx0] * (I0[idx0] - 1))

  # Full zero-inflated probability
  prob <- pr * I0
  prob[idx0] <- (1 - pr[idx0]) + pr[idx0] * I0[idx0]

  list(
    prob = prob,
    d_be = d_be,
    d_phi = d_phi,
    d_al = d_al
  )
}

# del <- 1e-6
# ind <- 1:5
# re <- ZI_lnorm(y = Y[ind], X = X0[ind,], Z = Z0[ind,], be = be1, al = al1, phi = phi1)
# re1 <- ZI_lnorm_grad(y = Y[ind], X = X0[ind,], Z = Z0[ind,], be = be1, al = al1, phi = phi1)
# 
# (I_lnorm(z = Y[ind], X = X0[ind,], be = be1 + c(0,0,0), phi = phi1+del, d = 1)-re)/del
# re1$d_phi
# 

