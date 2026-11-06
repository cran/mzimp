# ============================================================
# Gauss-Laguerre quadrature rule
# ============================================================

.loglogistic_laguerre_rule <- function(
    lag_nodes = NULL,
    lag_weights = NULL,
    n_lag = 32L
) {

  # ----------------------------------------------------------
  # Use user-supplied quadrature rule
  # ----------------------------------------------------------

  if (!is.null(lag_nodes) || !is.null(lag_weights)) {

    if (is.null(lag_nodes) || is.null(lag_weights)) {
      stop(
        "lag_nodes and lag_weights must either both be supplied ",
        "or both be NULL."
      )
    }

    lag_nodes <- as.numeric(lag_nodes)
    lag_weights <- as.numeric(lag_weights)

    if (length(lag_nodes) == 0L ||
        length(lag_nodes) != length(lag_weights)) {
      stop(
        "lag_nodes and lag_weights must have the same ",
        "positive length."
      )
    }

    if (any(!is.finite(lag_nodes)) ||
        any(lag_nodes < 0)) {
      stop(
        "lag_nodes must contain finite nonnegative values."
      )
    }

    if (any(!is.finite(lag_weights)) ||
        any(lag_weights <= 0)) {
      stop(
        "lag_weights must contain finite positive values."
      )
    }

    return(list(
      nodes = lag_nodes,
      weights = lag_weights
    ))
  }

  # ----------------------------------------------------------
  # Generate Gauss-Laguerre quadrature rule
  # ----------------------------------------------------------

  if (!is.numeric(n_lag) ||
      length(n_lag) != 1L ||
      !is.finite(n_lag) ||
      n_lag <= 0 ||
      abs(n_lag - round(n_lag)) > 1e-12) {
    stop("n_lag must be one positive integer.")
  }

  n_lag <- as.integer(round(n_lag))

  if (!requireNamespace("statmod", quietly = TRUE)) {
    stop(
      "Package 'statmod' is required to generate ",
      "Gauss-Laguerre nodes and weights. Install it using ",
      "install.packages('statmod'), or supply lag_nodes ",
      "and lag_weights directly."
    )
  }

  rule <- statmod::gauss.quad(
    n = n_lag,
    kind = "laguerre"
  )

  list(
    nodes = as.numeric(rule$nodes),
    weights = as.numeric(rule$weights)
  )
}


# ============================================================
# I_loglogistic
# Log-logistic mixed Poisson interval probabilities
# ============================================================

I_loglogistic <- function(
    z,
    X,
    be,
    phi,
    d,
    var_power = 2,
    lag_nodes = NULL,
    lag_weights = NULL,
    n_lag = 32L,
    log_out = FALSE,
    eps_u = 1e-15
) {

  X <- as.matrix(X)
  be <- as.numeric(be)

  # ----------------------------------------------------------
  # Input validation
  # ----------------------------------------------------------

  if (nrow(X) != length(z)) {
    stop("nrow(X) must equal length(z).")
  }

  if (length(be) != ncol(X)) {
    stop("length(be) must equal ncol(X).")
  }

  if (!is.numeric(z) ||
      any(!is.finite(z)) ||
      any(z < 0) ||
      any(abs(z - round(z)) > 1e-12)) {
    stop(
      "z must be a numeric vector of nonnegative integers."
    )
  }

  z <- as.integer(round(z))

  if (!is.numeric(phi) ||
      length(phi) != 1L ||
      !is.finite(phi) ||
      phi <= 0) {
    stop("phi must be a positive scalar.")
  }

  if (!is.numeric(d) ||
      length(d) != 1L ||
      !is.finite(d) ||
      d <= 0 ||
      abs(d - round(d)) > 1e-12) {
    stop("d must be one positive integer.")
  }

  d <- as.integer(round(d))

  if (!is.numeric(var_power) ||
      length(var_power) != 1L ||
      !is.finite(var_power) ||
      abs(var_power - round(var_power)) > 1e-12 ||
      !(as.integer(round(var_power)) %in% c(1L, 2L, 3L))) {
    stop("var_power must be one of 1, 2, or 3.")
  }

  var_power <- as.integer(round(var_power))

  if (!is.logical(log_out) ||
      length(log_out) != 1L ||
      is.na(log_out)) {
    stop("log_out must be TRUE or FALSE.")
  }

  if (!is.numeric(eps_u) ||
      length(eps_u) != 1L ||
      !is.finite(eps_u) ||
      eps_u <= 0 ||
      eps_u >= 0.5) {
    stop("eps_u must lie in (0, 0.5).")
  }

  # ----------------------------------------------------------
  # Conditional means
  # ----------------------------------------------------------

  eta <- as.vector(X %*% be)
  mu <- exp(eta)

  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop(
      "mu contains non-finite or non-positive values."
    )
  }

  # ----------------------------------------------------------
  # Gauss-Laguerre rule
  # ----------------------------------------------------------

  lag_rule <- .loglogistic_laguerre_rule(
    lag_nodes = lag_nodes,
    lag_weights = lag_weights,
    n_lag = n_lag
  )

  # ----------------------------------------------------------
  # Evaluate interval probabilities
  # ----------------------------------------------------------

  if (!exists(
    "I_loglogistic_mu_phi_cpp",
    mode = "function"
  )) {
    stop(
      "The C++ backend requires a compiled function named ",
      "'I_loglogistic_mu_phi_cpp'."
    )
  }

  out <- I_loglogistic_mu_phi_cpp(
    z = z,
    mu = mu,
    phi = phi,
    d = d,
    var_power = var_power,
    lag_nodes = lag_rule$nodes,
    lag_weights = lag_rule$weights,
    log_out = log_out,
    eps_u = eps_u
  )

  out <- as.matrix(out)

  rownames(out) <- paste0(
    "obs_",
    seq_along(z)
  )

  colnames(out) <- paste0(
    "I_",
    seq_len(d)
  )

  out
}


# ============================================================
# I_loglogistic_grad
# Gradient with respect to beta and phi
# ============================================================

I_loglogistic_grad <- function(
    z,
    X,
    be,
    phi,
    d,
    var_power = 2,
    lag_nodes = NULL,
    lag_weights = NULL,
    n_lag = 32L,
    eps_u = 1e-15
) {

  X <- as.matrix(X)
  be <- as.numeric(be)

  # ----------------------------------------------------------
  # Input validation
  # ----------------------------------------------------------

  if (nrow(X) != length(z)) {
    stop("nrow(X) must equal length(z).")
  }

  if (length(be) != ncol(X)) {
    stop("length(be) must equal ncol(X).")
  }

  if (!is.numeric(z) ||
      any(!is.finite(z)) ||
      any(z < 0) ||
      any(abs(z - round(z)) > 1e-12)) {
    stop(
      "z must be a numeric vector of nonnegative integers."
    )
  }

  z <- as.integer(round(z))

  if (!is.numeric(phi) ||
      length(phi) != 1L ||
      !is.finite(phi) ||
      phi <= 0) {
    stop("phi must be a positive scalar.")
  }

  if (!is.numeric(d) ||
      length(d) != 1L ||
      !is.finite(d) ||
      d <= 0 ||
      abs(d - round(d)) > 1e-12) {
    stop("d must be one positive integer.")
  }

  d <- as.integer(round(d))

  if (!is.numeric(var_power) ||
      length(var_power) != 1L ||
      !is.finite(var_power) ||
      abs(var_power - round(var_power)) > 1e-12 ||
      !(as.integer(round(var_power)) %in% c(1L, 2L, 3L))) {
    stop("var_power must be one of 1, 2, or 3.")
  }

  var_power <- as.integer(round(var_power))

  if (!is.numeric(eps_u) ||
      length(eps_u) != 1L ||
      !is.finite(eps_u) ||
      eps_u <= 0 ||
      eps_u >= 0.5) {
    stop("eps_u must lie in (0, 0.5).")
  }

  # ----------------------------------------------------------
  # Conditional means
  # ----------------------------------------------------------

  eta <- as.vector(X %*% be)
  mu <- exp(eta)

  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop(
      "mu contains non-finite or non-positive values."
    )
  }

  n <- length(z)
  p <- ncol(X)

  # ----------------------------------------------------------
  # Gauss-Laguerre rule
  # ----------------------------------------------------------

  lag_rule <- .loglogistic_laguerre_rule(
    lag_nodes = lag_nodes,
    lag_weights = lag_weights,
    n_lag = n_lag
  )

  # ----------------------------------------------------------
  # Gradient with respect to mu and phi
  # ----------------------------------------------------------

  if (!exists(
    "I_loglogistic_mu_phi_grad_cpp",
    mode = "function"
  )) {
    stop(
      "The C++ backend requires a compiled function named ",
      "'I_loglogistic_mu_phi_grad_cpp'."
    )
  }

  g <- I_loglogistic_mu_phi_grad_cpp(
    z = z,
    mu = mu,
    phi = phi,
    d = d,
    var_power = var_power,
    lag_nodes = lag_rule$nodes,
    lag_weights = lag_rule$weights,
    eps_u = eps_u
  )

  prob <- if (!is.null(g$prob)) {
    as.matrix(g$prob)
  } else if (!is.null(g$I)) {
    as.matrix(g$I)
  } else {
    as.matrix(g$value)
  }
  d_mu <- as.matrix(g$d_mu)
  d_phi <- as.matrix(g$d_phi)

  expected_dimension <- c(n, d)

  if (!identical(dim(prob), expected_dimension)) {
    stop(
      "The dimension of probability returned by the backend is invalid."
    )
  }

  if (!identical(dim(d_mu), expected_dimension)) {
    stop(
      "The dimension of g$d_mu returned by the backend is invalid."
    )
  }

  if (!identical(dim(d_phi), expected_dimension)) {
    stop(
      "The dimension of g$d_phi returned by the backend is invalid."
    )
  }

  # ----------------------------------------------------------
  # Chain rule:
  #
  # mu_i = exp(x_i' beta)
  #
  # d I_i^(m) / d beta_j
  # = d I_i^(m) / d mu_i * mu_i * x_ij
  # ----------------------------------------------------------

  d_be <- array(
    0,
    dim = c(n, p, d)
  )

  for (m in seq_len(d)) {
    d_be[, , m] <-
      X * as.numeric(mu * d_mu[, m])
  }

  # ----------------------------------------------------------
  # Names
  # ----------------------------------------------------------

  observation_names <- paste0(
    "obs_",
    seq_len(n)
  )

  interval_names <- paste0(
    "I_",
    seq_len(d)
  )

  coefficient_names <- colnames(X)

  if (is.null(coefficient_names)) {
    coefficient_names <- paste0(
      "be_",
      seq_len(p)
    )
  }

  rownames(prob) <- observation_names
  rownames(d_mu) <- observation_names
  rownames(d_phi) <- observation_names

  colnames(prob) <- interval_names
  colnames(d_mu) <- interval_names
  colnames(d_phi) <- interval_names

  dimnames(d_be) <- list(
    observation_names,
    coefficient_names,
    interval_names
  )

  list(
    prob = prob,
    d_mu = d_mu,
    d_phi = d_phi,
    d_be = d_be
  )
}


# ============================================================
# ZI_loglogistic
# Zero-inflated log-logistic mixed Poisson
#
# Interpretation:
#
#   pr_i = P(observation belongs to count component)
#
# Thus:
#
#   P(Y_i = y_i)
#   = (1 - pr_i) 1(y_i = 0)
#     + pr_i P_MP(Y_i = y_i).
# ============================================================

ZI_loglogistic <- function(
    y,
    X,
    be,
    Z,
    al,
    phi,
    var_power = 2,
    lag_nodes = NULL,
    lag_weights = NULL,
    n_lag = 32L,
    eps_u = 1e-15
) {

  X <- as.matrix(X)
  Z <- as.matrix(Z)

  be <- as.numeric(be)
  al <- as.numeric(al)

  # ----------------------------------------------------------
  # Input validation
  # ----------------------------------------------------------

  if (!is.numeric(y) ||
      any(!is.finite(y)) ||
      any(y < 0) ||
      any(abs(y - round(y)) > 1e-12)) {
    stop(
      "y must be a numeric vector of nonnegative integers."
    )
  }

  y <- as.integer(round(y))

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

  if (!is.numeric(phi) ||
      length(phi) != 1L ||
      !is.finite(phi) ||
      phi <= 0) {
    stop("phi must be a positive scalar.")
  }

  if (!is.numeric(var_power) ||
      length(var_power) != 1L ||
      !is.finite(var_power) ||
      abs(var_power - round(var_power)) > 1e-12 ||
      !(as.integer(round(var_power)) %in% c(1L, 2L, 3L))) {
    stop("var_power must be one of 1, 2, or 3.")
  }

  var_power <- as.integer(round(var_power))

  # ----------------------------------------------------------
  # Probability of belonging to count component
  # ----------------------------------------------------------

  eta_pr <- as.numeric(Z %*% al)
  pr <- plogis(eta_pr)

  # ----------------------------------------------------------
  # Log-logistic mixed Poisson marginal probability
  # ----------------------------------------------------------

  base <- I_loglogistic(
    z = y,
    X = X,
    be = be,
    phi = phi,
    d = 1L,
    var_power = var_power,
    lag_nodes = lag_nodes,
    lag_weights = lag_weights,
    n_lag = n_lag,
    log_out = FALSE,
    eps_u = eps_u
  )

  base <- as.numeric(base)

  # ----------------------------------------------------------
  # Zero-inflated probability
  # ----------------------------------------------------------

  out <-
    (1 - pr) * as.numeric(y == 0L) +
    pr * base

  out
}


# ============================================================
# ZI_loglogistic_grad
# Gradient of zero-inflated log-logistic mixed Poisson
# ============================================================

ZI_loglogistic_grad <- function(
    y,
    X,
    be,
    Z,
    al,
    phi,
    var_power = 2,
    lag_nodes = NULL,
    lag_weights = NULL,
    n_lag = 32L,
    eps_u = 1e-15
) {

  X <- as.matrix(X)
  Z <- as.matrix(Z)

  be <- as.numeric(be)
  al <- as.numeric(al)

  # ----------------------------------------------------------
  # Input validation
  # ----------------------------------------------------------

  if (!is.numeric(y) ||
      any(!is.finite(y)) ||
      any(y < 0) ||
      any(abs(y - round(y)) > 1e-12)) {
    stop(
      "y must be a numeric vector of nonnegative integers."
    )
  }

  y <- as.integer(round(y))

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

  if (!is.numeric(phi) ||
      length(phi) != 1L ||
      !is.finite(phi) ||
      phi <= 0) {
    stop("phi must be a positive scalar.")
  }

  if (!is.numeric(var_power) ||
      length(var_power) != 1L ||
      !is.finite(var_power) ||
      abs(var_power - round(var_power)) > 1e-12 ||
      !(as.integer(round(var_power)) %in% c(1L, 2L, 3L))) {
    stop("var_power must be one of 1, 2, or 3.")
  }

  var_power <- as.integer(round(var_power))

  # ----------------------------------------------------------
  # Conditional mean and count-component probability
  # ----------------------------------------------------------

  mu <- as.numeric(
    exp(X %*% be)
  )

  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop(
      "mu contains non-finite or non-positive values."
    )
  }

  eta_pr <- as.numeric(Z %*% al)
  pr <- plogis(eta_pr)

  # ----------------------------------------------------------
  # Log-logistic mixed Poisson probability and derivatives
  # ----------------------------------------------------------

  gr <- I_loglogistic_grad(
    z = y,
    X = X,
    be = be,
    phi = phi,
    d = 1L,
    var_power = var_power,
    lag_nodes = lag_nodes,
    lag_weights = lag_weights,
    n_lag = n_lag,
    eps_u = eps_u
  )

  base <- as.numeric(gr$prob[, 1L])

  db <- gr$d_be
  if (length(dim(db)) == 3L) {
    db <- db[, , 1L]
  }

  base_d_phi <- as.numeric(gr$d_phi[, 1L])

  # ----------------------------------------------------------
  # Zero-inflated probability and derivatives
  # ----------------------------------------------------------

  is_zero <- as.numeric(y == 0L)

  prob <-
    (1 - pr) * is_zero +
    pr * base

  d_be <- db * pr
  d_phi <- pr * base_d_phi

  d_pr <- base - is_zero

  d_al <-
    Z * as.numeric(
      d_pr * pr * (1 - pr)
    )

  list(
    prob = prob,
    d_al = d_al,
    d_be = d_be,
    d_phi = d_phi
  )
}
