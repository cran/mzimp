# ============================================================
# Gauss-Laguerre quadrature rule for symmetric beta-prime
# ============================================================

.betaprime_laguerre_rule <- function(
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

    return(
      list(
        nodes = lag_nodes,
        weights = lag_weights
      )
    )
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
# I_betaprime
# symmetric beta-prime mixed Poisson interval probabilities
# ============================================================

I_betaprime <- function(
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

  if (!is.numeric(var_power) ||
      length(var_power) != 1L ||
      !is.finite(var_power)) {
    stop("var_power must be one finite numeric scalar.")
  }

  if (!is.numeric(d) ||
      length(d) != 1L ||
      !is.finite(d) ||
      d <= 0 ||
      abs(d - round(d)) > 1e-12) {
    stop("d must be one positive integer.")
  }

  d <- as.integer(round(d))
  
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

  if (any(!is.finite(mu)) ||
      any(mu <= 0)) {
    stop(
      "mu contains non-finite or non-positive values."
    )
  }

  # ----------------------------------------------------------
  # Gauss-Laguerre rule
  # ----------------------------------------------------------

  lag_rule <- .betaprime_laguerre_rule(
    lag_nodes = lag_nodes,
    lag_weights = lag_weights,
    n_lag = n_lag
  )

  # ----------------------------------------------------------
  # Evaluate interval probabilities
  # ----------------------------------------------------------

  if (!exists(
      "I_betaprime_mu_phi_cpp",
      mode = "function",
      inherits = TRUE
    )) {
      stop(
        "Function 'I_betaprime_mu_phi_cpp' was not found. ",
        "Compile the beta-prime C++ code first."
      )
    }

    out <- I_betaprime_mu_phi_cpp(
      z = z,
      mu = mu,
      phi = phi,
      d = d,
      lag_nodes = lag_rule$nodes,
      lag_weights = lag_rule$weights,
      var_power = var_power,

      log_out = log_out,
      eps_u = eps_u
    )


  out <- as.matrix(out)

  expected_dimension <- c(
    length(z),
    d
  )

  if (!identical(dim(out), expected_dimension)) {
    stop(
      "The backend returned an object with an invalid dimension."
    )
  }

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
# I_betaprime_grad
# Gradient with respect to regression coefficients and phi
# ============================================================

I_betaprime_grad <- function(
    z,
    X,
    be,
    phi,
    d,
    var_power = 2,
    lag_nodes = NULL,
    lag_weights = NULL,
    n_lag = 24L,
    eps_u = 1e-15,
    fd_rel_step = 1e-5
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

  if (!is.numeric(var_power) ||
      length(var_power) != 1L ||
      !is.finite(var_power)) {
    stop("var_power must be one finite numeric scalar.")
  }

  if (!is.numeric(d) ||
      length(d) != 1L ||
      !is.finite(d) ||
      d <= 0 ||
      abs(d - round(d)) > 1e-12) {
    stop("d must be one positive integer.")
  }

  d <- as.integer(round(d))

  if (!is.numeric(eps_u) ||
      length(eps_u) != 1L ||
      !is.finite(eps_u) ||
      eps_u <= 0 ||
      eps_u >= 0.5) {
    stop("eps_u must lie in (0, 0.5).")
  }

  if (!is.numeric(fd_rel_step) ||
      length(fd_rel_step) != 1L ||
      !is.finite(fd_rel_step) ||
      fd_rel_step <= 0) {
    stop("fd_rel_step must be one positive finite scalar.")
  }
# ----------------------------------------------------------
  # Conditional means
  # ----------------------------------------------------------

  eta <- as.vector(X %*% be)
  mu <- exp(eta)

  if (any(!is.finite(mu)) ||
      any(mu <= 0)) {
    stop(
      "mu contains non-finite or non-positive values."
    )
  }

  n <- length(z)
  p <- ncol(X)

  # ----------------------------------------------------------
  # Gauss-Laguerre rule
  # ----------------------------------------------------------

  lag_rule <- .betaprime_laguerre_rule(
    lag_nodes = lag_nodes,
    lag_weights = lag_weights,
    n_lag = n_lag
  )

  # ----------------------------------------------------------
  # Gradient with respect to mu and phi
  # ----------------------------------------------------------

  if (!exists(
      "I_betaprime_mu_phi_grad_cpp",
      mode = "function",
      inherits = TRUE
    )) {
      stop(
        "Function 'I_betaprime_mu_phi_grad_cpp' was not found. ",
        "Compile the beta-prime C++ code first."
      )
    }

    g <- I_betaprime_mu_phi_grad_cpp(
      z = z,
      mu = mu,
      phi = phi,
      d = d,
      lag_nodes = lag_rule$nodes,
      lag_weights = lag_rule$weights,
      var_power = var_power,
      eps_u = eps_u,
      fd_rel_step = fd_rel_step
    )


  # ----------------------------------------------------------
  # Extract results
  # ----------------------------------------------------------

  if (is.null(g$I) &&
      !is.null(g$value)) {
    g$I <- g$value
  }

  if (is.null(g$I)) {
    stop(
      "The gradient backend must return an element named 'I' ",
      "or 'value'."
    )
  }

  if (is.null(g$d_mu)) {
    stop(
      "The gradient backend must return an element named 'd_mu'."
    )
  }

  if (is.null(g$d_phi)) {
    stop(
      "The gradient backend must return an element named 'd_phi'."
    )
  }

  I_value <- as.matrix(g$I)
  d_mu <- as.matrix(g$d_mu)
  d_phi <- as.matrix(g$d_phi)

  expected_dimension <- c(n, d)

  if (!identical(dim(I_value), expected_dimension)) {
    stop(
      "The dimension of g$I returned by the backend is invalid."
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
  # Chain rule
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
      X *
      as.numeric(
        mu * d_mu[, m]
      )
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

  rownames(I_value) <- observation_names
  rownames(d_mu) <- observation_names
  rownames(d_phi) <- observation_names

  colnames(I_value) <- interval_names
  colnames(d_mu) <- interval_names
  colnames(d_phi) <- interval_names

  dimnames(d_be) <- list(
    observation_names,
    coefficient_names,
    interval_names
  )

  # ----------------------------------------------------------
  # Aggregated beta gradients
  # ----------------------------------------------------------

  
  list(
  	I = I_value,
  	value = I_value,
  	d_mu = d_mu,
  	d_phi = d_phi,
  	d_be = d_be,
  	mu = mu
  )
}

# ============================================================
# ZI_betaprime
# Zero-inflated symmetric beta-prime mixed Poisson
#
# Interpretation:
#
#   pr_i = P(observation belongs to the count component)
#
# Therefore:
#
#   P(Y_i = y_i)
#   = (1 - pr_i) 1(y_i = 0)
#     + pr_i P_MP(Y_i = y_i).
# ============================================================

ZI_betaprime <- function(
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

  # ----------------------------------------------------------
  # Probability of belonging to count component
  # ----------------------------------------------------------

  eta_pr <- as.numeric(Z %*% al)
  pr <- plogis(eta_pr)

  # ----------------------------------------------------------
  # beta-prime mixed Poisson marginal probability
  # ----------------------------------------------------------

  base <- I_betaprime(
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
    (1 - pr) *
    as.numeric(y == 0L) +
    pr *
    base

  out
}

# ============================================================
# ZI_betaprime_grad
# Gradient of zero-inflated symmetric beta-prime mixed Poisson
# ============================================================

ZI_betaprime_grad <- function(
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
    eps_u = 1e-15,
    fd_rel_step = 1e-5
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

  if (!is.numeric(fd_rel_step) ||
      length(fd_rel_step) != 1L ||
      !is.finite(fd_rel_step) ||
      fd_rel_step <= 0) {
    stop("fd_rel_step must be one positive finite scalar.")
  }

  # ----------------------------------------------------------
  # Conditional mean
  # ----------------------------------------------------------

  mu <- as.numeric(
    exp(X %*% be)
  )

  if (any(!is.finite(mu)) ||
      any(mu <= 0)) {
    stop(
      "mu contains non-finite or non-positive values."
    )
  }

  # ----------------------------------------------------------
  # Probability of belonging to count component
  # ----------------------------------------------------------

  eta_pr <- as.numeric(Z %*% al)
  pr <- plogis(eta_pr)

  # ----------------------------------------------------------
  # beta-prime mixed Poisson probability and derivatives
  # ----------------------------------------------------------

  gr <- I_betaprime_grad(
    z = y,
    X = X,
    be = be,
    phi = phi,
    d = 1L,
    var_power = var_power,
    lag_nodes = lag_nodes,
    lag_weights = lag_weights,
    n_lag = n_lag,
    eps_u = eps_u,
    fd_rel_step = fd_rel_step
  )

  base <- as.numeric(gr$I)
  d_base_mu <- as.numeric(gr$d_mu)
  d_base_phi <- as.numeric(gr$d_phi)

  # ----------------------------------------------------------
  # Derivatives of
  #
  # f_i = (1-pr_i) 1(y_i=0) + pr_i base_i
  # ----------------------------------------------------------

  is_zero <- as.numeric(
    y == 0L
  )

  d_mu <-
    pr *
    d_base_mu

  d_phi <-
    pr *
    d_base_phi

  d_pr <-
    base -
    is_zero

  # ----------------------------------------------------------
  # Chain rule for regression coefficients
  # ----------------------------------------------------------

  d_be <-
    X *
    as.numeric(
      mu * d_mu
    )

  d_al <-
    Z *
    as.numeric(
      d_pr *
      pr *
      (1 - pr)
    )

  value <-
    (1 - pr) *
    is_zero +
    pr *
    base

  list(
    value = value,
    probability = value,
    d_al = d_al,
    d_be = d_be,
    d_phi = d_phi,
    d_mu = d_mu,
    d_pr = d_pr,
    d_base_mu = d_base_mu,
    d_base_phi = d_base_phi,
    base = base,
    pr = pr,
    mu = mu
  )
}