# ============================================================
# Multivariate ZIMP simulation
#
# This simulator reuses simu_zipm() margin-by-margin so that the
# marginal parameterizations are exactly the same as in the
# univariate simulator (including native/default variance powers).
#
# Dependence:
#   Count/mixing layer: Gaussian factor copula with loading psi_N (k x J)
#   Bernoulli layer:    Gaussian one-factor copula with loading psi_B (k x 1)
#
# Identification:
#   psi_N[j, l] = 0 whenever l > j, with psi_N[j,j] > 0 for j <= J.
#   For the one-factor Bernoulli layer, the first active psi_B loading > 0.
#   Off-diagonal / remaining loadings may have either sign.
# ============================================================


# ------------------------------------------------------------
# Random loading generator for the count/mixing layer
# ------------------------------------------------------------
simu_psi_N <- function(
    k,
    J,
    loading.min = 0.65,
    loading.max = 0.85,
    seed = NULL
) {

    if (!is.null(seed)) set.seed(seed)

    if (length(k) != 1L || length(J) != 1L ||
        !is.finite(k) || !is.finite(J) ||
        k < 1L || J < 1L ||
        k != as.integer(k) || J != as.integer(J) || J > k) {
        stop("k and J must be positive integers with J <= k.")
    }

    if (!is.finite(loading.min) || !is.finite(loading.max) ||
        loading.min <= 0 || loading.max >= 1 ||
        loading.min >= loading.max) {
        stop("Require 0 < loading.min < loading.max < 1.")
    }

    k <- as.integer(k)
    J <- as.integer(J)
    psi <- matrix(0, k, J)

    for (j in seq_len(k)) {

        active <- seq_len(min(j, J))
        qj <- length(active)

        # Random direction. Off-diagonal entries may have either sign,
        # which allows both positive and negative implied correlations.
        v <- runif(qj, 0.8, 1.0)
        signs <- sample(c(-1, 1), qj, replace = TRUE)
        v <- v * signs

        # Identification constraint:
        # for rows j <= J, psi[j,j] is the identifying diagonal
        # and is forced to be strictly positive.
        if (j <= J) {
            v[j] <- abs(v[j])
        }

        v <- v / sqrt(sum(v^2))

        # Keep every row norm below one.
        r <- runif(1, loading.min, loading.max)
        psi[j, active] <- r * v
    }

    colnames(psi) <- paste0("N_factor", seq_len(J))
    rownames(psi) <- paste0("margin", seq_len(k))

    psi
}


# ------------------------------------------------------------
# Random loading generator for the Bernoulli layer
# ------------------------------------------------------------
simu_psi_B <- function(
    k,
    loading.min = 0.65,
    loading.max = 0.85,
    seed = NULL
) {

    if (!is.null(seed)) set.seed(seed)

    if (length(k) != 1L || !is.finite(k) ||
        k < 1L || k != as.integer(k)) {
        stop("k must be a positive integer.")
    }

    if (!is.finite(loading.min) || !is.finite(loading.max) ||
        loading.min <= 0 || loading.max >= 1 ||
        loading.min >= loading.max) {
        stop("Require 0 < loading.min < loading.max < 1.")
    }

    k <- as.integer(k)

    magnitude <- runif(k, loading.min, loading.max)
    signs <- sample(c(-1, 1), k, replace = TRUE)

    # Identification for the one-factor Bernoulli layer:
    # psi_B and -psi_B are observationally equivalent, so fix
    # the first zero-inflated margin loading to be positive.
    signs[1L] <- 1

    psi <- matrix(
        magnitude * signs,
        nrow = k,
        ncol = 1L
    )

    colnames(psi) <- "B_factor1"
    rownames(psi) <- paste0("ZI_margin", seq_len(k))

    psi
}


# ------------------------------------------------------------
# Main multivariate simulator
# ------------------------------------------------------------
simu_multizimp <- function(
    X,
    beta,
    Z = NULL,
    gamma = NULL,
    zero_inflation = NULL,
    mixing.distribution,
    psi_N,
    psi_B = NULL,
    phi = 1,
    seed = NULL,
    return_latent = TRUE
) {

  if (!is.null(seed)) {
    if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed)) {
      stop("seed must be NULL or a finite numeric scalar.")
    }
    set.seed(seed)
  }

  if (!exists("simu_zipm", mode = "function", inherits = TRUE)) {
    stop(
      "simu_multizimp() requires simu_zipm() to be loaded. ",
      "Load/source the univariate simulator first."
    )
  }

  # ----------------------------------------------------------
  # X and beta
  # ----------------------------------------------------------
  if (is.data.frame(X)) X <- as.matrix(X)

  if (!is.matrix(X) || !is.numeric(X) ||
      anyNA(X) || any(!is.finite(X))) {
    stop("X must be a finite numeric matrix.")
  }

  n <- nrow(X)
  k1 <- ncol(X)

  if (n < 1L || k1 < 1L) {
    stop("X must have at least one row and one column.")
  }

  if (!is.numeric(beta) || length(beta) != k1 ||
      anyNA(beta) || any(!is.finite(beta))) {
    stop("beta must be a finite numeric vector of length ncol(X).")
  }
  beta <- as.numeric(beta)

  # ----------------------------------------------------------
  # psi_N determines k and J
  # ----------------------------------------------------------
  if (is.data.frame(psi_N)) psi_N <- as.matrix(psi_N)

  if (!is.matrix(psi_N) || !is.numeric(psi_N) ||
      anyNA(psi_N) || any(!is.finite(psi_N))) {
    stop("psi_N must be a finite numeric matrix.")
  }

  k <- nrow(psi_N)
  J <- ncol(psi_N)

  if (k < 1L || J < 1L) {
    stop("psi_N must be a non-empty k x J matrix.")
  }
  if (J > k) {
    stop("ncol(psi_N) = J must not exceed k.")
  }

  for (j in seq_len(k)) {
    if (j < J) {
      upper_index <- seq.int(j + 1L, J)
      if (any(abs(psi_N[j, upper_index]) > 1e-12)) {
        stop(
          "psi_N must satisfy psi_N[j,l] = 0 for l > j. ",
          "Violation found in row ", j, "."
        )
      }
    }
  }

  # Canonical factor identification used by both simulation and fitting:
  # lower triangular leading J x J block with strictly positive diagonal.
  if (any(diag(psi_N[seq_len(J), seq_len(J), drop = FALSE]) <= 0)) {
    stop(
      "psi_N must have a strictly positive identifying diagonal: ",
      "psi_N[j,j] > 0 for j = 1,...,J."
    )
  }

  psi_N_norm2 <- rowSums(psi_N^2)
  if (any(psi_N_norm2 >= 1)) {
    stop("Each row of psi_N must have squared norm strictly smaller than 1.")
  }

  # ----------------------------------------------------------
  # Zero-inflation specification
  #
  # Z = NULL:
  #   no margin is zero-inflated; gamma and psi_B are unused.
  #
  # Z != NULL:
  #   zero_inflation must be logical length k.
  #   psi_B contains loadings ONLY for TRUE margins, in their
  #   natural margin order.
  # ----------------------------------------------------------
  if (is.null(Z)) {

    zero_inflation <- rep(FALSE, k)
    zi_index <- integer(0)
    n_zi <- 0L

    gamma <- NULL
    psi_B <- NULL
    prob <- rep(1, n)

  } else {

    if (is.data.frame(Z)) Z <- as.matrix(Z)

    if (!is.matrix(Z) || !is.numeric(Z) ||
        anyNA(Z) || any(!is.finite(Z))) {
      stop("Z must be NULL or a finite numeric matrix.")
    }

    if (nrow(Z) != n) {
      stop("X and Z must have the same number of rows.")
    }

    k2 <- ncol(Z)
    if (k2 < 1L) {
      stop("Z must have at least one column.")
    }

    if (!is.numeric(gamma) || length(gamma) != k2 ||
        anyNA(gamma) || any(!is.finite(gamma))) {
      stop("gamma must be a finite numeric vector of length ncol(Z).")
    }
    gamma <- as.numeric(gamma)

    if (is.null(zero_inflation) ||
        !is.logical(zero_inflation) ||
        length(zero_inflation) != k ||
        anyNA(zero_inflation)) {
      stop(
        "When Z is not NULL, zero_inflation must be a non-missing ",
        "logical vector of length k."
      )
    }

    zi_index <- which(zero_inflation)
    n_zi <- length(zi_index)

    eta_B <- drop(Z %*% gamma)
    prob <- plogis(eta_B)

    if (any(!is.finite(prob))) {
      stop(
        "plogis(Z %*% gamma) contains invalid values. ",
        "Check Z and gamma."
      )
    }

    if (n_zi == 0L) {
      # Allowed: Z supplied, but all margins are PM.
      if (!is.null(psi_B) && length(psi_B) != 0L) {
        stop(
          "sum(zero_inflation) = 0, so psi_B must be NULL or have length 0."
        )
      }
      psi_B <- NULL

    } else {
      if (is.data.frame(psi_B)) psi_B <- as.matrix(psi_B)

      # Accept either a vector or an n_zi x 1 matrix.
      if (is.matrix(psi_B)) {
        if (!is.numeric(psi_B) || ncol(psi_B) != 1L ||
            nrow(psi_B) != n_zi ||
            anyNA(psi_B) || any(!is.finite(psi_B))) {
          stop(
            "psi_B must be a finite numeric vector of length ",
            "sum(zero_inflation), or a sum(zero_inflation) x 1 matrix."
          )
        }
        psi_B_vector <- as.numeric(psi_B[, 1L])
      } else {
        if (!is.numeric(psi_B) || length(psi_B) != n_zi ||
            anyNA(psi_B) || any(!is.finite(psi_B))) {
          stop(
            "psi_B must have length sum(zero_inflation) = ", n_zi, "."
          )
        }
        psi_B_vector <- as.numeric(psi_B)
      }

      if (any(abs(psi_B_vector) >= 1)) {
        stop("All entries of psi_B must satisfy abs(psi_B) < 1.")
      }

      if (psi_B_vector[1L] <= 0) {
        stop(
          "For one-factor identification, the first active loading in ",
          "psi_B must be strictly positive."
        )
      }

      psi_B <- matrix(psi_B_vector, ncol = 1L)
    }
  }

  # ----------------------------------------------------------
  # Mixing distributions and phi
  # ----------------------------------------------------------
  allowed_distribution <- c(
    "gamma", "lnorm", "pareto", "weibull",
    "ig", "loglogistic", "betaprime"
  )

  if (!is.character(mixing.distribution) ||
      length(mixing.distribution) != k) {
    stop("mixing.distribution must be a character vector of length k.")
  }

  mixing.distribution <- tolower(mixing.distribution)
  bad_distribution <- !(mixing.distribution %in% allowed_distribution)

  if (any(bad_distribution)) {
    stop(
      "Unsupported mixing.distribution: ",
      paste(unique(mixing.distribution[bad_distribution]), collapse = ", "),
      "."
    )
  }

  if (length(phi) == 1L) phi <- rep(phi, k)

  if (!is.numeric(phi) || length(phi) != k ||
      anyNA(phi) || any(!is.finite(phi)) || any(phi <= 0)) {
    stop(
      "phi must be a positive scalar or a positive numeric vector ",
      "of length k."
    )
  }
  phi <- as.numeric(phi)

  # ----------------------------------------------------------
  # Marginal count mean
  # ----------------------------------------------------------
  eta_N <- drop(X %*% beta)
  mu <- exp(eta_N)

  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop("exp(X %*% beta) contains invalid values. Check X and beta.")
  }

  # ----------------------------------------------------------
  # Gaussian factor copula for mixing variables
  # ----------------------------------------------------------
  F_N <- matrix(rnorm(n * J), nrow = n, ncol = J)
  epsilon_N <- matrix(rnorm(n * k), nrow = n, ncol = k)

  G_N <- F_N %*% t(psi_N)

  residual_sd_N <- sqrt(pmax(0, 1 - psi_N_norm2))
  G_N <- G_N + sweep(
    epsilon_N, MARGIN = 2L,
    STATS = residual_sd_N, FUN = "*"
  )

  U_N <- pnorm(G_N)
  eps <- 1e-12
  U_N[] <- pmin(1 - eps, pmax(eps, as.vector(U_N)))

  # ----------------------------------------------------------
  # Bernoulli layer
  #
  # B is n x k for convenience. PM margins are identically 1.
  # Only zero-inflated margins participate in the B factor.
  # ----------------------------------------------------------
  B <- matrix(1L, nrow = n, ncol = k)
  U_B <- NULL
  F_B <- NULL
  correlation_B <- NULL

  if (n_zi > 0L) {

    F_B <- rnorm(n)

    epsilon_B <- matrix(
      rnorm(n * n_zi),
      nrow = n,
      ncol = n_zi
    )

    psi_B_vector <- as.numeric(psi_B[, 1L])

    G_B <-
      outer(F_B, psi_B_vector) +
      sweep(
        epsilon_B,
        MARGIN = 2L,
        STATS = sqrt(1 - psi_B_vector^2),
        FUN = "*"
      )

    U_B <- pnorm(G_B)

    # prob_i is P(B_ij = 1), i.e. Poisson-component probability.
    B[, zi_index] <- matrix(
      as.integer(U_B <= prob),
      nrow = n,
      ncol = n_zi
    )

    correlation_B <- psi_B %*% t(psi_B)
    diag(correlation_B) <- 1
  }

  # ----------------------------------------------------------
  # Simulate each mixed-Poisson margin using simu_zipm()
  # ----------------------------------------------------------
  N_count <- matrix(0L, nrow = n, ncol = k)
  Lambda <- matrix(NA_real_, nrow = n, ncol = k)
  theoretical_variance <- matrix(NA_real_, nrow = n, ncol = k)
  var_power <- rep(NA_integer_, k)

  for (j in seq_len(k)) {

    margin_j <- simu_zipm(
      n = n,
      mu = mu,
      phi = phi[j],
      latent.distribution = mixing.distribution[j],
      use_custom_var_power = FALSE,
      zero_inflation = FALSE,
      prob = NULL,
      seed = NULL,
      uniform.import = TRUE,
      U.import = U_N[, j],
      return_latent = TRUE
    )

    N_count[, j] <- as.integer(margin_j$Y)
    Lambda[, j] <- as.numeric(margin_j$Lambda)
    theoretical_variance[, j] <- as.numeric(
      margin_j$theoretical_variance
    )

    vp_j <- attr(margin_j, "var_power")
    if (length(vp_j) == 1L) {
      var_power[j] <- as.integer(vp_j)
    }
  }

  # ----------------------------------------------------------
  # Observed response
  # ----------------------------------------------------------
  Y <- B * N_count
  storage.mode(Y) <- "integer"

  response_names <- paste0("Y", seq_len(k))
  colnames(Y) <- response_names
  colnames(N_count) <- response_names
  colnames(Lambda) <- response_names
  colnames(B) <- response_names
  colnames(U_N) <- response_names

  rownames(psi_N) <- response_names

  if (n_zi > 0L) {
    zi_names <- response_names[zi_index]
    colnames(U_B) <- zi_names
    rownames(psi_B) <- zi_names
    dimnames(correlation_B) <- list(zi_names, zi_names)
  }

  correlation_N <- psi_N %*% t(psi_N)
  diag(correlation_N) <- 1
  dimnames(correlation_N) <- list(response_names, response_names)

  # ----------------------------------------------------------
  # Return
  # ----------------------------------------------------------
  out <- list(
    Y = Y,
    X = X,
    Z = Z,
    beta = beta,
    gamma = gamma,
    zero_inflation = zero_inflation,
    zero_inflated_index = zi_index,
    phi = phi,
    mixing.distribution = mixing.distribution,
    var.power = var_power,
    psi_N = psi_N,
    psi_B = psi_B,
    correlation_N = correlation_N,
    correlation_B = correlation_B,
    mu = mu,
    prob = prob,
    n = n,
    k = k,
    J = J,
    n_zi = n_zi
  )

  if (isTRUE(return_latent)) {
    out$N <- N_count
    out$B <- B
    out$Lambda <- Lambda
    out$U_N <- U_N
    out$U_B <- U_B
    out$F_N <- F_N
    out$F_B <- F_B
    out$theoretical_variance <- theoretical_variance
  }

  class(out) <- "simu_multizimp"
  out
}

# ============================================================
# Examples
# ============================================================
#
# ## 1. No zero-inflated margins
# sim0 <- simu_multizimp(
#   X = X, beta = beta,
#   Z = NULL, gamma = NULL,
#   zero_inflation = NULL,
#   mixing.distribution = mixing.distribution,
#   psi_N = psi_N,
#   psi_B = NULL,
#   phi = 1, seed = 100
# )
#
# ## 2. Some zero-inflated margins
# zi <- c(TRUE, FALSE, TRUE, FALSE, TRUE)
# psi_B <- simu_psi_B(sum(zi), seed = 2)
#
# sim1 <- simu_multizimp(
#   X = X, beta = beta,
#   Z = Z, gamma = gamma,
#   zero_inflation = zi,
#   mixing.distribution = mixing.distribution,
#   psi_N = psi_N,
#   psi_B = psi_B,
#   phi = 1, seed = 100
# )
#
# ## 3. All margins zero-inflated
# zi <- rep(TRUE, k)
# psi_B <- simu_psi_B(k, seed = 2)
#
# sim2 <- simu_multizimp(
#   X = X, beta = beta,
#   Z = Z, gamma = gamma,
#   zero_inflation = zi,
#   mixing.distribution = mixing.distribution,
#   psi_N = psi_N,
#   psi_B = psi_B,
#   phi = 1, seed = 100
# )
