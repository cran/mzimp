simu_biv_pm <- function(
    n,
    mu1, phi1,
    mu2, phi2,
    dist1 = c(
      "gamma", "lnorm", "pareto",
      "weibull", "ig", "loglogistic", "betaprime"
    ),
    dist2 = c(
      "gamma", "lnorm", "pareto",
      "weibull", "ig", "loglogistic", "betaprime"
    ),
    use_custom_var_power = FALSE,
    var_power = c(2, 2),
    family = 1,
    par = 0.0,
    par2 = 3,
    d = 10,
    seed = NULL,
    return_latent = TRUE,
    checkboard = TRUE
) {

  # -------------------------------------------------------
  # Seed
  # -------------------------------------------------------
  if (!is.null(seed)) {
    if (!is.numeric(seed) ||
        length(seed) != 1L ||
        !is.finite(seed)) {
      stop("seed must be NULL or a finite numeric scalar.")
    }
    set.seed(seed)
  }

  # -------------------------------------------------------
  # Match marginal distributions
  # -------------------------------------------------------
  dist1 <- match.arg(dist1)
  dist2 <- match.arg(dist2)

  # -------------------------------------------------------
  # Basic checks
  # -------------------------------------------------------
  if (!is.numeric(n) ||
      length(n) != 1L ||
      !is.finite(n) ||
      n <= 0 ||
      abs(n - round(n)) > 1e-12) {
    stop("n must be a positive integer.")
  }
  n <- as.integer(round(n))

  if (length(mu1) == 1L) mu1 <- rep(mu1, n)
  if (!is.numeric(mu1) ||
      length(mu1) != n ||
      any(!is.finite(mu1)) ||
      any(mu1 <= 0)) {
    stop(
      "mu1 must be a positive numeric vector of length n, ",
      "or a positive scalar."
    )
  }

  if (length(mu2) == 1L) mu2 <- rep(mu2, n)
  if (!is.numeric(mu2) ||
      length(mu2) != n ||
      any(!is.finite(mu2)) ||
      any(mu2 <= 0)) {
    stop(
      "mu2 must be a positive numeric vector of length n, ",
      "or a positive scalar."
    )
  }

  if (!is.numeric(phi1) ||
      length(phi1) != 1L ||
      !is.finite(phi1) ||
      phi1 <= 0) {
    stop("phi1 must be a positive numeric scalar.")
  }

  if (!is.numeric(phi2) ||
      length(phi2) != 1L ||
      !is.finite(phi2) ||
      phi2 <= 0) {
    stop("phi2 must be a positive numeric scalar.")
  }

  if (!is.logical(checkboard) ||
      length(checkboard) != 1L ||
      is.na(checkboard)) {
    stop("checkboard must be TRUE or FALSE.")
  }

  if (!is.logical(return_latent) ||
      length(return_latent) != 1L ||
      is.na(return_latent)) {
    stop("return_latent must be TRUE or FALSE.")
  }

  if (!is.numeric(d) ||
      length(d) != 1L ||
      !is.finite(d) ||
      d <= 0 ||
      abs(d - round(d)) > 1e-12) {
    stop("d must be a positive integer.")
  }
  d <- as.integer(round(d))

  # -------------------------------------------------------
  # Variance-power specification
  #
  # Supported distributions:
  # gamma, lnorm, weibull, ig, loglogistic, betaprime
  # -------------------------------------------------------
  if (!is.logical(use_custom_var_power) ||
      length(use_custom_var_power) != 1L ||
      is.na(use_custom_var_power)) {
    stop("use_custom_var_power must be TRUE or FALSE.")
  }

  distributions_with_var_power <- c(
    "gamma",
    "lnorm",
    "weibull",
    "ig",
    "loglogistic",
    "betaprime"
  )

  default_var_power <- function(dist) {
    switch(
      dist,
      gamma = 2L,
      lnorm = 2L,
      weibull = 2L,
      ig = 3L,
      loglogistic = 2L,
      betaprime = 2L,
      NA_integer_
    )
  }

  if (use_custom_var_power) {
    if (length(var_power) == 1L) {
      var_power <- rep(var_power, 2L)
    }

    if (!is.numeric(var_power) ||
        length(var_power) != 2L ||
        any(!is.finite(var_power)) ||
        any(abs(var_power - round(var_power)) > 1e-12) ||
        any(!(as.integer(round(var_power)) %in% c(1L, 2L, 3L)))) {
      stop(
        "When use_custom_var_power = TRUE, var_power must be ",
        "a scalar or a numeric vector of length 2 containing ",
        "only 1, 2, or 3."
      )
    }

    var_power <- as.integer(round(var_power))

    if (!(dist1 %in% distributions_with_var_power)) {
      warning(
        "var_power[1] is not used for dist1 = '",
        dist1,
        "'. The distribution-specific parameterization is used."
      )
    }

    if (!(dist2 %in% distributions_with_var_power)) {
      warning(
        "var_power[2] is not used for dist2 = '",
        dist2,
        "'. The distribution-specific parameterization is used."
      )
    }
  } else {
    var_power <- c(
      default_var_power(dist1),
      default_var_power(dist2)
    )
  }

  # -------------------------------------------------------
  # Copula parameter checks
  # -------------------------------------------------------
  allowed_family <- c(
    0L, 1L, 2L,
    3L, 13L, 23L, 33L,
    4L, 14L, 24L, 34L,
    5L, 25L
  )

  if (!is.numeric(family) ||
      length(family) != 1L ||
      !is.finite(family) ||
      abs(family - round(family)) > 1e-12) {
    stop("family must be an integer copula family code.")
  }
  family <- as.integer(round(family))

  if (!(family %in% allowed_family)) {
    stop(
      "Unsupported copula family. Supported codes are: ",
      paste(allowed_family, collapse = ", "),
      "."
    )
  }

  if (family != 0L) {
    if (!is.numeric(par) ||
        length(par) != 1L ||
        !is.finite(par)) {
      stop("par must be a finite numeric scalar.")
    }
  }

  if (family %in% c(1L, 2L) && (par <= -1 || par >= 1)) {
    stop(
      "For Gaussian and Student-t copulas, ",
      "par must satisfy -1 < par < 1."
    )
  }

  if (family == 2L) {
    if (!is.numeric(par2) ||
        length(par2) != 1L ||
        !is.finite(par2) ||
        par2 <= 0) {
      stop(
        "For the Student-t copula, par2 must be ",
        "a positive degrees-of-freedom parameter."
      )
    }
  }

  if (family %in% c(3L, 13L, 23L, 33L) && par <= 0) {
    stop("For Clayton copulas, par must be positive.")
  }

  if (family %in% c(4L, 14L, 24L, 34L) && par < 1) {
    stop("For Gumbel copulas, par must be at least 1.")
  }

  if (family %in% c(5L, 25L) && abs(par) < 1e-10) {
    stop(
      "For Frank copulas, par must be nonzero; ",
      "use family = 0 for independence."
    )
  }

  # -------------------------------------------------------
  # Step 1: simulate bivariate uniforms
  # -------------------------------------------------------
  if (checkboard) {
    if (!exists(
      "simulate_checkerboard_bicop",
      mode = "function",
      inherits = TRUE
    )) {
      stop("Function 'simulate_checkerboard_bicop' was not found.")
    }

    U <- simulate_checkerboard_bicop(
      n = n,
      family = family,
      par = par,
      par2 = par2,
      d = d
    )
  } else {
    if (!requireNamespace("VineCopula", quietly = TRUE)) {
      stop(
        "Package 'VineCopula' is required when ",
        "checkboard = FALSE."
      )
    }

    U <- VineCopula::BiCopSim(
      N = n,
      family = family,
      par = par,
      par2 = par2
    )
  }

  # -------------------------------------------------------
  # Validate simulated uniforms
  # -------------------------------------------------------
  U <- as.matrix(U)

  if (nrow(U) != n || ncol(U) < 2L) {
    stop("The copula simulation must return an n x 2 matrix.")
  }

  U <- U[, 1:2, drop = FALSE]

  if (any(!is.finite(U))) {
    stop("The simulated copula uniforms contain non-finite values.")
  }

  eps <- 1e-12
  U1 <- pmin(1 - eps, pmax(eps, U[, 1L]))
  U2 <- pmin(1 - eps, pmax(eps, U[, 2L]))

  # -------------------------------------------------------
  # Step 2: simulate mixed-Poisson margins
  # -------------------------------------------------------
  margin1 <- simu_zipm(
    n = n,
    mu = mu1,
    phi = phi1,
    latent.distribution = dist1,
    use_custom_var_power = use_custom_var_power,
    var_power = var_power[1L],
    zero_inflation = FALSE,
    prob = NULL,
    seed = NULL,
    uniform.import = TRUE,
    U.import = U1,
    return_latent = TRUE
  )

  margin2 <- simu_zipm(
    n = n,
    mu = mu2,
    phi = phi2,
    latent.distribution = dist2,
    use_custom_var_power = use_custom_var_power,
    var_power = var_power[2L],
    zero_inflation = FALSE,
    prob = NULL,
    seed = NULL,
    uniform.import = TRUE,
    U.import = U2,
    return_latent = TRUE
  )

  effective_var_power1 <- attr(margin1, "var_power")
  effective_var_power2 <- attr(margin2, "var_power")

  # -------------------------------------------------------
  # Step 3: construct output
  # -------------------------------------------------------
  out <- data.frame(
    N1 = margin1$Y,
    N2 = margin2$Y
  )

  if (return_latent) {
    out$G1 <- margin1$Lambda
    out$G2 <- margin2$Lambda
    out$U1 <- margin1$U
    out$U2 <- margin2$U
    out$mu1 <- margin1$mu
    out$mu2 <- margin2$mu
    out$theoretical_variance1 <- margin1$theoretical_variance
    out$theoretical_variance2 <- margin2$theoretical_variance
  }

  # -------------------------------------------------------
  # Attributes
  # -------------------------------------------------------
  attr(out, "copula_family") <- family
  attr(out, "par") <- par
  attr(out, "par2") <- par2
  attr(out, "checkboard") <- checkboard
  attr(out, "checkerboard_d") <- if (checkboard) d else NULL
  attr(out, "copula_type") <- if (checkboard) {
    "checkerboard"
  } else {
    "continuous"
  }
  attr(out, "use_custom_var_power") <- use_custom_var_power
  attr(out, "var_power") <- c(
    effective_var_power1,
    effective_var_power2
  )

  attr(out, "latent_margins") <- list(
    margin1 = list(
      distribution = dist1,
      mu = mu1,
      phi = phi1,
      use_custom_var_power = use_custom_var_power,
      var_power = effective_var_power1
    ),
    margin2 = list(
      distribution = dist2,
      mu = mu2,
      phi = phi2,
      use_custom_var_power = use_custom_var_power,
      var_power = effective_var_power2
    )
  )

  out
}
