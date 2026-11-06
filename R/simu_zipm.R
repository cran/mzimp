simu_zipm <- function(
    n,
    mu,
    phi,
    latent.distribution = c(
      "gamma",
      "lnorm",
      "pareto",
      "weibull",
      "ig",
      "loglogistic",
      "betaprime"
    ),
    use_custom_var_power = FALSE,
    var_power = 2,
    zero_inflation = FALSE,
    prob = NULL,
    seed = NULL,
    uniform.import = FALSE,
    U.import = NULL,
    return_latent = TRUE
) {
  
  # ==========================================================
  # Seed
  # ==========================================================
  
  if (!is.null(seed)) {
    
    if (!is.numeric(seed) ||
        length(seed) != 1L ||
        !is.finite(seed)) {
      stop("seed must be NULL or a finite numeric scalar.")
    }
    
    set.seed(seed)
  }
  
  latent.distribution <- match.arg(latent.distribution)
  
  # ==========================================================
  # Validate n
  # ==========================================================
  
  if (!is.numeric(n) ||
      length(n) != 1L ||
      !is.finite(n) ||
      n <= 0 ||
      abs(n - round(n)) > 1e-12) {
    stop("n must be a positive integer.")
  }
  
  n <- as.integer(round(n))
  
  # ==========================================================
  # Validate mu
  # ==========================================================
  
  if (length(mu) == 1L) {
    mu <- rep(mu, n)
  }
  
  if (!is.numeric(mu) ||
      length(mu) != n ||
      any(!is.finite(mu)) ||
      any(mu <= 0)) {
    stop(
      "mu must be a positive numeric vector of length n, ",
      "or a positive scalar."
    )
  }
  
  mu <- as.numeric(mu)
  
  # ==========================================================
  # Validate phi
  # ==========================================================
  
  if (!is.numeric(phi) ||
      length(phi) != 1L ||
      !is.finite(phi) ||
      phi <= 0) {
    stop("phi must be a positive scalar.")
  }
  
  phi <- as.numeric(phi)
  
  # ==========================================================
  # Variance-power specification
  # ==========================================================
  
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
  
  if (use_custom_var_power) {
    
    if (!(latent.distribution %in%
          distributions_with_var_power)) {
      warning(
        "var_power is not used for latent.distribution = '",
        latent.distribution,
        "'. The distribution-specific parameterization ",
        "is used instead."
      )
    }
    
    if (!is.numeric(var_power) ||
        length(var_power) != 1L ||
        !is.finite(var_power) ||
        abs(var_power - round(var_power)) > 1e-12 ||
        !(as.integer(round(var_power)) %in% c(1L, 2L, 3L))) {
      stop(
        "When use_custom_var_power = TRUE, ",
        "var_power must be 1, 2, or 3."
      )
    }
    
    var_power <- as.integer(round(var_power))
  }
  
  # ----------------------------------------------------------
  # Default variance powers
  #
  # gamma:
  #   Var(Lambda) = phi * mu^2
  #
  # lnorm:
  #   Var(Lambda) = phi * mu^2
  #
  # weibull:
  #   Var(Lambda) = phi * mu^2
  #
  # inverse Gaussian:
  #   Var(Lambda) = phi * mu^3
  #
  # log-logistic:
  #   Var(Lambda) = phi * mu^2
  #
  # beta-prime:
  #   Var(Lambda) = phi * mu^2
  # ----------------------------------------------------------
  
  default_var_power <- switch(
    latent.distribution,
    gamma = 2L,
    lnorm = 2L,
    weibull = 2L,
    ig = 3L,
    loglogistic = 2L,
    betaprime = 2L,
    NA_integer_
)
  
  if (latent.distribution %in%
      distributions_with_var_power) {
    
    effective_var_power <- if (use_custom_var_power) {
      as.integer(var_power)
    } else {
      default_var_power
    }
    
  } else {
    effective_var_power <- NA_integer_
  }
  
  # ==========================================================
  # Zero-inflation probability
  #
  # prob = probability of belonging to the Poisson component
  # ==========================================================
  
  if (!is.logical(zero_inflation) ||
      length(zero_inflation) != 1L ||
      is.na(zero_inflation)) {
    stop("zero_inflation must be TRUE or FALSE.")
  }
  
  if (zero_inflation) {
    
    if (is.null(prob)) {
      stop(
        "When zero_inflation = TRUE, prob must be provided."
      )
    }
    
    if (length(prob) == 1L) {
      prob <- rep(prob, n)
    }
    
    if (!is.numeric(prob) ||
        length(prob) != n ||
        any(!is.finite(prob)) ||
        any(prob < 0) ||
        any(prob > 1)) {
      stop(
        "prob must be a numeric vector in [0,1] of length n, ",
        "or a scalar."
      )
    }
    
    prob <- as.numeric(prob)
    
  } else {
    prob <- rep(1, n)
  }
  
  # ==========================================================
  # Generate or import uniform random variables
  # ==========================================================
  
  if (!is.logical(uniform.import) ||
      length(uniform.import) != 1L ||
      is.na(uniform.import)) {
    stop("uniform.import must be TRUE or FALSE.")
  }
  
  eps <- 1e-12
  
  if (uniform.import) {
    
    if (is.null(U.import)) {
      stop(
        "U.import must be supplied when uniform.import = TRUE."
      )
    }
    
    if (!is.numeric(U.import) ||
        length(U.import) != n) {
      stop(
        "U.import must be a numeric vector of length n."
      )
    }
    
    if (anyNA(U.import) ||
        any(!is.finite(U.import))) {
      stop(
        "U.import contains missing or non-finite values."
      )
    }
    
    if (any(U.import <= 0 | U.import >= 1)) {
      stop(
        "All U.import values must satisfy 0 < U.import < 1."
      )
    }
    
    U <- as.numeric(U.import)
    
  } else {
    U <- runif(n)
  }
  
  U <- pmin(
    1 - eps,
    pmax(eps, U)
  )
  
  # ==========================================================
  # Gamma quantile
  # ==========================================================
  
  qgamma_mu_phi_vp <- function(
    u,
    mu,
    phi,
    var_power
  ) {
    
    variance <- phi * mu^var_power
    
    shape <- mu^2 / variance
    rate <- mu / variance
    
    qgamma(
      p = u,
      shape = shape,
      rate = rate
    )
  }
  
  # ==========================================================
  # Lognormal quantile
  # ==========================================================
  
  qlnorm_mu_phi_vp <- function(
    u,
    mu,
    phi,
    var_power
  ) {
    
    variance <- phi * mu^var_power
    
    sigma2 <- log1p(
      variance / mu^2
    )
    
    meanlog <- log(mu) - 0.5 * sigma2
    
    qlnorm(
      p = u,
      meanlog = meanlog,
      sdlog = sqrt(sigma2)
    )
  }
  
  # ==========================================================
  # Inverse Gaussian quantile
  # ==========================================================
  
  qig_mu_phi_vp <- function(
    u,
    mu,
    phi,
    var_power
  ) {
    
    if (!requireNamespace(
      "statmod",
      quietly = TRUE
    )) {
      stop(
        "Package 'statmod' is required for ",
        "inverse Gaussian simulation."
      )
    }
    
    variance <- phi * mu^var_power
    
    shape <- mu^3 / variance
    
    statmod::qinvgauss(
      p = u,
      mean = mu,
      shape = shape
    )
  }
  
  # ==========================================================
  # Weibull parameter conversion
  #
  # E(Lambda)   = mu
  # Var(Lambda) = phi * mu^var_power
  #
  # Let
  #
  #   a = 1 / shape.
  #
  # Solve
  #
  #   lgamma(1 + 2a) - 2 lgamma(1 + a)
  #   = log(1 + phi * mu^(var_power - 2)).
  # ==========================================================
  
  weibull_log_moment_ratio <- function(a) {
    lgamma(1 + 2 * a) -
      2 * lgamma(1 + a)
  }
  
  weibull_log_moment_ratio_deriv <- function(a) {
    2 * (
      digamma(1 + 2 * a) -
        digamma(1 + a)
    )
  }
  
  solve_weibull_a <- function(
    target_cv2,
    tol = 1e-12,
    maxit = 100L
  ) {
    
    if (!is.numeric(target_cv2) ||
        length(target_cv2) != 1L ||
        !is.finite(target_cv2) ||
        target_cv2 <= 0) {
      stop(
        "The target squared coefficient of variation ",
        "must be positive."
      )
    }
    
    target_log_ratio <- log1p(target_cv2)
    
    lower <- 0
    upper <- 1
    
    g_upper <- weibull_log_moment_ratio(upper)
    
    while (g_upper < target_log_ratio) {
      
      upper <- 2 * upper
      
      if (!is.finite(upper) || upper > 1e6) {
        stop(
          "Unable to bracket the Weibull shape parameter."
        )
      }
      
      g_upper <- weibull_log_moment_ratio(upper)
      
      if (!is.finite(g_upper)) {
        stop(
          "Non-finite value encountered while solving ",
          "the Weibull parameters."
        )
      }
    }
    
    a <- 0.5 * (lower + upper)
    
    for (iter in seq_len(maxit)) {
      
      g <- weibull_log_moment_ratio(a)
      
      difference <- g - target_log_ratio
      
      if (abs(difference) <
          tol * (1 + abs(target_log_ratio))) {
        return(a)
      }
      
      if (difference < 0) {
        lower <- a
      } else {
        upper <- a
      }
      
      gp <- weibull_log_moment_ratio_deriv(a)
      
      if (is.finite(gp) && gp > 0) {
        a_new <- a - difference / gp
      } else {
        a_new <- 0.5 * (lower + upper)
      }
      
      if (!is.finite(a_new) ||
          a_new <= lower ||
          a_new >= upper) {
        a_new <- 0.5 * (lower + upper)
      }
      
      if (abs(a_new - a) <
          tol * (1 + abs(a))) {
        return(a_new)
      }
      
      a <- a_new
    }
    
    a
  }
  
  weibull_mu_phi_parameters_R <- function(
    mu,
    phi,
    var_power
  ) {
    
    cv2 <- phi * mu^(var_power - 2)
    
    if (any(!is.finite(cv2)) ||
        any(cv2 <= 0)) {
      stop(
        "Invalid Weibull squared coefficient of variation."
      )
    }
    
    # Avoid repeatedly solving identical parameter values
    unique_cv2 <- unique(cv2)
    
    unique_a <- vapply(
      unique_cv2,
      solve_weibull_a,
      numeric(1)
    )
    
    a <- unique_a[
      match(cv2, unique_cv2)
    ]
    
    shape <- 1 / a
    
    log_scale <-
      log(mu) -
      lgamma(1 + a)
    
    scale <- exp(log_scale)
    
    data.frame(
      mu = mu,
      variance = phi * mu^var_power,
      cv2 = cv2,
      a = a,
      shape = shape,
      scale = scale,
      log_scale = log_scale
    )
  }
  
  # ==========================================================
  # Weibull quantile
  #
  # Prefer the existing C++ parameter conversion when it has
  # already been compiled. Otherwise use the equivalent R code.
  # ==========================================================
  
  qweibull_mu_phi_vp <- function(
    u,
    mu,
    phi,
    var_power
  ) {
    
    if (exists(
      "weibull_mu_phi_parameters_cpp",
      mode = "function",
      inherits = TRUE
    )) {
      
      parameters <- weibull_mu_phi_parameters_cpp(
        mu = mu,
        phi = phi,
        var_power = var_power
      )
      
      shape <- as.numeric(
        parameters[, "shape"]
      )
      
      scale <- as.numeric(
        parameters[, "scale"]
      )
      
    } else {
      
      parameters <- weibull_mu_phi_parameters_R(
        mu = mu,
        phi = phi,
        var_power = var_power
      )
      
      shape <- parameters$shape
      scale <- parameters$scale
    }
    
    qweibull(
      p = u,
      shape = shape,
      scale = scale
    )
  }
  
  # ==========================================================
  # Log-logistic parameter conversion
  #
  # E(Lambda)   = mu
  # Var(Lambda) = phi * mu^var_power
  #
  # Let a = 1/kappa, with 0 < a < 1/2. Then
  #
  #   A1(a) = pi*a/sin(pi*a),
  #   A2(a) = 2*pi*a/sin(2*pi*a),
  #
  # and a solves
  #
  #   log A2(a) - 2 log A1(a)
  #   = log(1 + phi * mu^(var_power - 2)).
  #
  # The scale parameter is
  #
  #   beta = mu * sin(pi*a)/(pi*a).
  # ==========================================================
  
  loglogistic_log_moment_ratio <- function(a) {

    if (any(!is.finite(a)) ||
        any(a <= 0) ||
        any(a >= 0.5)) {
      stop("a must lie in (0, 0.5).")
    }

    x <- pi * a

    # Algebraically,
    #
    #   A2(a) / A1(a)^2 = tan(pi*a) / (pi*a).
    #
    # For small x, use a series to avoid cancellation:
    #
    #   log{tan(x)/x}
    #   = x^2/3 + 7*x^4/90 + 62*x^6/2835 + O(x^8).

    small <- abs(x) < 1e-4
    out <- numeric(length(x))

    out[small] <-
      x[small]^2 / 3 +
      7 * x[small]^4 / 90 +
      62 * x[small]^6 / 2835

    out[!small] <-
      log(tan(x[!small])) -
      log(x[!small])

    out
  }

  loglogistic_log_moment_ratio_deriv <- function(a) {

    if (!is.numeric(a) ||
        length(a) != 1L ||
        !is.finite(a) ||
        a <= 0 ||
        a >= 0.5) {
      stop("a must lie in (0, 0.5).")
    }

    x <- pi * a

    if (abs(x) < 1e-4) {
      return(
        2 * pi^2 * a / 3 +
        14 * pi^4 * a^3 / 90 +
        372 * pi^6 * a^5 / 2835
      )
    }

    # d/da log{tan(pi*a)/(pi*a)}
    pi / (sin(x) * cos(x)) - 1 / a
  }
  
  solve_loglogistic_a <- function(target_cv2, tol = 1e-12, maxit = 100L) {
    if (!is.numeric(target_cv2) || length(target_cv2) != 1L ||
        !is.finite(target_cv2) || target_cv2 <= 0) {
      stop("The target squared coefficient of variation must be positive.")
    }
    target_log_ratio <- log1p(target_cv2)
    lower <- 1e-10
    upper <- 0.5 - 1e-10
    f_lower <- loglogistic_log_moment_ratio(lower) - target_log_ratio
    f_upper <- loglogistic_log_moment_ratio(upper) - target_log_ratio
    if (!is.finite(f_lower) || !is.finite(f_upper) ||
        f_lower > 0 || f_upper < 0) {
      stop("Unable to bracket the log-logistic shape parameter.")
    }
    a <- 0.25
    for (iter in seq_len(maxit)) {
      f <- loglogistic_log_moment_ratio(a) - target_log_ratio
      if (abs(f) < tol * (1 + abs(target_log_ratio))) return(a)
      if (f < 0) lower <- a else upper <- a
      fp <- loglogistic_log_moment_ratio_deriv(a)
      if (is.finite(fp) && fp > 0) {
        a_new <- a - f / fp
      } else {
        a_new <- 0.5 * (lower + upper)
      }
      if (!is.finite(a_new) || a_new <= lower || a_new >= upper) {
        a_new <- 0.5 * (lower + upper)
      }
      if (abs(a_new - a) < tol * (1 + abs(a))) return(a_new)
      a <- a_new
    }
    a
  }
  
  loglogistic_mu_phi_parameters_R <- function(mu, phi, var_power) {
    cv2 <- phi * mu^(var_power - 2)
    if (any(!is.finite(cv2)) || any(cv2 <= 0)) {
      stop("Invalid log-logistic squared coefficient of variation.")
    }
    unique_cv2 <- unique(cv2)
    unique_a <- vapply(unique_cv2, solve_loglogistic_a, numeric(1))
    a <- unique_a[match(cv2, unique_cv2)]
    kappa <- 1 / a
    log_beta <- log(mu) + log(sin(pi * a)) - log(pi) - log(a)
    beta <- exp(log_beta)
    data.frame(mu=mu, variance=phi*mu^var_power, cv2=cv2,
               a=a, kappa=kappa, beta=beta, log_beta=log_beta)
  }
  
  qloglogistic_mu_phi_vp <- function(u, mu, phi, var_power) {
    if (exists("loglogistic_mu_phi_parameters_cpp", mode="function", inherits=TRUE)) {
      parameters <- loglogistic_mu_phi_parameters_cpp(
        mu=mu, phi=phi, var_power=var_power
      )
      a <- as.numeric(parameters[, "a"])
      beta <- as.numeric(parameters[, "beta"])
    } else {
      parameters <- loglogistic_mu_phi_parameters_R(
        mu=mu, phi=phi, var_power=var_power
      )
      a <- parameters$a
      beta <- parameters$beta
    }
    log_quantile <- log(beta) + a * (log(u) - log1p(-u))
    exp(log_quantile)
  }
  


  # ==========================================================
  # Symmetric Beta-prime quantile
  #
  # Lambda ~ BetaPrime(a, a, scale = s), with
  #
  #   E(Lambda)   = mu,
  #   Var(Lambda) = phi * mu^var_power.
  #
  # Let
  #
  #   cv2 = Var(Lambda) / mu^2
  #       = phi * mu^(var_power - 2).
  #
  # For a symmetric Beta-prime distribution,
  #
  #   cv2 = (2a - 1) / {a(a - 2)},   a > 2.
  #
  # Solving for a gives
  #
  #   a = {cv2 + 1 + sqrt(cv2^2 + cv2 + 1)} / cv2.
  #
  # The scale is then
  #
  #   s = mu * (a - 1) / a.
  #
  # If T ~ Beta(a,a), then
  #
  #   Lambda = s * T / (1 - T).
  # ==========================================================

  qbetaprime_mu_phi_vp <- function(
    u,
    mu,
    phi,
    var_power
  ) {

    cv2 <- phi * mu^(var_power - 2)

    if (any(!is.finite(cv2)) ||
        any(cv2 <= 0)) {
      stop(
        "Invalid Beta-prime squared coefficient of variation."
      )
    }

    shape <- (
      cv2 + 1 + sqrt(cv2^2 + cv2 + 1)
    ) / cv2

    if (any(!is.finite(shape)) ||
        any(shape <= 2)) {
      stop(
        "Invalid symmetric Beta-prime shape parameter."
      )
    }

    scale <- mu * (shape - 1) / shape

    beta_quantile <- qbeta(
      p = u,
      shape1 = shape,
      shape2 = shape
    )

    beta_quantile <- pmin(
      1 - eps,
      pmax(eps, beta_quantile)
    )

    log_quantile <-
      log(scale) +
      log(beta_quantile) -
      log1p(-beta_quantile)

    exp(log_quantile)
  }

  # ==========================================================
  # Pareto Type II / Lomax quantile
  #
  # E(Lambda)   = mu
  # Var(Lambda) = mu^2 * (1 + 2 * phi)
  #
  # Define
  #
  #   cv2 = (1 + 2 * phi) * mu^(2 - 2).
  #
  # For a Lomax distribution,
  #
  #   cv2 = alpha / (alpha - 2).
  #
  # Hence
  #
  #   alpha = 2 * cv2 / (cv2 - 1),
  #   beta  = mu * (alpha - 1).
  #
  # The quantile is
  #
  #   Q(u) = beta * {(1-u)^(-1/alpha) - 1}.
  #
  # A finite-variance Lomax distribution requires cv2 > 1.
  # ==========================================================
  
  qpareto_mu_phi_vp <- function(
    u,
    mu,
    phi
  ) {
    alpha <- 2 + 1 / phi
    beta  <- mu * (1 + 1 / phi)
    beta * ((1 - u)^(-1 / alpha) - 1)
  }
  
  # ==========================================================
  # Generate latent intensity
  # ==========================================================
  
  Lambda <- switch(
    latent.distribution,
    
    gamma = qgamma_mu_phi_vp(
      u = U,
      mu = mu,
      phi = phi,
      var_power = effective_var_power
    ),
    
    lnorm = qlnorm_mu_phi_vp(
      u = U,
      mu = mu,
      phi = phi,
      var_power = effective_var_power
    ),
    
    ig = qig_mu_phi_vp(
      u = U,
      mu = mu,
      phi = phi,
      var_power = effective_var_power
    ),
    
    weibull = qweibull_mu_phi_vp(
      u = U,
      mu = mu,
      phi = phi,
      var_power = effective_var_power
    ),
    
    pareto = qpareto_mu_phi_vp(
      u = U,
      mu = mu,
      phi = phi
    ),
    
    loglogistic = qloglogistic_mu_phi_vp(
      u = U,
      mu = mu,
      phi = phi,
      var_power = effective_var_power
    ),

    betaprime = qbetaprime_mu_phi_vp(
      u = U,
      mu = mu,
      phi = phi,
      var_power = effective_var_power
    )
  )
  
  Lambda <- as.numeric(Lambda)
  
  if (any(!is.finite(Lambda)) ||
      any(Lambda < 0)) {
    stop(
      "The generated latent intensities contain invalid values. ",
      "Check mu, phi, U.import, var_power, and the ",
      "latent distribution."
    )
  }
  
  # ==========================================================
  # Poisson and zero-inflation layers
  # ==========================================================
  
  Z <- rpois(
    n = n,
    lambda = Lambda
  )
  
  if (zero_inflation) {
    
    B <- rbinom(
      n = n,
      size = 1L,
      prob = prob
    )
    
    Y <- B * Z
    
  } else {
    
    B <- rep(1L, n)
    Y <- Z
  }
  
  # ==========================================================
  # Return result
  # ==========================================================
  
  theoretical_variance <- switch(
    latent.distribution,
    
    gamma =
      phi * mu^effective_var_power,
    
    lnorm =
      phi * mu^effective_var_power,
    
    pareto =
      mu^2 * (1 + 2 * phi),
    
    weibull =
      phi * mu^effective_var_power,
    
    ig =
      phi * mu^effective_var_power,
    
    loglogistic =
      phi * mu^effective_var_power,

    betaprime =
      phi * mu^effective_var_power
  )
  
  
  
  
  if (return_latent) {
    
    out <- data.frame(
      Y = Y,
      N = Z,
      B = B,
      Lambda = Lambda,
      U = U,
      mu = mu,
      theoretical_variance = theoretical_variance,
      phi = rep(phi, n),
      prob = prob,
      latent.distribution = rep(
        latent.distribution,
        n
      ),
      use_custom_var_power = rep(
        use_custom_var_power,
        n
      ),
      var_power = rep(
        effective_var_power,
        n
      ),
      zero_inflation = rep(
        zero_inflation,
        n
      ),
      stringsAsFactors = FALSE
    )
    
  } else {
    out <- Y
  }
  
  attr(out, "latent.distribution") <-
    latent.distribution
  
  attr(out, "phi") <- phi
  
  attr(out, "use_custom_var_power") <-
    use_custom_var_power
  
  attr(out, "var_power") <-
    effective_var_power
  
  attr(out, "default_var_power") <-
    default_var_power
  
  attr(out, "zero_inflation") <-
    zero_inflation
  
  out
}