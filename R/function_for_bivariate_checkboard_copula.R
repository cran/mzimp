# ============================================================
# Bivariate copula families using VineCopula family codes
# ============================================================
# family = 0  : Independence copula; no parameter; independent case
# family = 1  : Gaussian copula; par = rho in (-1,1); independence at par = 0
# family = 2  : Student-t copula; par = rho in (-1,1), par2 = nu > 0 fixed;
#               independence at par = 0
# family = 3  : Clayton copula; par > 0; independence as par -> 0
# family = 13 : Clayton copula rotated 180 degrees; par > 0; independence as par -> 0
# family = 23 : Clayton copula rotated 90 degrees; par > 0; independence as par -> 0
# family = 33 : Clayton copula rotated 270 degrees; par > 0; independence as par -> 0
# family = 4  : Gumbel copula; par >= 1; independence at par = 1
# family = 14 : Gumbel copula rotated 180 degrees; par >= 1; independence at par = 1
# family = 24 : Gumbel copula rotated 90 degrees; par >= 1; independence at par = 1
# family = 34 : Gumbel copula rotated 270 degrees; par >= 1; independence at par = 1
# family = 5  : Frank copula; par in R \ {0}; independence as par -> 0
# family = 25 : Frank copula rotated 90 degrees; par in R \ {0}; independence as par -> 0
#
# Output:
#   W        : d x d matrix of checkerboard cell probabilities
#   density  : d x d matrix of checkerboard densities, equal to d^2 W
#   dW       : derivative of W with respect to par
#   ddensity : derivative of density with respect to par
# ============================================================


# ============================================================
# 1. Helpers
# ============================================================




one_par <- c(
  0,
  1,
  2,
  3, 13, 23, 33,
  4, 14, 24, 34,
  5,
  6, 16, 26, 36
)

two_par <- c(
  7, 17, 27, 37,
  8, 18, 28, 38,
  9, 19, 29, 39,
  10, 20, 30, 40,
  104, 114, 124, 134,
  204, 214, 224, 234
)


copula.list <- data.frame( copula.famity = c(one_par, two_par), number.par = c(rep(1, length(one_par)), rep(2, length(two_par))))




.safe_u <- function(u, eps = 1e-12) {
  pmin(1 - eps, pmax(eps, u))
}

.base_family <- function(family) {
  if (family %in% c(3, 13, 23, 33)) return(3)
  if (family %in% c(4, 14, 24, 34)) return(4)
  if (family %in% c(5, 25)) return(5)
  family
}

.rotation_type <- function(family) {
  if (family %in% c(13, 14)) return(180)
  if (family %in% c(23, 24, 25)) return(90)
  if (family %in% c(33, 34)) return(270)
  return(0)
}

.check_par <- function(family, par = NULL, par2 = 3) {
  if (family == 0) return(invisible(TRUE))
  
  base <- .base_family(family)
  
  if (base %in% c(1, 2)) {
    if (is.null(par) || !is.finite(par) || abs(par) >= 1) {
      stop("For Gaussian/Student copula, par must be rho in (-1,1).")
    }
  }
  
  if (base == 2) {
    if (is.null(par2) || !is.finite(par2) || par2 <= 0) {
      stop("For Student-t copula, par2 must be the fixed degrees of freedom nu > 0.")
    }
  }
  
  if (base == 3) {
    if (is.null(par) || !is.finite(par) || par <= 0) {
      stop("For Clayton copula, par must be positive.")
    }
  }
  
  if (base == 4) {
    if (is.null(par) || !is.finite(par) || par < 1) {
      stop("For Gumbel copula, par must be >= 1.")
    }
  }
  
  if (base == 5) {
    if (is.null(par) || !is.finite(par)) {
      stop("For Frank copula, par must be finite.")
    }
  }
  
  invisible(TRUE)
}


# ============================================================
# 2. Base CDF and analytic gradient
# ============================================================

.base_copula_cdf <- function(u, v, family, par = NULL, par2 = 3) {
  .check_par(family, par, par2)
  
  if (u <= 0 || v <= 0) return(0)
  if (u >= 1 && v >= 1) return(1)
  if (family == 0) return(u * v)
  
  u <- .safe_u(u)
  v <- .safe_u(v)
  
  if (family == 1) {
    rho <- par
    x <- qnorm(u)
    y <- qnorm(v)
    
    return(as.numeric(mvtnorm::pmvnorm(
      lower = c(-Inf, -Inf),
      upper = c(x, y),
      sigma = matrix(c(1, rho, rho, 1), 2, 2)
    )))
  }
  
  if (family == 2) {
    rho <- par
    nu  <- par2
    
    x <- qt(u, df = nu)
    y <- qt(v, df = nu)
    
    Sigma <- matrix(
      c(1, rho,
        rho, 1),
      2, 2
    )
    
    return(as.numeric(mvtnorm::pmvt(
      lower = c(-Inf, -Inf),
      upper = c(x, y),
      sigma = Sigma,
      df = nu
    )))
  }
  
  if (family == 3) {
    theta <- par
    S <- u^(-theta) + v^(-theta) - 1
    return(S^(-1 / theta))
  }
  
  if (family == 4) {
    theta <- par
    a <- -log(u)
    b <- -log(v)
    S <- a^theta + b^theta
    return(exp(-S^(1 / theta)))
  }
  
  if (family == 5) {
    theta <- par
    if (abs(theta) < 1e-8) return(u * v)
    
    A <- expm1(-theta * u)
    B <- expm1(-theta * v)
    D <- expm1(-theta)
    
    H <- 1 + A * B / D
    return(-log(H) / theta)
  }
  
  stop("Unsupported base family.")
}

.base_copula_grad <- function(u, v, family, par = NULL, par2 = 3) {
  .check_par(family, par, par2)
  
  if (u <= 0 || v <= 0) return(0)
  if (u >= 1 && v >= 1) return(0)
  if (family == 0) return(0)
  
  u <- .safe_u(u)
  v <- .safe_u(v)
  
  if (family == 1) {
    rho <- par
    x <- qnorm(u)
    y <- qnorm(v)
    
    den <- 2 * pi * sqrt(1 - rho^2)
    expo <- -0.5 / (1 - rho^2) *
      (x^2 - 2 * rho * x * y + y^2)
    
    return(exp(expo) / den)
  }
  
  if (family == 2) {
    rho <- par
    nu  <- par2
    
    x <- qt(u, df = nu)
    y <- qt(v, df = nu)
    
    Sigma <- matrix(
      c(1, rho,
        rho, 1),
      2, 2
    )
    
    den <- mvtnorm::dmvt(
      x = cbind(x, y),
      sigma = Sigma,
      df = nu,
      log = FALSE
    )
    
    return(as.numeric(den))
  }
  
  if (family == 3) {
    theta <- par
    S <- u^(-theta) + v^(-theta) - 1
    C <- S^(-1 / theta)
    
    S_prime <- -u^(-theta) * log(u) -
      v^(-theta) * log(v)
    
    return(C * (log(S) / theta^2 - S_prime / (theta * S)))
  }
  
  if (family == 4) {
    theta <- par
    
    a <- -log(u)
    b <- -log(v)
    
    S <- a^theta + b^theta
    T <- S^(1 / theta)
    C <- exp(-T)
    
    term_a <- ifelse(a > 0, a^theta * log(a), 0)
    term_b <- ifelse(b > 0, b^theta * log(b), 0)
    
    S_prime <- term_a + term_b
    
    dT <- T * (
      S_prime / (theta * S) -
        log(S) / theta^2
    )
    
    return(-C * dT)
  }
  
  if (family == 5) {
    theta <- par
    
    if (abs(theta) < 1e-6) {
      h <- 1e-5
      return((
        .base_copula_cdf(u, v, 5, theta + h) -
          .base_copula_cdf(u, v, 5, theta - h)
      ) / (2 * h))
    }
    
    A <- expm1(-theta * u)
    B <- expm1(-theta * v)
    D <- expm1(-theta)
    
    dA <- -u * exp(-theta * u)
    dB <- -v * exp(-theta * v)
    dD <- -exp(-theta)
    
    Q <- A * B / D
    H <- 1 + Q
    
    dQ <- (dA * B + A * dB) / D -
      A * B * dD / D^2
    
    dC <- log(H) / theta^2 -
      dQ / (theta * H)
    
    return(dC)
  }
  
  stop("Unsupported base family.")
}


# ============================================================
# 3. Rotated CDF and gradient
# ============================================================

bicop_cdf_analytic <- function(u, v, family, par = NULL, par2 = 3) {
  if (family == 0) return(u * v)
  
  .check_par(family, par, par2)
  
  base <- .base_family(family)
  rot  <- .rotation_type(family)
  
  if (rot == 0) {
    return(.base_copula_cdf(u, v, base, par, par2))
  }
  
  if (rot == 90) {
    return(v - .base_copula_cdf(1 - u, v, base, par, par2))
  }
  
  if (rot == 180) {
    return(u + v - 1 + .base_copula_cdf(1 - u, 1 - v, base, par, par2))
  }
  
  if (rot == 270) {
    return(u - .base_copula_cdf(u, 1 - v, base, par, par2))
  }
  
  stop("Unsupported rotation.")
}

bicop_cdf_grad_analytic <- function(u, v, family, par = NULL, par2 = 3) {
  if (family == 0) return(0)
  
  .check_par(family, par, par2)
  
  base <- .base_family(family)
  rot  <- .rotation_type(family)
  
  if (rot == 0) {
    return(.base_copula_grad(u, v, base, par, par2))
  }
  
  if (rot == 90) {
    return(-.base_copula_grad(1 - u, v, base, par, par2))
  }
  
  if (rot == 180) {
    return(.base_copula_grad(1 - u, 1 - v, base, par, par2))
  }
  
  if (rot == 270) {
    return(-.base_copula_grad(u, 1 - v, base, par, par2))
  }
  
  stop("Unsupported rotation.")
}


# ============================================================
# 4. Bivariate checkerboard density and gradient
# ============================================================

checkerboard_density_bicop <- function(family, par = NULL, par2 = 3, d = 10) {
  .check_par(family, par, par2)
  
  W  <- matrix(0, d, d)
  dW <- matrix(0, d, d)
  
  for (m in 1:d) {
    for (n in 1:d) {
      u0 <- (m - 1) / d
      u1 <- m / d
      v0 <- (n - 1) / d
      v1 <- n / d
      
      C11 <- bicop_cdf_analytic(u1, v1, family, par, par2)
      C01 <- bicop_cdf_analytic(u0, v1, family, par, par2)
      C10 <- bicop_cdf_analytic(u1, v0, family, par, par2)
      C00 <- bicop_cdf_analytic(u0, v0, family, par, par2)
      
      W[m, n] <- C11 - C01 - C10 + C00
      
      G11 <- bicop_cdf_grad_analytic(u1, v1, family, par, par2)
      G01 <- bicop_cdf_grad_analytic(u0, v1, family, par, par2)
      G10 <- bicop_cdf_grad_analytic(u1, v0, family, par, par2)
      G00 <- bicop_cdf_grad_analytic(u0, v0, family, par, par2)
      
      dW[m, n] <- G11 - G01 - G10 + G00
    }
  }
  
  W[W < 0 & W > -1e-12] <- 0
  
  list(
    W        = W,
    density  = d^2 * W,
    dW       = dW,
    ddensity = d^2 * dW,
    family   = family,
    par      = par,
    par2     = par2,
    d        = d
  )
}


# ============================================================
# 5. Simulate from checkerboard copula
# ============================================================

simulate_checkerboard_bicop <- function(n, family, par = NULL, par2 = 3, d = 10) {
  cb <- checkerboard_density_bicop(
    family = family,
    par = par,
    par2 = par2,
    d = d
  )
  
  prob <- as.vector(cb$W)
  prob <- pmax(prob, 0)
  prob <- prob / sum(prob)
  
  cell_id <- sample(seq_along(prob), size = n, replace = TRUE, prob = prob)
  
  ij <- arrayInd(cell_id, .dim = c(d, d))
  
  U <- (ij[, 1] - 1 + runif(n)) / d
  V <- (ij[, 2] - 1 + runif(n)) / d
  
  cbind(U = U, V = V)
}


# ============================================================
# 6. Gradient of checkerboard density
# ============================================================

checkerboard_density_grad_bicop <- function(family, par = NULL, par2 = 3, d = 10) {
  cb <- checkerboard_density_bicop(
    family = family,
    par = par,
    par2 = par2,
    d = d
  )
  
  list(
    dW       = cb$dW,
    ddensity = cb$ddensity,
    family   = family,
    par      = par,
    par2     = par2,
    d        = d
  )
}


simulate_biv_mixpoisson <- function(
    n,
    mu1, phi1,
    mu2, phi2,
    dist1 = c("gamma", "lnorm", "pareto", "weibull", "ig"),
    dist2 = c("gamma", "lnorm", "pareto", "weibull", "ig"),
    family = 1,
    par = 0.0,
    par2 = 3,
    d = 10,
    seed = NULL,
    return_latent = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  
  dist1 <- match.arg(dist1)
  dist2 <- match.arg(dist2)
  
  # -------------------------------------------------------
  # Basic checks
  # -------------------------------------------------------
  if (!is.numeric(n) || length(n) != 1 || !is.finite(n) ||
      n <= 0 || abs(n - round(n)) > 1e-12) {
    stop("n must be a positive integer.")
  }
  n <- as.integer(round(n))
  
  if (!is.numeric(mu1) || length(mu1) != n ||
      any(!is.finite(mu1)) || any(mu1 <= 0)) {
    stop("mu1 must be a positive numeric vector of length n.")
  }
  
  if (!is.numeric(mu2) || length(mu2) != n ||
      any(!is.finite(mu2)) || any(mu2 <= 0)) {
    stop("mu2 must be a positive numeric vector of length n.")
  }
  
  if (!is.numeric(phi1) || length(phi1) != 1 ||
      !is.finite(phi1) || phi1 <= 0) {
    stop("phi1 must be a positive number.")
  }
  
  if (!is.numeric(phi2) || length(phi2) != 1 ||
      !is.finite(phi2) || phi2 <= 0) {
    stop("phi2 must be a positive number.")
  }
  
  if (!is.numeric(d) || length(d) != 1 || !is.finite(d) ||
      d <= 0 || abs(d - round(d)) > 1e-12) {
    stop("d must be a positive integer.")
  }
  d <- as.integer(round(d))
  
  # -------------------------------------------------------
  # Copula parameter checks
  # The family codes follow VineCopula convention:
  # 0  independence
  # 1  Gaussian
  # 2  Student-t, par = rho, par2 = nu fixed
  # 3  Clayton
  # 13 Clayton 180
  # 23 Clayton 90
  # 33 Clayton 270
  # 4  Gumbel
  # 14 Gumbel 180
  # 24 Gumbel 90
  # 34 Gumbel 270
  # 5  Frank
  # 25 Frank 90
  # -------------------------------------------------------
  allowed_family <- c(
    0, 1, 2,
    3, 13, 23, 33,
    4, 14, 24, 34,
    5, 25
  )
  
  if (!(family %in% allowed_family)) {
    stop("Unsupported copula family.")
  }
  
  if (family != 0) {
    if (is.null(par) || length(par) != 1 || !is.finite(par)) {
      stop("par must be a finite scalar.")
    }
  }
  
  if (family %in% c(1, 2)) {
    if (par <= -1 || par >= 1) {
      stop("Gaussian/Student: par must be rho in (-1,1).")
    }
  }
  
  if (family == 2) {
    if (is.null(par2) || length(par2) != 1 ||
        !is.finite(par2) || par2 <= 0) {
      stop("Student-t: par2 must be fixed degrees of freedom nu > 0.")
    }
  }
  
  if (family %in% c(3, 13, 23, 33)) {
    if (par <= 0) stop("Clayton: par must be positive.")
  }
  
  if (family %in% c(4, 14, 24, 34)) {
    if (par < 1) stop("Gumbel: par must be >= 1.")
  }
  
  if (family %in% c(5, 25)) {
    if (abs(par) < 1e-10) {
      stop("Frank: par must be nonzero; use family = 0 for independence.")
    }
  }
  
  # -------------------------------------------------------
  # Step 1: simulate latent uniforms from checkerboard copula
  # -------------------------------------------------------
  U <- simulate_checkerboard_bicop(
    n = n,
    family = family,
    par = par,
    par2 = par2,
    d = d
  )
  
  U1 <- U[, 1]
  U2 <- U[, 2]
  
  # Avoid exact numerical boundaries
  U1 <- pmin(1 - 1e-12, pmax(1e-12, U1))
  U2 <- pmin(1 - 1e-12, pmax(1e-12, U2))
  
  # -------------------------------------------------------
  # Helper: Weibull quantile under (mu, phi)
  # Q(u; mu, phi) = mu / Gamma(1+phi) * [-log(1-u)]^phi
  # -------------------------------------------------------
  qweibull_mu_phi_simu <- function(u, mu, phi) {
    c0 <- mu / gamma(1 + phi)
    c0 * (-log1p(-u))^phi
  }
  
  # -------------------------------------------------------
  # Step 2: latent variables
  # -------------------------------------------------------
  
  # ===== Margin 1 =====
  if (dist1 == "gamma") {
    shape1 <- 1 / phi1
    rate1  <- 1 / (phi1 * mu1)
    G1 <- qgamma(U1, shape = shape1, rate = rate1)
    
  } else if (dist1 == "lnorm") {
    sdlog1   <- log(1 + phi1)
    meanlog1 <- log(mu1) - 0.5 * sdlog1^2
    G1 <- qlnorm(U1, meanlog = meanlog1, sdlog = sdlog1)
    
  } else if (dist1 == "pareto") {
    G1 <- qpareto_mu_phi_simu(U1, mu = mu1, phi = phi1)
    
  } else if (dist1 == "weibull") {
    G1 <- qweibull_mu_phi_simu(U1, mu = mu1, phi = phi1)
    
  } else if (dist1 == "ig") {
    G1 <- statmod::qinvgauss(
      p = U1,
      mean = mu1,
      shape = 1 / phi1
    )
  }
  
  # ===== Margin 2 =====
  if (dist2 == "gamma") {
    shape2 <- 1 / phi2
    rate2  <- 1 / (phi2 * mu2)
    G2 <- qgamma(U2, shape = shape2, rate = rate2)
    
  } else if (dist2 == "lnorm") {
    sdlog2   <- log(1 + phi2)
    meanlog2 <- log(mu2) - 0.5 * sdlog2^2
    G2 <- qlnorm(U2, meanlog = meanlog2, sdlog = sdlog2)
    
  } else if (dist2 == "pareto") {
    G2 <- qpareto_mu_phi_simu(U2, mu = mu2, phi = phi2)
    
  } else if (dist2 == "weibull") {
    G2 <- qweibull_mu_phi_simu(U2, mu = mu2, phi = phi2)
    
  } else if (dist2 == "ig") {
    G2 <- statmod::qinvgauss(
      p = U2,
      mean = mu2,
      shape = 1 / phi2
    )
  }
  
  # -------------------------------------------------------
  # Step 3: mixed Poisson
  # -------------------------------------------------------
  X1 <- rpois(n, lambda = G1)
  X2 <- rpois(n, lambda = G2)
  
  out <- data.frame(X1 = X1, X2 = X2)
  
  if (return_latent) {
    out$G1 <- G1
    out$G2 <- G2
    out$U1 <- U1
    out$U2 <- U2
    out$mu1 <- mu1
    out$mu2 <- mu2
  }
  
  attr(out, "copula_family") <- family
  attr(out, "par") <- par
  attr(out, "par2") <- par2
  attr(out, "checkerboard_d") <- d
  attr(out, "latent_margins") <- list(
    dist1 = dist1, mu1 = mu1, phi1 = phi1,
    dist2 = dist2, mu2 = mu2, phi2 = phi2
  )
  
  out
}


simu_biv_bernoulli <- function(n, p1, p2, rho){
  
  if(abs(rho) >= 1){
    stop("rho must be in (-1,1)")
  }
  
  q1 <- qnorm(p1)
  q2 <- qnorm(p2)
  
  Z1 <- rnorm(n)
  Z2 <- rho * Z1 + sqrt(1 - rho^2) * rnorm(n)
  
  B1 <- as.integer(Z1 <= q1)
  B2 <- as.integer(Z2 <= q2)
  
  cbind(B1 = B1, B2 = B2)
}
