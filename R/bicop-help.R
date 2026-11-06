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


