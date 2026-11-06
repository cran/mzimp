# Cần package mvtnorm cho Gaussian copula
# install.packages("mvtnorm")

# =========================================================
# 1. CDF của từng copula
# 2. Đạo hàm theo theta của từng copula
# 3. Hàm tổng quát trả về 2 ma trận d x d
# =========================================================

one_par <- c(
  0,
  1,
  3, 13, 23, 33,
  4, 14, 24, 34,
  5,
  6, 16, 26, 36
)

two_par <- c(
  2,
  7, 17, 27, 37,
  8, 18, 28, 38,
  9, 19, 29, 39,
  10, 20, 30, 40,
  104, 114, 124, 134,
  204, 214, 224, 234
)


copula.list <- data.frame( copula.famity = c(one_par, two_par), number.par = c(rep(1, length(one_par)), rep(2, length(two_par))))

gaussian_copula_cdf_vec <- function(u, v, theta) {
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta) ||
      theta <= -1 || theta >= 1) {
    stop("For Gaussian copula, theta must be a single number in (-1, 1).")
  }
  
  if (!is.numeric(u) || !is.numeric(v)) {
    stop("u and v must be numeric vectors.")
  }
  
  if (length(u) != length(v)) {
    stop("u and v must have the same length.")
  }
  
  n <- length(u)
  out <- numeric(n)
  
  # các vị trí không hữu hạn
  bad <- !is.finite(u) | !is.finite(v)
  out[bad] <- NA_real_
  
  # boundary conditions
  idx <- !bad
  
  out[idx & (u <= 0 | v <= 0)] <- 0
  out[idx & (u >= 1 & v >= 1)] <- 1
  out[idx & (u >= 1 & v < 1)]  <- v[idx & (u >= 1 & v < 1)]
  out[idx & (v >= 1 & u < 1)]  <- u[idx & (v >= 1 & u < 1)]
  
  # các điểm nằm trong (0,1)^2
  interior <- idx & (u > 0 & u < 1) & (v > 0 & v < 1)
  
  if (any(interior)) {
    x <- qnorm(u[interior])
    y <- qnorm(v[interior])
    
    Sigma <- matrix(c(1, theta, theta, 1), 2, 2)
    
    out[interior] <- vapply(
      seq_along(x),
      function(i) {
        mvtnorm::pmvnorm(
          lower = c(-Inf, -Inf),
          upper = c(x[i], y[i]),
          mean  = c(0, 0),
          sigma = Sigma
        )[1]
      },
      numeric(1)
    )
  }
  
  return(out)
}

gaussian_copula_dtheta_vec <- function(u, v, theta, log = FALSE) {
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta) ||
      theta <= -1 || theta >= 1) {
    stop("For Gaussian copula, theta must be a single number in (-1, 1).")
  }
  
  if (!is.numeric(u) || !is.numeric(v)) {
    stop("u and v must be numeric vectors.")
  }
  
  if (length(u) != length(v)) {
    stop("u and v must have the same length.")
  }
  
  n <- length(u)
  out <- numeric(n)
  
  # non-finite positions
  bad <- !is.finite(u) | !is.finite(v)
  out[bad] <- NA_real_
  
  # boundary positions: derivative = 0
  idx <- !bad
  boundary <- idx & ((u <= 0) | (v <= 0) | (u >= 1) | (v >= 1))
  out[boundary] <- if (log) -Inf else 0
  
  # interior positions
  interior <- idx & (u > 0 & u < 1) & (v > 0 & v < 1)
  
  if (any(interior)) {
    x <- qnorm(u[interior])
    y <- qnorm(v[interior])
    
    # log density of bivariate normal with correlation theta
    logdet <- log1p(-theta^2)              # log(1 - theta^2)
    quad <- (x^2 - 2 * theta * x * y + y^2) / (1 - theta^2)
    
    logdens <- -log(2 * pi) - 0.5 * logdet - 0.5 * quad
    
    out[interior] <- if (log) logdens else exp(logdens)
  }
  
  return(out)
}

# ---------- Clayton copula ----------
clayton_copula_cdf <- function(u, v, theta) {
  if (theta <= -1 || theta == 0) {
    stop("For Clayton copula, theta must be > -1 and != 0.")
  }
  
  if (u <= 0 || v <= 0) return(0)
  if (u >= 1 && v >= 1) return(1)
  if (u >= 1) return(v)
  if (v >= 1) return(u)
  
  A <- u^(-theta) + v^(-theta) - 1
  
  if (theta > 0) {
    return(A^(-1 / theta))
  } else {
    # for -1 < theta < 0, use max(.,0)
    return(max(A, 0)^(-1 / theta))
  }
}

clayton_copula_dtheta <- function(u, v, theta) {
  if (theta <= -1 || theta == 0) {
    stop("For Clayton copula, theta must be > -1 and != 0.")
  }
  
  if (u <= 0 || v <= 0 || u >= 1 || v >= 1) return(0)
  
  A <- u^(-theta) + v^(-theta) - 1
  
  if (theta < 0 && A <= 0) return(0)
  
  Cuv <- A^(-1 / theta)
  
  # d/dtheta C = C * [ log(A)/theta^2 + (log(u)u^{-theta}+log(v)v^{-theta})/(theta A) ]
  term <- log(A) / theta^2 +
    (log(u) * u^(-theta) + log(v) * v^(-theta)) / (theta * A)
  
  Cuv * term
}

# ---------- Frank copula ----------
frank_copula_cdf <- function(u, v, theta) {
  if (theta == 0) {
    stop("For Frank copula, theta must be != 0.")
  }
  
  if (u <= 0 || v <= 0) return(0)
  if (u >= 1 && v >= 1) return(1)
  if (u >= 1) return(v)
  if (v >= 1) return(u)
  
  num <- (exp(-theta * u) - 1) * (exp(-theta * v) - 1)
  den <- exp(-theta) - 1
  
  -(1 / theta) * log(1 + num / den)
}

frank_copula_cdf_vec <- function(u, v, theta) {
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta) || theta == 0) {
    stop("For Frank copula, theta must be a single finite number different from 0.")
  }
  
  if (!is.numeric(u) || !is.numeric(v)) {
    stop("u and v must be numeric vectors.")
  }
  
  if (length(u) != length(v)) {
    stop("u and v must have the same length.")
  }
  
  n <- length(u)
  out <- numeric(n)
  
  bad <- !is.finite(u) | !is.finite(v)
  out[bad] <- NA_real_
  
  idx <- !bad
  
  # boundary conditions
  out[idx & (u <= 0 | v <= 0)] <- 0
  out[idx & (u >= 1 & v >= 1)] <- 1
  out[idx & (u >= 1 & v < 1)]  <- v[idx & (u >= 1 & v < 1)]
  out[idx & (v >= 1 & u < 1)]  <- u[idx & (v >= 1 & u < 1)]
  
  # interior points
  interior <- idx & (u > 0 & u < 1) & (v > 0 & v < 1)
  
  if (any(interior)) {
    ui <- u[interior]
    vi <- v[interior]
    
    num <- (exp(-theta * ui) - 1) * (exp(-theta * vi) - 1)
    den <- exp(-theta) - 1
    
    out[interior] <- -(1 / theta) * log(1 + num / den)
  }
  
  return(out)
}


frank_copula_cdf_vec_all <- function(u, v, theta, eps = 1e-10) {
  # kiểm tra input
  if (!is.numeric(u) || !is.numeric(v) || !is.numeric(theta)) {
    stop("u, v, theta must be numeric.")
  }
  
  if (!(length(u) == length(v) && length(v) == length(theta))) {
    stop("u, v, theta must have the same length.")
  }
  
  n <- length(u)
  out <- numeric(n)
  
  bad <- !is.finite(u) | !is.finite(v) | !is.finite(theta)
  out[bad] <- NA_real_
  
  idx <- !bad
  
  # boundary conditions
  out[idx & (u <= 0 | v <= 0)] <- 0
  out[idx & (u >= 1 & v >= 1)] <- 1
  out[idx & (u >= 1 & v < 1)]  <- v[idx & (u >= 1 & v < 1)]
  out[idx & (v >= 1 & u < 1)]  <- u[idx & (v >= 1 & u < 1)]
  
  # interior
  interior <- idx & (u > 0 & u < 1) & (v > 0 & v < 1)
  
  if (any(interior)) {
    ui <- u[interior]
    vi <- v[interior]
    ti <- theta[interior]
    
    # xử lý gần theta = 0 (độc lập)
    near0 <- abs(ti) < eps
    
    # case theta ≈ 0
    if (any(near0)) {
      out[interior][near0] <- ui[near0] * vi[near0]
    }
    
    # case theta ≠ 0
    if (any(!near0)) {
      t2 <- ti[!near0]
      u2 <- ui[!near0]
      v2 <- vi[!near0]
      
      num <- (exp(-t2 * u2) - 1) * (exp(-t2 * v2) - 1)
      den <- exp(-t2) - 1
      
      out[interior][!near0] <- -(1 / t2) * log(1 + num / den)
    }
  }
  
  return(out)
}

frank_copula_dtheta <- function(u, v, theta) {
  if (theta == 0) {
    stop("For Frank copula, theta must be != 0.")
  }
  
  if (u <= 0 || v <= 0 || u >= 1 || v >= 1) return(0)
  
  a <- exp(-theta * u)
  b <- exp(-theta * v)
  D <- exp(-theta) - 1
  B <- (a - 1) * (b - 1)
  A <- 1 + B / D
  
  # derivatives
  B_prime <- (-u * a) * (b - 1) + (a - 1) * (-v * b)
  D_prime <- -exp(-theta)
  A_prime <- (B_prime * D - B * D_prime) / D^2
  
  # C = -(1/theta) log(A)
  # dC/dtheta = log(A)/theta^2 - (1/theta)*(A'/A)
  log(A) / theta^2 - (A_prime / A) / theta
}
# ---------- Gumbel copula ----------
gumbel_copula_cdf <- function(u, v, theta) {
  if (theta < 1) {
    stop("For Gumbel copula, theta must be >= 1.")
  }
  
  if (u <= 0 || v <= 0) return(0)
  if (u >= 1 && v >= 1) return(1)
  if (u >= 1) return(v)
  if (v >= 1) return(u)
  
  x <- -log(u)
  y <- -log(v)
  A <- x^theta + y^theta
  
  exp(-A^(1 / theta))
}

gumbel_copula_dtheta <- function(u, v, theta) {
  if (theta < 1) {
    stop("For Gumbel copula, theta must be >= 1.")
  }
  
  if (u <= 0 || v <= 0 || u >= 1 || v >= 1) return(0)
  
  x <- -log(u)
  y <- -log(v)
  A <- x^theta + y^theta
  s <- A^(1 / theta)
  Cuv <- exp(-s)
  
  A_prime <- x^theta * log(x) + y^theta * log(y)
  
  # s' = s * [ -log(A)/theta^2 + (A'/A)/theta ]
  # dC/dtheta = -C * s'
  Cuv * s * (log(A) / theta^2 - (A_prime / A) / theta)
}

# =========================================================
# Hàm chọn family giống VineCopula
# =========================================================
copula_cdf_theta <- function(family, u, v, theta) {
  switch(
    as.character(family),
    "1" = gaussian_copula_cdf_vec(u, v, theta), # Gaussian
    "3" = clayton_copula_cdf(u, v, theta),  # Clayton
    "4" = gumbel_copula_cdf(u, v, theta),   # Gumbel
    "5" = frank_copula_cdf(u, v, theta),    # Frank
    stop("Unsupported family. Use 1=Gaussian, 3=Clayton, 4=Gumbel, 5=Frank.")
  )
}

copula_dtheta_theta <- function(family, u, v, theta) {
  switch(
    as.character(family),
    "1" = gaussian_copula_dtheta_vec(u, v, theta), # Gaussian
    "3" = clayton_copula_dtheta(u, v, theta),  # Clayton
    "4" = gumbel_copula_dtheta(u, v, theta),   # Gumbel
    "5" = frank_copula_dtheta(u, v, theta),    # Frank
    stop("Unsupported family. Use 1=Gaussian, 3=Clayton, 4=Gumbel, 5=Frank.")
  )
}

# =========================================================
# Hàm chính
# =========================================================
copula_rectangle_masses <- function(family, theta, d) {
  if (!is.numeric(d) || length(d) != 1 || d <= 0 || d != as.integer(d)) {
    stop("d must be a positive integer.")
  }
  
  grid <- (0:d) / d
  
  Cmat  <- matrix(0, nrow = d + 1, ncol = d + 1)
  dCmat <- matrix(0, nrow = d + 1, ncol = d + 1)
  
  for (i in 1:(d + 1)) {
    for (j in 1:(d + 1)) {
      u <- grid[i]
      v <- grid[j]
      Cmat[i, j]  <- copula_cdf_theta(family, u, v, theta)
      dCmat[i, j] <- copula_dtheta_theta(family, u, v, theta)
    }
  }
  
  mass <- matrix(0, nrow = d, ncol = d)
  d_mass <- matrix(0, nrow = d, ncol = d)
  
  for (i in 1:d) {
    for (j in 1:d) {
      mass[i, j] <-
        Cmat[i + 1, j + 1] + Cmat[i, j] - Cmat[i, j + 1] - Cmat[i + 1, j]
      
      d_mass[i, j] <-
        dCmat[i + 1, j + 1] + dCmat[i, j] - dCmat[i, j + 1] - dCmat[i + 1, j]
    }
  }
  
  list(
    mass = mass,
    d_mass = d_mass
  )
}

copula_rectangle_masses_tawn <- function(theta, delta, d) {
  if (!is.numeric(d) || length(d) != 1 || d <= 0 || d != as.integer(d)) {
    stop("d must be a positive integer.")
  }
  
  if (!is.numeric(theta) || length(theta) != 1 || !is.finite(theta) || theta < 1) {
    stop("theta must be a finite scalar >= 1.")
  }
  
  if (!is.numeric(delta) || length(delta) != 1 || !is.finite(delta) ||
      delta < 0 || delta > 1) {
    stop("delta must be a finite scalar in [0, 1].")
  }
  
  grid <- (0:d) / d
  
  Cmat      <- matrix(0, nrow = d + 1, ncol = d + 1)
  dthetaMat <- matrix(0, nrow = d + 1, ncol = d + 1)
  ddeltaMat <- matrix(0, nrow = d + 1, ncol = d + 1)
  
  for (i in 1:(d + 1)) {
    for (j in 1:(d + 1)) {
      u <- grid[i]
      v <- grid[j]
      
      Cmat[i, j]      <- tawn_copula_cdf(u, v, theta, delta)
      # dthetaMat[i, j] <- tawn_copula_dtheta(u, v, theta, delta)
      # ddeltaMat[i, j] <- tawn_copula_ddelta(u, v, theta, delta)
    }
  }
  
  mass         <- matrix(0, nrow = d, ncol = d)
  d_mass_theta <- matrix(0, nrow = d, ncol = d)
  d_mass_delta <- matrix(0, nrow = d, ncol = d)
  
  for (i in 1:d) {
    for (j in 1:d) {
      mass[i, j] <-
        Cmat[i + 1, j + 1] + Cmat[i, j] - Cmat[i, j + 1] - Cmat[i + 1, j]

      # d_mass_theta[i, j] <-
      #   dthetaMat[i + 1, j + 1] + dthetaMat[i, j] -
      #   dthetaMat[i, j + 1] - dthetaMat[i + 1, j]
      # 
      # d_mass_delta[i, j] <-
      #   ddeltaMat[i + 1, j + 1] + ddeltaMat[i, j] -
      #   ddeltaMat[i, j + 1] - ddeltaMat[i + 1, j]
    }
  }
  
  list(
    mass         = mass,
    d_mass_theta = d_mass_theta,
    d_mass_delta = d_mass_delta
  )
}

tawn_copula_cdf <- function(u, v, par, par2, type = 1, check.pars = TRUE,
                            eps = 1e-12) {
  if (!requireNamespace("VineCopula", quietly = TRUE)) {
    stop("Package 'VineCopula' is required.")
  }
  
  fam <- switch(as.character(type),
                "1" = 104L,
                "2" = 204L,
                stop("type must be 1 or 2."))
  
  if (!is.numeric(par) || length(par) != 1 || !is.finite(par)) {
    stop("par must be a finite scalar.")
  }
  if (!is.numeric(par2) || length(par2) != 1 || !is.finite(par2)) {
    stop("par2 must be a finite scalar.")
  }
  if (check.pars) {
    if (par < 1) stop("For Tawn copula, par must be >= 1.")
    if (par2 < 0 || par2 > 1) stop("For Tawn copula, par2 must be in [0, 1].")
  }
  
  ## ---- xử lý biên bằng tay ----
  if (u <= eps || v <= eps) return(0)
  if (u >= 1 - eps && v >= 1 - eps) return(1)
  if (u >= 1 - eps) return(v)
  if (v >= 1 - eps) return(u)
  
  out <- VineCopula::BiCopCDF(
    u1 = u, u2 = v,
    family = fam,
    par = par,
    par2 = par2,
    check.pars = check.pars
  )
  
  out <- as.numeric(out)
  
  if (!is.finite(out)) {
    stop(sprintf("BiCopCDF returned non-finite value at u=%g, v=%g, par=%g, par2=%g",
                 u, v, par, par2))
  }
  
  out
}
# =========================================================
# Helper: simulate (U1,U2) from chosen copula
# =========================================================

simulate_copula_u <- function(n, family = 1, theta = 0.5) {
  if (family == 1) {
    return(simulate_gaussian_copula_u(n, rho = theta))
  } else if (family == 3) {
    return(simulate_clayton_copula_u(n, theta = theta))
  } else if (family == 4) {
    return(simulate_gumbel_copula_u(n, theta = theta))
  } else if (family == 5) {
    return(simulate_frank_copula_u(n, theta = theta))
  } else {
    stop("Unsupported family. Use 1=Gaussian, 3=Clayton, 4=Gumbel, 5=Frank.")
  }
}

# =========================================================
# 1. Gaussian copula
# =========================================================
simulate_gaussian_copula_u <- function(n, rho) {
  Z1 <- rnorm(n)
  Z2 <- rnorm(n)
  
  Y1 <- Z1
  Y2 <- rho * Z1 + sqrt(1 - rho^2) * Z2
  
  cbind(pnorm(Y1), pnorm(Y2))
}

# =========================================================
# 2. Clayton copula
#
# Marshall-Olkin representation for Archimedean copula:
# if V ~ Gamma(shape = 1/theta, rate = 1)
# and E1,E2 ~ Exp(1), independent,
# then Uj = (1 + Ej/V)^(-1/theta)
# =========================================================
simulate_clayton_copula_u <- function(n, theta) {
  V  <- rgamma(n, shape = 1 / theta, rate = 1)
  E1 <- rexp(n, rate = 1)
  E2 <- rexp(n, rate = 1)
  
  U1 <- (1 + E1 / V)^(-1 / theta)
  U2 <- (1 + E2 / V)^(-1 / theta)
  
  cbind(U1, U2)
}

# =========================================================
# 3. Gumbel copula
#
# Uses positive stable mixing variable:
# If alpha = 1/theta in (0,1],
# S ~ positive stable(alpha) with Laplace transform
# E exp(-t S) = exp(-t^alpha),
# and E1,E2 ~ Exp(1),
# then Uj = exp(-(Ej / S)^alpha)
#
# Here we simulate positive stable using
# Chambers-Mallows-Stuck type formula.
# =========================================================
simulate_gumbel_copula_u <- function(n, theta) {
  alpha <- 1 / theta
  
  if (abs(alpha - 1) < 1e-14) {
    # independence case theta = 1
    return(cbind(runif(n), runif(n)))
  }
  
  S <- rstable_positive(n, alpha = alpha)
  E1 <- rexp(n, rate = 1)
  E2 <- rexp(n, rate = 1)
  
  U1 <- exp(-(E1 / S)^alpha)
  U2 <- exp(-(E2 / S)^alpha)
  
  cbind(U1, U2)
}

# =========================================================
# Positive stable random variable with LT exp(-t^alpha)
# for alpha in (0,1)
# =========================================================
rstable_positive <- function(n, alpha) {
  if (alpha <= 0 || alpha > 1) {
    stop("alpha must be in (0,1].")
  }
  
  if (abs(alpha - 1) < 1e-14) {
    return(rep(1, n))
  }
  
  U <- runif(n, 0, pi)
  W <- rexp(n, rate = 1)
  
  # Kanter representation
  X <- (sin(alpha * U) / (sin(U))^(alpha)) *
    (sin((1 - alpha) * U) / W)^((1 - alpha) / alpha)
  
  X
}

# =========================================================
# 4. Frank copula
#
# For simplicity and robustness, use conditional inversion:
# 1) U ~ Uniform(0,1)
# 2) W ~ Uniform(0,1)
# 3) solve C_{2|1}(v|u) = w for v
#
# Frank copula:
# C(u,v) = -(1/theta) log(1 + ((e^{-theta u}-1)(e^{-theta v}-1))/(e^{-theta}-1))
#
# Conditional cdf:
# d/du C(u,v)
# =========================================================
simulate_frank_copula_u <- function(n, theta) {
  U1 <- runif(n)
  W  <- runif(n)
  U2 <- numeric(n)
  
  for (i in seq_len(n)) {
    u <- U1[i]
    w <- W[i]
    
    froot <- function(v) frank_conditional_cdf(v = v, u = u, theta = theta) - w
    
    U2[i] <- uniroot(froot, interval = c(1e-12, 1 - 1e-12), tol = 1e-10)$root
  }
  
  cbind(U1, U2)
}

frank_conditional_cdf <- function(v, u, theta) {
  A <- exp(-theta * u)
  B <- exp(-theta * v)
  D <- exp(-theta) - 1
  
  # derivative of Frank copula wrt u
  num <- A * (B - 1)
  den <- D + (A - 1) * (B - 1)
  
  num / den
}

qpareto_mu_phi_simu <- function(p, mu, phi) {
  if (any(!is.finite(p)) || any(p < 0 | p > 1)) {
    stop("p must be in [0,1].")
  }
  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop("mu must be positive.")
  }
  if (!is.numeric(phi) || length(phi) != 1 || !is.finite(phi) || phi <= 0) {
    stop("phi must be a positive scalar.")
  }
  
  alpha  <- 2 + 1 / phi
  lambda <- mu * (1 + 1 / phi)
  
  lambda * ((1 - p)^(-1 / alpha) - 1)
}

simulate_biv_mixpoisson_v3 <- function(
    n,
    mu1, phi1,
    mu2, phi2,
    dist1 = c("gamma", "lnorm", "pareto", "weibull", "ig"),
    dist2 = c("gamma", "lnorm", "pareto", "weibull", "ig"),
    family = 1,
    theta = 0.0,
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
  
  if (!is.numeric(mu1) || length(mu1) != n || any(!is.finite(mu1)) || any(mu1 <= 0)) {
    stop("mu1 must be a positive numeric vector of length n.")
  }
  if (!is.numeric(mu2) || length(mu2) != n || any(!is.finite(mu2)) || any(mu2 <= 0)) {
    stop("mu2 must be a positive numeric vector of length n.")
  }
  
  if (!is.numeric(phi1) || length(phi1) != 1 || !is.finite(phi1) || phi1 <= 0) {
    stop("phi1 must be a positive number.")
  }
  if (!is.numeric(phi2) || length(phi2) != 1 || !is.finite(phi2) || phi2 <= 0) {
    stop("phi2 must be a positive number.")
  }
  
  # -------------------------------------------------------
  # Copula parameter checks
  # -------------------------------------------------------
  if (family == 1) {
    if (theta <= -1 || theta >= 1) stop("Gaussian: theta in (-1,1)")
  } else if (family == 3) {
    if (theta <= 0) stop("Clayton: theta > 0")
  } else if (family == 4) {
    if (theta < 1) stop("Gumbel: theta >= 1")
  } else if (family == 5) {
    if (theta == 0) stop("Frank: theta != 0")
  } else {
    stop("Unsupported family.")
  }
  
  # -------------------------------------------------------
  # Step 1: simulate copula
  # -------------------------------------------------------
  U <- simulate_copula_u(n = n, family = family, theta = theta)
  U1 <- U[, 1]
  U2 <- U[, 2]
  
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
    G1 <- statmod::qinvgauss(p = U1, mean = mu1, shape = 1/phi1)
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
    G2 <- statmod::qinvgauss(p = U2, mean = mu2, shape = 1/phi2)
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
  attr(out, "theta") <- theta
  attr(out, "latent_margins") <- list(
    dist1 = dist1, mu1 = mu1, phi1 = phi1,
    dist2 = dist2, mu2 = mu2, phi2 = phi2
  )
  
  out
}

# =========================================================
# Tổng quát cho rectangle masses của copula trong VineCopula
# - Hỗ trợ nhiều family theo mã của VineCopula
# - Tính mass trên các ô vuông [(i-1)/d, i/d] x [(j-1)/d, j/d]
# - Tính gradient theo từng tham số copula bằng finite difference
# =========================================================

copula_npars_vc <- function(family) {
  fam <- as.integer(family)
  
  one_par <- c(
    0,   # independence
    1,   # Gaussian
    3, 13, 23, 33,   # Clayton + rotations
    4, 14, 24, 34,   # Gumbel + rotations
    5,               # Frank
    6, 16, 26, 36    # Joe + rotations
  )
  
  two_par <- c(
    2,               # Student-t
    7, 17, 27, 37,   # BB1 + rotations
    8, 18, 28, 38,   # BB6 + rotations
    9, 19, 29, 39,   # BB7 + rotations
    10, 20, 30, 40,  # BB8 + rotations
    104, 114, 124, 134,  # Tawn type 1 + rotations
    204, 214, 224, 234   # Tawn type 2 + rotations
  )
  
  if (fam %in% one_par) return(1L)
  if (fam %in% two_par) return(2L)
  
  stop("Unsupported family code for this wrapper.")
}

copula_param_names_vc <- function(family) {
  k <- copula_npars_vc(family)
  if (k == 1L) return("par")
  c("par", "par2")
}

# ---------------------------------------------------------
# CDF trên vector điểm (u,v)
# theta:
#   - length 1 với họ 1 tham số
#   - length 2 với họ 2 tham số
# ---------------------------------------------------------
copula_cdf_vc_vec <- function(u, v, family, theta,
                              check.pars = TRUE,
                              eps = 1e-12) {
  if (!requireNamespace("VineCopula", quietly = TRUE)) {
    stop("Package 'VineCopula' is required.")
  }
  
  if (!is.numeric(u) || !is.numeric(v) || length(u) != length(v)) {
    stop("u and v must be numeric vectors of the same length.")
  }
  
  fam <- as.integer(family)
  k <- copula_npars_vc(fam)
  
  if (!is.numeric(theta) || length(theta) != k || any(!is.finite(theta))) {
    stop(sprintf("theta must be a numeric vector of length %d.", k))
  }
  
  n <- length(u)
  out <- numeric(n)
  
  bad <- !is.finite(u) | !is.finite(v)
  out[bad] <- NA_real_
  
  idx <- !bad
  if (!any(idx)) return(out)
  
  ui <- u[idx]
  vi <- v[idx]
  
  # Chuẩn hóa nhẹ để chống lỗi số kiểu 1 + 1e-16 hay -1e-16
  ui <- pmin(1, pmax(0, ui))
  vi <- pmin(1, pmax(0, vi))
  
  # --------- xử lý biên bằng công thức chính xác ----------
  res <- numeric(length(ui))
  
  at_u0 <- ui <= eps
  at_v0 <- vi <= eps
  at_u1 <- ui >= 1 - eps
  at_v1 <- vi >= 1 - eps
  
  # C(0,v)=0 hoặc C(u,0)=0
  res[at_u0 | at_v0] <- 0
  
  # C(1,1)=1
  both1 <- at_u1 & at_v1
  res[both1] <- 1
  
  # C(1,v)=v, với v<1
  only_u1 <- at_u1 & !at_v1 & !at_v0
  res[only_u1] <- vi[only_u1]
  
  # C(u,1)=u, với u<1
  only_v1 <- at_v1 & !at_u1 & !at_u0
  res[only_v1] <- ui[only_v1]
  
  # --------- Independence copula ----------
  interior <- !(at_u0 | at_v0 | at_u1 | at_v1)
  
  if (fam == 0L) {
    res[interior] <- ui[interior] * vi[interior]
    out[idx] <- res
    return(out)
  }
  
  # --------- chỉ gọi VineCopula ở miền trong (0,1)^2 ----------
  if (any(interior)) {
    par  <- theta[1]
    par2 <- if (k >= 2L) theta[2] else 0
    
    tmp <- VineCopula::BiCopCDF(
      u1 = ui[interior],
      u2 = vi[interior],
      family = fam,
      par = par,
      par2 = par2,
      check.pars = check.pars
    )
    
    res[interior] <- as.numeric(tmp)
  }
  
  out[idx] <- res
  out
}

# ---------------------------------------------------------
# Sinh lưới các đỉnh của các ô
# ---------------------------------------------------------
copula_grid_vertices_vc <- function(d) {
  if (!is.numeric(d) || length(d) != 1 || !is.finite(d) ||
      d <= 0 || abs(d - round(d)) > 1e-12) {
    stop("d must be a positive integer.")
  }
  
  grid <- (0:d) / d
  pts  <- expand.grid(u = grid, v = grid)
  
  list(
    grid = grid,
    u = pts$u,
    v = pts$v
  )
}

# ---------------------------------------------------------
# Chuyển ma trận CDF (d+1)x(d+1) -> mass dxd
# ---------------------------------------------------------
copula_cdfmat_to_mass_vc <- function(Cmat) {
  nr <- nrow(Cmat)
  nc <- ncol(Cmat)
  
  if (nr != nc || nr < 2L) {
    stop("Cmat must be a square matrix of size at least 2x2.")
  }
  
  d <- nr - 1L
  
  Cmat[2:(d + 1L), 2:(d + 1L), drop = FALSE] +
    Cmat[1:d, 1:d, drop = FALSE] -
    Cmat[1:d, 2:(d + 1L), drop = FALSE] -
    Cmat[2:(d + 1L), 1:d, drop = FALSE]
}

# ---------------------------------------------------------
# Tạo bước sai phân số theo từng tham số
# ---------------------------------------------------------
copula_fd_step_vc <- function(theta_j, eps = 1e-6) {
  eps * max(1, abs(theta_j))
}

# ---------------------------------------------------------
# Hàm thử tính CDF với một theta mới.
# Nếu lỗi hoặc không hữu hạn thì trả NULL.
# ---------------------------------------------------------
copula_try_cdfmat_vc <- function(family, theta, d, check.pars = TRUE) {
  gv <- copula_grid_vertices_vc(d)
  
  val <- try(
    copula_cdf_vc_vec(
      u = gv$u, v = gv$v,
      family = family,
      theta = theta,
      check.pars = check.pars
    ),
    silent = TRUE
  )
  
  if (inherits(val, "try-error")) return(NULL)
  if (any(!is.finite(val))) return(NULL)
  
  matrix(val, nrow = d + 1L, ncol = d + 1L, byrow = FALSE)
}

# ---------------------------------------------------------
# Gradient CDF theo từng tham số bằng sai phân số
# Trả về list các ma trận đạo hàm của CDF trên lưới
# ---------------------------------------------------------
copula_cdf_gradmats_vc <- function(family, theta, d,
                                   eps = 1e-6,
                                   check.pars = TRUE) {
  fam <- as.integer(family)
  k <- copula_npars_vc(fam)
  
  if (fam == 0L) {
    return(vector("list", 0L))
  }
  
  base_C <- copula_try_cdfmat_vc(
    family = fam, theta = theta, d = d, check.pars = check.pars
  )
  if (is.null(base_C)) {
    stop("Base copula CDF could not be evaluated.")
  }
  
  out <- vector("list", k)
  names(out) <- copula_param_names_vc(fam)
  
  for (j in seq_len(k)) {
    h <- copula_fd_step_vc(theta[j], eps = eps)
    
    th_plus  <- theta
    th_minus <- theta
    th_plus[j]  <- th_plus[j] + h
    th_minus[j] <- th_minus[j] - h
    
    C_plus <- copula_try_cdfmat_vc(
      family = fam, theta = th_plus, d = d, check.pars = check.pars
    )
    C_minus <- copula_try_cdfmat_vc(
      family = fam, theta = th_minus, d = d, check.pars = check.pars
    )
    
    # Ưu tiên central difference
    if (!is.null(C_plus) && !is.null(C_minus)) {
      out[[j]] <- (C_plus - C_minus) / (2 * h)
      next
    }
    
    # fallback: forward difference
    if (!is.null(C_plus)) {
      out[[j]] <- (C_plus - base_C) / h
      next
    }
    
    # fallback: backward difference
    if (!is.null(C_minus)) {
      out[[j]] <- (base_C - C_minus) / h
      next
    }
    
    stop(sprintf("Could not compute numerical derivative for parameter %d.", j))
  }
  
  out
}

# ---------------------------------------------------------
# Hàm chính
# - mass: ma trận d x d
# - d_mass:
#     + matrix d x d nếu copula có 1 tham số
#     + array d x d x k nếu copula có 2 tham số
# ---------------------------------------------------------
copula_rectangle_masses_vc <- function(family, theta, d,
                                       eps = 1e-6,
                                       check.pars = TRUE,
                                       clamp_negative = TRUE,
                                       renormalize = FALSE) {
  fam <- as.integer(family)
  k <- copula_npars_vc(fam)
  
  if (!is.numeric(theta) || length(theta) != k || any(!is.finite(theta))) {
    stop(sprintf("theta must be a numeric vector of length %d.", k))
  }
  
  Cmat <- copula_try_cdfmat_vc(
    family = fam, theta = theta, d = d, check.pars = check.pars
  )
  if (is.null(Cmat)) {
    stop("Could not evaluate copula CDF matrix at supplied parameter(s).")
  }
  
  mass <- copula_cdfmat_to_mass_vc(Cmat)
  
  if (clamp_negative) {
    mass[mass < 0 & mass > -1e-12] <- 0
  }
  
  if (renormalize) {
    s <- sum(mass)
    if (is.finite(s) && s > 0) mass <- mass / s
  }
  
  grad_C_list <- copula_cdf_gradmats_vc(
    family = fam,
    theta = theta,
    d = d,
    eps = eps,
    check.pars = check.pars
  )
  
  if (length(grad_C_list) == 0L) {
    d_mass <- NULL
  } else if (length(grad_C_list) == 1L) {
    d_mass <- copula_cdfmat_to_mass_vc(grad_C_list[[1]])
  } else {
    d_mass <- array(
      0,
      dim = c(d, d, length(grad_C_list)),
      dimnames = list(
        i = seq_len(d),
        j = seq_len(d),
        param = names(grad_C_list)
      )
    )
    for (m in seq_along(grad_C_list)) {
      d_mass[, , m] <- copula_cdfmat_to_mass_vc(grad_C_list[[m]])
    }
  }
  
  list(
    family = fam,
    theta = theta,
    mass = mass,
    d_mass = d_mass,
    Cmat = Cmat
  )
}

# =========================================================
# Map tham số tự do par -> theta hợp lệ cho VineCopula family
# - input:
#     par    : scalar hoặc vector tự do từ optimizer
#     family : mã family theo VineCopula
# - output:
#     scalar nếu family có 1 tham số
#     vector length 2 nếu family có 2 tham số
# =========================================================

copula_npars_vc <- function(family) {
  fam <- as.integer(family)
  
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
  
  if (fam %in% one_par) return(1L)
  if (fam %in% two_par) return(2L)
  
  stop("Unsupported family code.")
}

map_theta_copula <- function(par, family, eps = 1e-6) {
  fam <- as.integer(family)
  k <- copula_npars_vc(fam)
  
  if (!is.numeric(par) || any(!is.finite(par))) {
    stop("par must be numeric and finite.")
  }
  if (length(par) != k) {
    stop(sprintf("For family %s, par must have length %d.", fam, k))
  }
  
  # convenience
  softplus <- function(x) log1p(exp(-abs(x))) + pmax(x, 0)
  
  # independence
  if (fam == 0L) {
    return(numeric(0))
  }
  
  # ---------- 1-parameter families ----------
  if (k == 1L) {
    x <- par[1]
    
    # Gaussian: (-1, 1)
    if (fam == 1L) {
      return((1 - eps) * tanh(x))
    }
    
    # Clayton / survival Clayton: (0, Inf)
    if (fam %in% c(3, 13)) {
      return(softplus(x) + eps)
    }
    
    # Rotated Clayton 90/270: (-Inf, 0)
    if (fam %in% c(23, 33)) {
      return(-(softplus(x) + eps))
    }
    
    # Gumbel / survival Gumbel: [1, Inf)
    if (fam %in% c(4, 14)) {
      return(1 + softplus(x))
    }
    
    # Rotated Gumbel 90/270: (-Inf, -1]
    if (fam %in% c(24, 34)) {
      return(-(1 + softplus(x)))
    }
    
    # Frank: R \ {0}
    if (fam == 5L) {
      # tránh đúng 0
      return(ifelse(x >= 0, x + eps, x - eps))
    }
    
    # Joe / survival Joe: (1, Inf)
    if (fam %in% c(6, 16)) {
      return(1 + softplus(x) + eps)
    }
    
    # Rotated Joe 90/270: (-Inf, -1)
    if (fam %in% c(26, 36)) {
      return(-(1 + softplus(x) + eps))
    }
    
    if (fam == 2L) {
      rho <- (1 - eps) * tanh(x)
      return(rho)
    }
    
    stop("Unsupported 1-parameter family.")
  }
  
  
  
  # # ---------- 2-parameter families ----------
  # x1 <- par[1]
  # x2 <- par[2]
  # 
  # # Student t: par in (-1,1), par2 in (2, Inf)
  # 
  # 
  # # BB1 / survival BB1: par in (0,Inf), par2 in [1,Inf)
  # if (fam %in% c(7, 17)) {
  #   th1 <- softplus(x1) + eps
  #   th2 <- 1 + softplus(x2)
  #   return(c(th1, th2))
  # }
  # 
  # # Rotated BB1 90/270: par in (-Inf,0), par2 in (-Inf,-1]
  # if (fam %in% c(27, 37)) {
  #   th1 <- -(softplus(x1) + eps)
  #   th2 <- -(1 + softplus(x2))
  #   return(c(th1, th2))
  # }
  # 
  # # BB6 / survival BB6: par in [1,Inf), par2 in [1,Inf)
  # if (fam %in% c(8, 18)) {
  #   th1 <- 1 + softplus(x1)
  #   th2 <- 1 + softplus(x2)
  #   return(c(th1, th2))
  # }
  # 
  # # Rotated BB6 90/270: par in (-Inf,-1], par2 in (-Inf,-1]
  # if (fam %in% c(28, 38)) {
  #   th1 <- -(1 + softplus(x1))
  #   th2 <- -(1 + softplus(x2))
  #   return(c(th1, th2))
  # }
  # 
  # # BB7 / survival BB7: par in [1,Inf), par2 in (0,Inf)
  # if (fam %in% c(9, 19)) {
  #   th1 <- 1 + softplus(x1)
  #   th2 <- softplus(x2) + eps
  #   return(c(th1, th2))
  # }
  # 
  # # Rotated BB7 90/270: par in (-Inf,-1], par2 in (-Inf,0)
  # if (fam %in% c(29, 39)) {
  #   th1 <- -(1 + softplus(x1))
  #   th2 <- -(softplus(x2) + eps)
  #   return(c(th1, th2))
  # }
  # 
  # # BB8 / survival BB8: par in [1,Inf), par2 in (0,1]
  # if (fam %in% c(10, 20)) {
  #   th1 <- 1 + softplus(x1)
  #   th2 <- 1 / (1 + exp(-x2))          # logistic in (0,1)
  #   th2 <- pmin(1, pmax(eps, th2))
  #   return(c(th1, th2))
  # }
  # 
  # # Rotated BB8 90/270: par in (-Inf,-1], par2 in [-1,0)
  # if (fam %in% c(30, 40)) {
  #   th1 <- -(1 + softplus(x1))
  #   th2 <- -1 / (1 + exp(-x2))         # in (-1,0)
  #   th2 <- pmin(-eps, pmax(-1, th2))
  #   return(c(th1, th2))
  # }
  # 
  # # Tawn type 1/2, survival versions:
  # # par in [1,Inf), par2 in [0,1]
  # if (fam %in% c(104, 114, 204, 214)) {
  #   th1 <- 1 + softplus(x1)
  #   th2 <- 1 / (1 + exp(-x2))          # in (0,1)
  #   th2 <- pmin(1, pmax(0, th2))
  #   return(c(th1, th2))
  # }
  # 
  # # Rotated Tawn type 1/2:
  # # par in (-Inf,-1], par2 in [0,1]
  # if (fam %in% c(124, 134, 224, 234)) {
  #   th1 <- -(1 + softplus(x1))
  #   th2 <- 1 / (1 + exp(-x2))
  #   th2 <- pmin(1, pmax(0, th2))
  #   return(c(th1, th2))
  # }
  # 
  # stop("Unsupported 2-parameter family.")
}

init_par_copula_indep <- function(family) {
  fam <- as.integer(family)
  
  # 1 tham số
  if (fam == 1L) return(0)        # Gaussian: tanh(0)=0
  if (fam == 5L) return(0)        # Frank: trực tiếp
  
  if (fam %in% c(3,13)) return(-10)   # Clayton → softplus(-10) ~ 0
  if (fam %in% c(23,33)) return(10)   # rotated Clayton → ~0-
  
  if (fam %in% c(4,14)) return(-10)   # Gumbel → 1 + small ≈ 1
  if (fam %in% c(24,34)) return(10)
  
  if (fam %in% c(6,16)) return(-10)   # Joe → ≈1
  if (fam %in% c(26,36)) return(10)
  
  # 2 tham số
  if (fam == 2L) return(c(0, 0))      # t: rho=0, df ~ 2+
  
  if (fam %in% c(7,17)) return(c(-10, -10))  # BB1 → (0,1)
  if (fam %in% c(27,37)) return(c(10, 10))
  
  if (fam %in% c(8,18)) return(c(-10, -10))  # BB6 → (1,1)
  if (fam %in% c(28,38)) return(c(10, 10))
  
  if (fam %in% c(9,19)) return(c(-10, -10))  # BB7 → (1,0)
  if (fam %in% c(29,39)) return(c(10, 10))
  
  if (fam %in% c(10,20)) return(c(-10, 0))   # BB8 → (1,0.5)
  if (fam %in% c(30,40)) return(c(10, 0))
  
  if (fam %in% c(104,114,204,214)) return(c(-10, 0)) # Tawn
  if (fam %in% c(124,134,224,234)) return(c(10, 0))
  
  stop("Unsupported family")
}


############### FACTOR COPULA #####################

make_gh_grid <- function(J, gh_order = 15) {
  gh <- statmod::gauss.quad(gh_order, kind = "hermite")
  idx <- expand.grid(rep(list(seq_len(gh_order)), J))
  
  z <- sqrt(2) * as.matrix(
    sapply(seq_len(J), function(r) gh$nodes[idx[[r]]])
  )
  
  w <- apply(
    as.matrix(sapply(seq_len(J), function(r) gh$weights[idx[[r]]])),
    1, prod
  ) / (pi^(J / 2))
  
  list(z = z, w = w)
}


q_dq_one <- function(u_grid, psi_j, z) {
  J <- length(psi_j)
  d1 <- length(u_grid)
  
  qval <- numeric(d1)
  dq <- matrix(0, d1, J)
  
  qval[1] <- 0
  qval[d1] <- 1
  
  r2 <- sum(psi_j^2)
  if (r2 >= 1) stop("Each row of psi must satisfy ||psi_j|| < 1.")
  
  s <- sqrt(1 - r2)
  eta <- sum(psi_j * z)
  
  ind <- 2:(d1 - 1)
  tval <- qnorm(u_grid[ind])
  
  A <- tval - eta
  aval <- A / s
  
  qval[ind] <- pnorm(aval)
  
  phi_aval <- dnorm(aval)
  
  for (r in seq_len(J)) {
    da <- -z[r] / s + A * psi_j[r] / (s^3)
    dq[ind, r] <- phi_aval * da
  }
  
  list(q = qval, dq = dq)
}


negloglik_gauss_factor <- function(psi, Iarr, gh_grid, eps = 1e-300) {
  n <- dim(Iarr)[1]
  k <- dim(Iarr)[2]
  d <- dim(Iarr)[3]
  J <- ncol(psi)
  
  u_grid <- seq(0, 1, length.out = d + 1)
  
  Li <- numeric(n)
  
  for (ell in seq_along(gh_grid$w)) {
    z <- gh_grid$z[ell, ]
    weight <- gh_grid$w[ell]
    
    P <- matrix(1, n, k)
    
    for (j in seq_len(k)) {
      tmp <- q_dq_one(u_grid, psi[j, ], z)
      qdiff <- diff(tmp$q)
      
      P[, j] <- Iarr[, j, ] %*% qdiff
    }
    
    Li <- Li + weight * apply(P, 1, prod)
  }
  
  -sum(log(pmax(Li, eps)))
}


grad_negloglik_gauss_factor <- function(psi, Iarr, gh_grid, eps = 1e-300) {
  n <- dim(Iarr)[1]
  k <- dim(Iarr)[2]
  d <- dim(Iarr)[3]
  J <- ncol(psi)
  
  u_grid <- seq(0, 1, length.out = d + 1)
  
  Li <- numeric(n)
  dLi <- array(0, dim = c(n, k, J))
  
  for (ell in seq_along(gh_grid$w)) {
    z <- gh_grid$z[ell, ]
    weight <- gh_grid$w[ell]
    
    P <- matrix(1, n, k)
    DP <- array(0, dim = c(n, k, J))
    
    for (j in seq_len(k)) {
      tmp <- q_dq_one(u_grid, psi[j, ], z)
      
      qdiff <- diff(tmp$q)
      P[, j] <- Iarr[, j, ] %*% qdiff
      
      for (r in seq_len(J)) {
        dqdiff <- diff(tmp$dq[, r])
        DP[, j, r] <- Iarr[, j, ] %*% dqdiff
      }
    }
    
    prefix <- matrix(1, n, k)
    suffix <- matrix(1, n, k)
    
    if (k > 1) {
      for (j in 2:k) {
        prefix[, j] <- prefix[, j - 1] * P[, j - 1]
      }
      for (j in (k - 1):1) {
        suffix[, j] <- suffix[, j + 1] * P[, j + 1]
      }
    }
    
    prod_all <- apply(P, 1, prod)
    Li <- Li + weight * prod_all
    
    for (j in seq_len(k)) {
      prod_except_j <- prefix[, j] * suffix[, j]
      
      for (r in seq_len(J)) {
        dLi[, j, r] <- dLi[, j, r] +
          weight * prod_except_j * DP[, j, r]
      }
    }
  }
  
  Li_safe <- pmax(Li, eps)
  
  grad <- matrix(0, k, J)
  
  for (j in seq_len(k)) {
    for (r in seq_len(J)) {
      grad[j, r] <- -sum(dLi[, j, r] / Li_safe)
    }
  }
  
  grad
}


simulate_factor_checkerboard_gaussian <- function(n, psi, d, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  psi <- as.matrix(psi)
  
  k <- nrow(psi)
  J <- ncol(psi)
  
  norm2 <- rowSums(psi^2)
  
  if (any(norm2 >= 1)) {
    stop("Each row of psi must satisfy ||psi_j||^2 < 1.")
  }
  
  # Common factors: n x J
  Z <- matrix(rnorm(n * J), nrow = n, ncol = J)
  
  # Idiosyncratic noises: n x k
  Eps <- matrix(rnorm(n * k), nrow = n, ncol = k)
  
  # Gaussian latent variables G_j
  G <- Z %*% t(psi) +
    sweep(Eps, 2, sqrt(1 - norm2), "*")
  
  # Uniform variables from Gaussian copula
  U <- pnorm(G)
  
  # Checkerboard cell indices: 1,...,d
  M <- ceiling(d * U)
  M[M < 1] <- 1
  M[M > d] <- d
  
  colnames(U) <- paste0("U", seq_len(k))
  colnames(M) <- paste0("M", seq_len(k))
  
  list(
    U = U,     # simulated Gaussian copula uniforms
    M = M,     # checkerboard cells
    Z = Z,     # common factors
    G = G      # latent Gaussian variables
  )
}


simulate_k_mixpoisson_factor_checkerboard <- function(
    n,
    mu,
    phi,
    dist,
    psi,
    d = 10,
    seed = NULL,
    return_latent = TRUE
) {
  if (!is.null(seed)) set.seed(seed)
  
  psi <- as.matrix(psi)
  k <- nrow(psi)
  J <- ncol(psi)
  
  if (length(dist) != k) {
    stop("dist must be a character vector of length k.")
  }
  
  dist <- match.arg(
    dist,
    choices = c("gamma", "lnorm", "pareto", "weibull", "ig"),
    several.ok = TRUE
  )
  
  if (length(phi) != k || any(!is.finite(phi)) || any(phi <= 0)) {
    stop("phi must be a positive numeric vector of length k.")
  }
  
  if (is.vector(mu)) {
    if (k != 1 || length(mu) != n) {
      stop("If mu is a vector, then k must be 1 and length(mu) = n.")
    }
    mu <- matrix(mu, nrow = n, ncol = 1)
  } else {
    mu <- as.matrix(mu)
  }
  
  if (nrow(mu) != n || ncol(mu) != k) {
    stop("mu must be an n x k matrix.")
  }
  
  if (any(!is.finite(mu)) || any(mu <= 0)) {
    stop("All entries of mu must be positive and finite.")
  }
  
  norm2 <- rowSums(psi^2)
  if (any(norm2 >= 1)) {
    stop("Each row of psi must satisfy ||psi_j||^2 < 1.")
  }
  
  # -------------------------------------------------------
  # Step 1: simulate factor Gaussian copula
  # -------------------------------------------------------
  Z <- matrix(rnorm(n * J), nrow = n, ncol = J)
  Eps <- matrix(rnorm(n * k), nrow = n, ncol = k)
  
  G_gauss <- Z %*% t(psi) +
    sweep(Eps, 2, sqrt(1 - norm2), "*")
  
  U_gauss <- pnorm(G_gauss)
  
  # -------------------------------------------------------
  # Step 2: convert Gaussian copula to checkerboard copula
  # M_j is the checkerboard cell index
  # U_cb is uniform inside the selected cell
  # -------------------------------------------------------
  M <- ceiling(d * U_gauss)
  M[M < 1] <- 1
  M[M > d] <- d
  
  V <- matrix(runif(n * k), nrow = n, ncol = k)
  
  U_cb <- (M - 1 + V) / d
  # U_cb <- pmin(1 - 1e-12, pmax(1e-12, U_cb))
  
  # -------------------------------------------------------
  # Helper: quantile functions for latent Lambda
  # -------------------------------------------------------
  qweibull_mu_phi_simu <- function(u, mu, phi) {
    c0 <- mu / gamma(1 + phi)
    c0 * (-log1p(-u))^phi
  }
  
  qlatent_one_margin <- function(u, mu_j, phi_j, dist_j) {
    if (dist_j == "gamma") {
      shape <- 1 / phi_j
      rate  <- 1 / (phi_j * mu_j)
      out <- qgamma(u, shape = shape, rate = rate)
      
    } else if (dist_j == "lnorm") {
      sdlog   <- sqrt(log(1 + phi_j))
      meanlog <- log(mu_j) - 0.5 * sdlog^2
      out <- qlnorm(u, meanlog = meanlog, sdlog = sdlog)
      
    } else if (dist_j == "weibull") {
      out <- qweibull_mu_phi_simu(u, mu = mu_j, phi = phi_j)
      
    } else if (dist_j == "ig") {
      out <- statmod::qinvgauss(
        p = u,
        mean = mu_j,
        shape = 1 / phi_j
      )
      
    } else if (dist_j == "pareto") {
      if (!exists("qpareto_mu_phi_simu")) {
        stop("Please define qpareto_mu_phi_simu(u, mu, phi) before using dist = 'pareto'.")
      }
      out <- qpareto_mu_phi_simu(u, mu = mu_j, phi = phi_j)
    }
    
    out
  }
  
  # -------------------------------------------------------
  # Step 3: simulate latent Lambda_j
  # -------------------------------------------------------
  Lambda <- matrix(NA_real_, nrow = n, ncol = k)
  
  for (j in seq_len(k)) {
    Lambda[, j] <- qlatent_one_margin(
      u = U_cb[, j],
      mu_j = mu[, j],
      phi_j = phi[j],
      dist_j = dist[j]
    )
  }
  
  # -------------------------------------------------------
  # Step 4: mixed Poisson
  # -------------------------------------------------------
  X <- matrix(rpois(n * k, lambda = as.vector(Lambda)), nrow = n, ncol = k)
  
  colnames(X) <- paste0("X", seq_len(k))
  colnames(Lambda) <- paste0("Lambda", seq_len(k))
  colnames(U_cb) <- paste0("U", seq_len(k))
  colnames(U_gauss) <- paste0("U_gauss", seq_len(k))
  colnames(M) <- paste0("M", seq_len(k))
  
  out <- as.data.frame(X)
  
  if (return_latent) {
    out <- cbind(
      out,
      as.data.frame(Lambda),
      as.data.frame(U_cb),
      as.data.frame(M)
    )
  }
  
  attr(out, "dist") <- dist
  attr(out, "phi") <- phi
  attr(out, "psi") <- psi
  attr(out, "checkerboard_d") <- d
  attr(out, "J") <- J
  
  out
}

make_psi_from_psi0 <- function(psi0, k, J, return_all = FALSE) {
  
  npar <- k * J - J * (J - 1) / 2
  
  if (length(psi0) != npar) {
    stop(paste0("psi0 must have length ", npar))
  }
  
  A <- matrix(0, k, J)
  map <- matrix(NA_integer_, npar, 2)
  
  id <- 1
  
  for (j in seq_len(k)) {
    for (r in seq_len(J)) {
      if (r <= j) {
        
        if (j == r) {
          A[j, r] <- exp(psi0[id])
        } else {
          A[j, r] <- psi0[id]
        }
        
        map[id, ] <- c(j, r)
        id <- id + 1
      }
    }
  }
  
  denom <- sqrt(1 + rowSums(A^2))
  Psi <- sweep(A, 1, denom, "/")
  
  if (return_all) {
    return(list(Psi = Psi, A = A, map = map))
  }
  
  Psi
}

simulate_k_mixpoisson_biexp_factor_checkerboard <- function(
    n,
    mu,
    phi,
    dist,
    psi,
    lambda,
    d = 10,
    seed = NULL,
    return_latent = TRUE,
    zero_tol = 1e-12,
    sep_tol = 1e-9
) {

  if (!is.null(seed)) set.seed(seed)

  psi <- as.matrix(psi)

  k <- nrow(psi)
  J <- ncol(psi)

  ## ============================================================
  ## 0. Input checks
  ## ============================================================

  if (length(dist) != k) {
    stop("dist must be a character vector of length k.")
  }

  dist <- match.arg(
    dist,
    choices = c(
      "gamma",
      "lnorm",
      "pareto",
      "weibull",
      "ig"
    ),
    several.ok = TRUE
  )

  if (length(phi) != k ||
      any(!is.finite(phi)) ||
      any(phi <= 0)) {
    stop("phi must be a positive numeric vector of length k.")
  }

  if (length(lambda) != 1L ||
      !is.finite(lambda) ||
      lambda <= 0 ||
      lambda >= 1) {
    stop("lambda must be in (0,1).")
  }

  if (is.vector(mu)) {

    if (k != 1 || length(mu) != n) {
      stop(
        "If mu is a vector, k must be 1 and length(mu) = n."
      )
    }

    mu <- matrix(
      mu,
      nrow = n,
      ncol = 1
    )

  } else {

    mu <- as.matrix(mu)
  }

  if (nrow(mu) != n ||
      ncol(mu) != k) {
    stop("mu must be an n x k matrix.")
  }

  if (any(!is.finite(mu)) ||
      any(mu <= 0)) {
    stop("All entries of mu must be positive and finite.")
  }

  norm2 <- rowSums(psi^2)

  if (any(norm2 >= 1)) {
    stop(
      "Each row of psi must satisfy ||psi_j||^2 < 1."
    )
  }

  sigma <- sqrt(1 - norm2)


  ## ============================================================
  ## 1. Standard bilateral exponential
  ##
  ## Z =
  ## a E1 - b E2 - (a-b)
  ##
  ## a = lambda
  ## b = sqrt(1-lambda^2)
  ##
  ## E(Z)=0, Var(Z)=1
  ## ============================================================

  a0 <- lambda
  b0 <- sqrt(1 - lambda^2)

  shift0 <- a0 - b0

  z_break <- -shift0

  p_break <- b0 / (a0 + b0)


  ## ------------------------------------------------------------
  ## RNG
  ## ------------------------------------------------------------

  r_sd_biexp <- function(n) {

    E1 <- rexp(n)
    E2 <- rexp(n)

    a0 * E1 -
      b0 * E2 -
      shift0
  }


  ## ------------------------------------------------------------
  ## Exact CDF of SdBiExp
  ## ------------------------------------------------------------

  F_sd_biexp <- function(z) {

    zz <- as.numeric(z)

    out <- numeric(length(zz))

    id_left <- zz < z_break
    id_right <- !id_left

    if (any(id_left)) {

      out[id_left] <-
        b0 / (a0 + b0) *
        exp(
          (zz[id_left] + shift0) / b0
        )
    }

    if (any(id_right)) {

      out[id_right] <-
        1 -
        a0 / (a0 + b0) *
        exp(
          -(zz[id_right] + shift0) / a0
        )
    }

    pmin(
      pmax(out, 0),
      1
    )
  }


  ## ------------------------------------------------------------
  ## Exact density of SdBiExp
  ## ------------------------------------------------------------

  f_sd_biexp <- function(z) {

    zz <- as.numeric(z)

    out <- numeric(length(zz))

    id_left <- zz < z_break
    id_right <- !id_left

    if (any(id_left)) {

      out[id_left] <-
        1 / (a0 + b0) *
        exp(
          (zz[id_left] + shift0) / b0
        )
    }

    if (any(id_right)) {

      out[id_right] <-
        1 / (a0 + b0) *
        exp(
          -(zz[id_right] + shift0) / a0
        )
    }

    out
  }


  ## ------------------------------------------------------------
  ## Exact quantile of SdBiExp
  ## ------------------------------------------------------------

  Q_sd_biexp <- function(u) {

    uu <- as.numeric(u)

    uu <- pmin(
      pmax(uu, 1e-15),
      1 - 1e-15
    )

    out <- numeric(length(uu))

    id_left <- uu < p_break
    id_right <- !id_left

    if (any(id_left)) {

      out[id_left] <-
        b0 *
        log(
          (a0 + b0) *
            uu[id_left] /
            b0
        ) -
        shift0
    }

    if (any(id_right)) {

      out[id_right] <-
        -a0 *
        log(
          (a0 + b0) *
            (1 - uu[id_right]) /
            a0
        ) -
        shift0
    }

    out
  }


  ## ============================================================
  ## 2. Partial-fraction coefficients
  ##
  ## A_i =
  ## prod_{h != i}
  ## alpha_i / (alpha_i - alpha_h)
  ##
  ## No derivative is needed for simulation.
  ## ============================================================

  get_pf_coef <- function(
      scale,
      sep_tol = 1e-9
  ) {

    L <- length(scale)

    if (L == 0L) {
      stop("No active transformed scales.")
    }

    if (any(!is.finite(scale)) ||
        any(scale <= 0)) {
      stop("All transformed scales must be positive.")
    }

    if (L == 1L) {
      return(1)
    }

    DD <- abs(
      outer(
        scale,
        scale,
        "-"
      )
    )

    diag(DD) <- Inf

    relative_sep <-
      DD /
      pmax(
        outer(
          scale,
          scale,
          pmax
        ),
        1e-300
      )

    if (min(relative_sep) < sep_tol) {

      stop(
        paste0(
          "Nearly repeated transformed scales detected. ",
          "The simple partial-fraction representation is ",
          "numerically unstable."
        )
      )
    }

    coef <- numeric(L)

    for (i in seq_len(L)) {

      coef[i] <-
        prod(
          scale[i] /
            (
              scale[i] -
                scale[-i]
            )
        )
    }

    coef
  }


  ## ============================================================
  ## 3. Exact distribution of
  ##
  ## S = sum_l v_l Z_l
  ##
  ## where Z_l iid SdBiExp(lambda)
  ##
  ## For margin j:
  ##
  ## v =
  ## (psi_j1,...,psi_jJ,sigma_j)
  ##
  ## Returns exact analytical CDF and density.
  ## ============================================================

  make_exact_sum_sd_biexp <- function(
      psi_j
  ) {

    psi_j <- as.numeric(psi_j)

    norm2_j <- sum(psi_j^2)

    if (!is.finite(norm2_j) ||
        norm2_j >= 1) {
      stop("||psi_j|| must be smaller than 1.")
    }

    sigma_j <-
      sqrt(
        1 - norm2_j
      )

    v_full <-
      c(
        psi_j,
        sigma_j
      )


    ## ----------------------------------------------------------
    ## Zero coefficients contribute nothing
    ## ----------------------------------------------------------

    active <-
      abs(v_full) >
      zero_tol

    v <-
      v_full[active]

    if (length(v) == 0L) {
      stop("No active coefficients.")
    }


    ## ----------------------------------------------------------
    ## Transformed positive/negative exponential scales
    ##
    ## v > 0:
    ##
    ## alpha = a v
    ## beta  = b v
    ##
    ## v < 0:
    ##
    ## alpha = -b v
    ## beta  = -a v
    ## ----------------------------------------------------------

    positive <- v > 0

    alpha <-
      ifelse(
        positive,
        a0 * v,
        -b0 * v
      )

    beta <-
      ifelse(
        positive,
        b0 * v,
        -a0 * v
      )


    Acoef <-
      get_pf_coef(
        alpha,
        sep_tol = sep_tol
      )

    Bcoef <-
      get_pf_coef(
        beta,
        sep_tol = sep_tol
      )


    Lp <- length(alpha)
    Lm <- length(beta)


    ## ----------------------------------------------------------
    ## S = Y - shift0 * sum(v)
    ##
    ## therefore
    ##
    ## y = x + shift0 * sum(v)
    ## ----------------------------------------------------------

    shift <-
      shift0 *
      sum(v)


    ## ==========================================================
    ## Exact CDF -- vectorized in x
    ## ==========================================================

    cdf <- function(x) {

      xx <- as.numeric(x)

      y <-
        xx +
        shift

      out <-
        numeric(
          length(y)
        )


      left <- y < 0
      right <- !left


      ## --------------------------------------------------------
      ## y < 0
      ## --------------------------------------------------------

      if (any(left)) {

        yy <- y[left]

        val <-
          numeric(
            length(yy)
          )

        for (ii in seq_len(Lp)) {

          ai <- alpha[ii]

          for (jj in seq_len(Lm)) {

            bj <- beta[jj]

            val <-
              val +
              Acoef[ii] *
              Bcoef[jj] *
              bj /
              (ai + bj) *
              exp(
                yy / bj
              )
          }
        }

        out[left] <- val
      }


      ## --------------------------------------------------------
      ## y >= 0
      ## --------------------------------------------------------

      if (any(right)) {

        yy <- y[right]

        tail <-
          numeric(
            length(yy)
          )

        for (ii in seq_len(Lp)) {

          ai <- alpha[ii]

          for (jj in seq_len(Lm)) {

            bj <- beta[jj]

            tail <-
              tail +
              Acoef[ii] *
              Bcoef[jj] *
              ai /
              (ai + bj) *
              exp(
                -yy / ai
              )
          }
        }

        out[right] <-
          1 -
          tail
      }


      ## round-off protection only
      out <-
        pmin(
          pmax(out, 0),
          1
        )

      out
    }


    ## ==========================================================
    ## Exact density -- vectorized in x
    ## ==========================================================

    density <- function(x) {

      xx <- as.numeric(x)

      y <-
        xx +
        shift

      out <-
        numeric(
          length(y)
        )

      left <- y < 0
      right <- !left


      ## --------------------------------------------------------
      ## y < 0
      ## --------------------------------------------------------

      if (any(left)) {

        yy <- y[left]

        val <-
          numeric(
            length(yy)
          )

        for (ii in seq_len(Lp)) {

          ai <- alpha[ii]

          for (jj in seq_len(Lm)) {

            bj <- beta[jj]

            val <-
              val +
              Acoef[ii] *
              Bcoef[jj] /
              (ai + bj) *
              exp(
                yy / bj
              )
          }
        }

        out[left] <- val
      }


      ## --------------------------------------------------------
      ## y >= 0
      ## --------------------------------------------------------

      if (any(right)) {

        yy <- y[right]

        val <-
          numeric(
            length(yy)
          )

        for (ii in seq_len(Lp)) {

          ai <- alpha[ii]

          for (jj in seq_len(Lm)) {

            bj <- beta[jj]

            val <-
              val +
              Acoef[ii] *
              Bcoef[jj] /
              (ai + bj) *
              exp(
                -yy / ai
              )
          }
        }

        out[right] <- val
      }

      out
    }


    list(
      psi = psi_j,
      sigma = sigma_j,
      v = v,
      alpha = alpha,
      beta = beta,
      Acoef = Acoef,
      Bcoef = Bcoef,
      shift = shift,
      cdf = cdf,
      density = density
    )
  }


  ## ============================================================
  ## 4. Simulate common and idiosyncratic BiExp factors
  ## ============================================================

  Z <- matrix(
    r_sd_biexp(n * J),
    nrow = n,
    ncol = J
  )

  Eps <- matrix(
    r_sd_biexp(n * k),
    nrow = n,
    ncol = k
  )


  ## ============================================================
  ## 5. Construct standardized latent variables
  ##
  ## S_j =
  ## psi_j' Z + sigma_j epsilon_j
  ##
  ## Var(S_j)=1
  ## ============================================================

  S <-
    Z %*%
    t(psi) +
    sweep(
      Eps,
      2,
      sigma,
      "*"
    )


    ## ============================================================
  ## 6. Exact probability integral transform
  ## ============================================================

  U_factor <- matrix(
    NA_real_,
    nrow = n,
    ncol = k
  )

  S_dist <- vector("list", k)

  for (j in seq_len(k)) {

    Sobj <- make_exact_sum_sd_biexp(
      psi_j = psi[j, ]
    )

    S_dist[[j]] <- Sobj

    U_factor[, j] <- Sobj$cdf(
      S[, j]
    )
  }


  ## ------------------------------------------------------------
  ## Numerical endpoint protection only
  ## Keep matrix dimensions
  ## ------------------------------------------------------------

  U_factor <- pmax(
    U_factor,
    1e-12
  )

  U_factor <- pmin(
    U_factor,
    1 - 1e-12
  )

  U_factor <- matrix(
    U_factor,
    nrow = n,
    ncol = k
  )  ## ============================================================
  ## 7. Checkerboard cell indices
  ## ============================================================

  M <-
    ceiling(
      d * U_factor
    )

  M[M < 1] <- 1
  M[M > d] <- d


  ## ============================================================
  ## 8. Uniform randomization inside each selected cell
  ## ============================================================

  V <-
    matrix(
      runif(n * k),
      nrow = n,
      ncol = k
    )

  U_cb <-
    (M - 1 + V) /
    d


  ## ============================================================
  ## 9. Marginal latent Lambda quantiles
  ## ============================================================

  qweibull_mu_phi_simu <- function(
      u,
      mu,
      phi
  ) {

    c0 <-
      mu /
      gamma(1 + phi)

    c0 *
      (-log1p(-u))^phi
  }


  qlatent_one_margin <- function(
      u,
      mu_j,
      phi_j,
      dist_j
  ) {

    if (dist_j == "gamma") {

      shape <-
        1 / phi_j

      rate <-
        1 /
        (phi_j * mu_j)

      out <-
        qgamma(
          u,
          shape = shape,
          rate = rate
        )


    } else if (dist_j == "lnorm") {

      sdlog <-
        sqrt(
          log(
            1 + phi_j
          )
        )

      meanlog <-
        log(mu_j) -
        0.5 * sdlog^2

      out <-
        qlnorm(
          u,
          meanlog = meanlog,
          sdlog = sdlog
        )


    } else if (dist_j == "weibull") {

      out <-
        qweibull_mu_phi_simu(
          u,
          mu = mu_j,
          phi = phi_j
        )


    } else if (dist_j == "ig") {

      out <-
        statmod::qinvgauss(
          p = u,
          mean = mu_j,
          shape = 1 / phi_j
        )


    } else if (dist_j == "pareto") {

      if (!exists(
        "qpareto_mu_phi_simu"
      )) {

        stop(
          paste0(
            "Please define ",
            "qpareto_mu_phi_simu(u, mu, phi) ",
            "before using dist = 'pareto'."
          )
        )
      }

      out <-
        qpareto_mu_phi_simu(
          u,
          mu = mu_j,
          phi = phi_j
        )
    }


    out
  }


  ## ============================================================
  ## 10. Simulate latent Lambda
  ## ============================================================

  Lambda <-
    matrix(
      NA_real_,
      nrow = n,
      ncol = k
    )


  for (j in seq_len(k)) {

    Lambda[, j] <-
      qlatent_one_margin(
        u = U_cb[, j],
        mu_j = mu[, j],
        phi_j = phi[j],
        dist_j = dist[j]
      )
  }


  ## ============================================================
  ## 11. Mixed Poisson counts
  ## ============================================================

  Xcount <-
    matrix(
      rpois(
        n * k,
        lambda = as.vector(Lambda)
      ),
      nrow = n,
      ncol = k
    )


  ## ============================================================
  ## 12. Names
  ## ============================================================

  colnames(Xcount) <-
    paste0(
      "X",
      seq_len(k)
    )

  colnames(Lambda) <-
    paste0(
      "Lambda",
      seq_len(k)
    )

  colnames(U_cb) <-
    paste0(
      "U",
      seq_len(k)
    )

  colnames(U_factor) <-
    paste0(
      "U_biexp",
      seq_len(k)
    )

  colnames(S) <-
    paste0(
      "S",
      seq_len(k)
    )

  colnames(M) <-
    paste0(
      "M",
      seq_len(k)
    )


  ## ============================================================
  ## 13. Output
  ## ============================================================

  out <-
    as.data.frame(
      Xcount
    )


  if (return_latent) {

    out <-
      cbind(
        out,
        as.data.frame(Lambda),
        as.data.frame(U_cb),
        as.data.frame(U_factor),
        as.data.frame(S),
        as.data.frame(M)
      )
  }


  attr(out, "dist") <- dist

  attr(out, "phi") <- phi

  attr(out, "psi") <- psi

  attr(out, "sigma") <- sigma

  attr(out, "lambda") <- lambda

  attr(out, "biexp_skewness") <-
    2 * (
      lambda^3 -
      (1 - lambda^2)^(3 / 2)
    )

  attr(out, "checkerboard_d") <- d

  attr(out, "J") <- J

  attr(out, "S_distribution") <- S_dist


  out
}



# các hàm phụ
map_theta_copula_in_ff <- function(par, family, eps = 1e-6) {
  fam <- as.integer(family)
  
  softplus <- function(x) log1p(exp(-abs(x))) + pmax(x, 0)
  logistic <- function(x) 1 / (1 + exp(-x))
  
  # 1 tham số
  if (fam == 1L) {                      # Gaussian
    return((1 - eps) * tanh(par[1]))
  }
  if (fam %in% c(3, 13)) {              # Clayton, survival Clayton
    return(softplus(par[1]) + eps)
  }
  if (fam %in% c(23, 33)) {             # rotated Clayton
    return((softplus(par[1]) + eps))
  }
  if (fam %in% c(4, 14)) {              # Gumbel, survival Gumbel
    return(1 + softplus(par[1]) + eps)
  }
  if (fam %in% c(24, 34)) {             # rotated Gumbel
    return((1 + softplus(par[1]) + eps))
  }
  if (fam == 5L) {                      # Frank
    x <- par[1]
    return(ifelse(x >= 0, x + eps, x - eps))
  }
  if (fam %in% c(6, 16)) {              # Joe, survival Joe
    return(1 + softplus(par[1]) + eps)
  }
  if (fam %in% c(26, 36)) {             # rotated Joe
    return(-(1 + softplus(par[1]) + eps))
  }
  
  # 2 tham số
  if (fam == 2L) {                      # Student-t
    rho <- (1 - eps) * tanh(par[1])
    return(c(rho))
  }
  if (fam %in% c(7, 17)) {              # BB1
    return(c(softplus(par[1]) + eps,
             1 + softplus(par[2]) + eps))
  }
  if (fam %in% c(27, 37)) {             # rotated BB1
    return(c(-(softplus(par[1]) + eps),
             -(1 + softplus(par[2]) + eps)))
  }
  if (fam %in% c(8, 18)) {              # BB6
    return(c(1 + softplus(par[1]) + eps,
             1 + softplus(par[2]) + eps))
  }
  if (fam %in% c(28, 38)) {             # rotated BB6
    return(c(-(1 + softplus(par[1]) + eps),
             -(1 + softplus(par[2]) + eps)))
  }
  if (fam %in% c(9, 19)) {              # BB7
    return(c(1 + softplus(par[1]) + eps,
             softplus(par[2]) + eps))
  }
  if (fam %in% c(29, 39)) {             # rotated BB7
    return(c(-(1 + softplus(par[1]) + eps),
             -(softplus(par[2]) + eps)))
  }
  if (fam %in% c(10, 20)) {             # BB8
    return(c(1 + softplus(par[1]) + eps,
             eps + (1 - 2 * eps) * logistic(par[2])))
  }
  if (fam %in% c(30, 40)) {             # rotated BB8
    return(c(-(1 + softplus(par[1]) + eps),
             -(eps + (1 - 2 * eps) * logistic(par[2]))))
  }
  if (fam %in% c(104, 114, 204, 214)) { # Tawn
    return(c(1 + softplus(par[1]) + eps,
             eps + (1 - 2 * eps) * logistic(par[2])))
  }
  if (fam %in% c(124, 134, 224, 234)) { # rotated Tawn
    return(c(-(1 + softplus(par[1]) + eps),
             eps + (1 - 2 * eps) * logistic(par[2])))
  }
  
  stop("Unsupported copula family.")
}

theta_is_valid_for_ff <- function(theta, family) {
  out <- try(
    VineCopula::BiCopCheck(
      family = as.integer(family),
      par    = theta[1],
      par2   = if (length(theta) >= 2) theta[2] else 0
    ),
    silent = TRUE
  )
  !inherits(out, "try-error")
}
