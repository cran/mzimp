
p11_from_rho <- function(rho, p1, p2){
  q1 <- qnorm(p1)
  q2 <- qnorm(p2)
  
  pbivnorm::pbivnorm(
    x = q1,
    y = q2,
    rho = rho
  )
}

# ============================================================
# Extract marginal model from fit_zimp()
# ============================================================

.extract_zimp_margin <- function(fit) {
  model <- if (!is.null(fit$best_model)) fit$best_model else fit

  if (is.null(fit$y) || is.null(fit$X) || is.null(fit$zero_inflation))
    stop("Invalid univariate marginal fit.")

  dist <- model$distribution
  if (is.null(dist)) dist <- model$mixing.distribution
  if (is.null(dist)) dist <- fit$mixing.distribution

  X <- as.matrix(fit$X)
  storage.mode(X) <- "double"

  Z <- fit$Z
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    storage.mode(Z) <- "double"
  }

  # Current fit_zimp stores the optimizer parameter vector in $par.
  # Its order is:
  #   beta (ncol(X)), log(phi), alpha (ncol(Z), if ZIMP).
  # Do not require duplicated $be/$al fields to exist.
  par <- model$par
  if (is.null(par)) par <- fit$par

  p <- ncol(X)
  q <- if (isTRUE(fit$zero_inflation)) ncol(Z) else 0L

  be <- model$be
  if (is.null(be)) be <- model$beta
  if (is.null(be) && !is.null(par) && length(par) >= p)
    be <- par[seq_len(p)]

  phi <- model$phi
  if (is.null(phi) && !is.null(par) && length(par) >= p + 1L)
    phi <- exp(par[p + 1L])

  al <- model$al
  if (is.null(al)) al <- model$alpha
  if (isTRUE(fit$zero_inflation) && is.null(al) &&
      !is.null(par) && length(par) >= p + 1L + q) {
    al <- par[p + 1L + seq_len(q)]
  }

  if (is.null(be) || !is.numeric(be) || any(!is.finite(be)))
    stop("Could not extract numeric count coefficients from fit_zimp().")

  if (is.null(phi) || !is.numeric(phi) || length(phi) != 1L ||
      !is.finite(phi) || phi <= 0)
    stop("Could not extract a valid dispersion parameter phi from fit_zimp().")

  if (isTRUE(fit$zero_inflation) &&
      (is.null(al) || !is.numeric(al) || length(al) != q ||
       any(!is.finite(al)))) {
    stop(
      "Could not extract zero-inflation coefficients from fit_zimp(). ",
      "Expected ", q, " coefficient(s)."
    )
  }

  list(
    fit = fit,
    y = as.integer(fit$y),
    X = X,
    Z = Z,
    zero_inflation = isTRUE(fit$zero_inflation),
    dist = dist,
    distribution = dist,
    var_power = if (!is.null(fit$var.power)) fit$var.power else fit$var_power,
    use_custom_var_power = !is.null(if (!is.null(fit$var.power)) fit$var.power else fit$var_power),
    be = as.numeric(be),
    beta = as.numeric(be),
    al = if (is.null(al)) NULL else as.numeric(al),
    alpha = if (is.null(al)) NULL else as.numeric(al),
    phi = as.numeric(phi),
    neg_log_lik_value = model$neg_log_lik_value,
    log_lik_value = model$log_lik_value
  )
}


# ============================================================
# Compute marginal I matrix
# ============================================================

.compute_I_margin <- function(y, X, be, phi, dist, d,
                              use_custom_var_power = FALSE,
                              var_power = NULL,
                              n_gh = 20,
                              gh = NULL) {
  
  arg_vp <- if (use_custom_var_power) list(var_power = var_power) else list()
  
  if (dist == "gamma") {
    out <- do.call(I_gamma, c(
      list(z = y, X = X, be = be, phi = phi, d = d),
      arg_vp
    ))
  } else if (dist == "lnorm") {
    out <- do.call(I_lnorm, c(
      list(z = y, X = X, be = be, phi = phi, d = d, n_gh = n_gh, gh = gh),
      arg_vp
    ))
  } else if (dist == "ig") {
    out <- do.call(I_ig, c(
      list(z = y, X = X, be = be, phi = phi, d = d),
      arg_vp
    ))
  } else if (dist == "weibull") {
    out <- do.call(I_weibull, c(
      list(z = y, X = X, be = be, phi = phi, d = d),
      arg_vp
    ))
  } else if (dist == "loglogistic") {
    out <- do.call(I_loglogistic, c(
      list(z = y, X = X, be = be, phi = phi, d = d),
      arg_vp
    ))
  } else if (dist == "betaprime") {
    if (!exists("I_betaprime", mode = "function")) {
      stop("Distribution 'betaprime' was selected, but I_betaprime() is not available.")
    }
    out <- do.call(I_betaprime, c(
      list(z = y, X = X, be = be, phi = phi, d = d),
      arg_vp
    ))
  } else {
    stop("Unsupported marginal distribution: ", dist)
  }
  
  out
}



# ============================================================
# Extract two selected margins from select_biv_margin()
# ============================================================

.extract_selected_biv_margin <- function(margin) {

  if (!inherits(margin, "biv_margin_fit")) {
    stop("`margin` must be an object returned by select_biv_margin().")
  }

  if (is.null(margin$k) || length(margin$k) != 1L ||
      !is.finite(margin$k) || as.integer(margin$k) != 2L) {
    stop("`margin` must contain exactly two margins (k = 2).")
  }

  if (is.null(margin$selected_margins) ||
      !is.list(margin$selected_margins) ||
      length(margin$selected_margins) != 2L ||
      any(vapply(margin$selected_margins, is.null, logical(1)))) {
    stop("`margin$selected_margins` must contain two selected margins.")
  }

  validate_one <- function(x, j) {

    required <- c(
      "y", "X", "zero_inflation",
      "dist", "var_power", "use_custom_var_power",
      "be", "phi", "neg_log_lik_value", "log_lik_value"
    )

    missing_fields <- required[
      !vapply(required, function(nm) !is.null(x[[nm]]), logical(1))
    ]

    if (length(missing_fields) > 0L) {
      stop(
        "Selected margin ", j, " is missing: ",
        paste(missing_fields, collapse = ", "), "."
      )
    }

    x$y <- as.integer(x$y)
    x$X <- as.matrix(x$X)
    if (!is.null(x$Z)) x$Z <- as.matrix(x$Z)

    if (length(x$y) != nrow(x$X)) {
      stop("Selected margin ", j, " has inconsistent y and X.")
    }

    if (!is.null(x$Z) && nrow(x$Z) != length(x$y)) {
      stop("Selected margin ", j, " has inconsistent y and Z.")
    }

    if (!is.numeric(x$be) || any(!is.finite(x$be))) {
      stop("Invalid beta estimate in selected margin ", j, ".")
    }

    if (!is.numeric(x$phi) || length(x$phi) != 1L ||
        !is.finite(x$phi) || x$phi <= 0) {
      stop("Invalid phi estimate in selected margin ", j, ".")
    }

    if (isTRUE(x$zero_inflation)) {
      if (is.null(x$Z)) {
        stop("Selected margin ", j, " is zero-inflated but has no Z.")
      }
      if (is.null(x$al) || !is.numeric(x$al) || any(!is.finite(x$al))) {
        stop("Selected margin ", j, " has invalid alpha estimates.")
      }
    }

    x
  }

  margin1 <- validate_one(margin$selected_margins[[1L]], 1L)
  margin2 <- validate_one(margin$selected_margins[[2L]], 2L)

  if (length(margin1$y) != length(margin2$y)) {
    stop("The two selected marginal fits must use the same number of observations.")
  }

  list(margin1 = margin1, margin2 = margin2)
}



# ============================================================
# VineCopula helpers for latent-intensity copula selection
# ============================================================

.vine_all_biv_families <- c(
  0L,1L,2L,3L,4L,5L,6L,7L,8L,9L,10L,
  13L,14L,16L,17L,18L,19L,20L,
  23L,24L,26L,27L,28L,29L,30L,
  33L,34L,36L,37L,38L,39L,40L,
  104L,114L,124L,134L,204L,214L,224L,234L
)

.vine_positive_families <- c(
  3L,4L,6L,7L,8L,9L,10L,
  13L,14L,16L,17L,18L,19L,20L,
  104L,114L,204L,214L
)

.vine_negative_families <- c(
  23L,24L,26L,27L,28L,29L,30L,
  33L,34L,36L,37L,38L,39L,40L,
  124L,134L,224L,234L
)

.vine_flexible_families <- c(0L,1L,2L,5L)

.vine_filter_by_gaussian_sign <- function(families, rho, sign.tol=1e-6) {
  families <- unique(as.integer(families))
  bad <- setdiff(families, .vine_all_biv_families)
  if (length(bad)) stop('Unsupported VineCopula family code(s): ', paste(bad, collapse=', '))
  families <- unique(c(1L, families))
  if (!is.finite(rho) || abs(rho) <= sign.tol) return(families)
  keep <- if (rho > 0) c(.vine_flexible_families, .vine_positive_families) else c(.vine_flexible_families, .vine_negative_families)
  families[families %in% keep]
}

.vine_family_spec <- function(family, gaussian.rho=0, eps=1e-5, max.par=50) {
  family <- as.integer(family)
  rho.pos <- gaussian.rho >= 0
  if (family==0L) return(list(npar=0L,start=numeric(0),lower=numeric(0),upper=numeric(0)))
  if (family==1L) return(list(npar=1L,start=pmin(pmax(gaussian.rho,-0.8),0.8),lower=-0.999,upper=0.999))
  if (family==2L) return(list(npar=2L,start=c(pmin(pmax(gaussian.rho,-0.8),0.8),4),lower=c(-0.999,3),upper=c(0.999,30)))
  if (family==5L) {
    if (rho.pos) return(list(npar=1L,start=1,lower=eps,upper=max.par))
    return(list(npar=1L,start=-1,lower=-max.par,upper=-eps))
  }
  if (family %in% c(3L,13L)) return(list(npar=1L,start=1,lower=eps,upper=max.par))
  if (family %in% c(23L,33L)) return(list(npar=1L,start=-1,lower=-max.par,upper=-eps))
  if (family %in% c(4L,14L)) return(list(npar=1L,start=1.2,lower=1,upper=max.par))
  if (family %in% c(24L,34L)) return(list(npar=1L,start=-1.2,lower=-max.par,upper=-1))
  if (family %in% c(6L,16L)) return(list(npar=1L,start=1.2,lower=1+eps,upper=max.par))
  if (family %in% c(26L,36L)) return(list(npar=1L,start=-1.2,lower=-max.par,upper=-1-eps))
  if (family %in% c(7L,17L)) return(list(npar=2L,start=c(1,1.5),lower=c(eps,1),upper=c(max.par,max.par)))
  if (family %in% c(27L,37L)) return(list(npar=2L,start=c(-1,-1.5),lower=c(-max.par,-max.par),upper=c(-eps,-1)))
  if (family %in% c(8L,18L)) return(list(npar=2L,start=c(1.2,1.2),lower=c(1,1),upper=c(max.par,max.par)))
  if (family %in% c(28L,38L)) return(list(npar=2L,start=c(-1.2,-1.2),lower=c(-max.par,-max.par),upper=c(-1,-1)))
  if (family %in% c(9L,19L)) return(list(npar=2L,start=c(1.2,1),lower=c(1,eps),upper=c(max.par,max.par)))
  if (family %in% c(29L,39L)) return(list(npar=2L,start=c(-1.2,-1),lower=c(-max.par,-max.par),upper=c(-1,-eps)))
  if (family %in% c(10L,20L)) return(list(npar=2L,start=c(1.2,0.5),lower=c(1,eps),upper=c(max.par,1)))
  if (family %in% c(30L,40L)) return(list(npar=2L,start=c(-1.2,-0.5),lower=c(-max.par,-1),upper=c(-1,-eps)))
  if (family %in% c(104L,114L,204L,214L)) return(list(npar=2L,start=c(1.2,0.5),lower=c(1,0),upper=c(max.par,1)))
  if (family %in% c(124L,134L,224L,234L)) return(list(npar=2L,start=c(-1.2,0.5),lower=c(-max.par,0),upper=c(-1,1)))
  stop('No parameter specification implemented for family ', family)
}

.vine_checkerboard_density <- function(family, par=0, par2=0, d=10L) {

  if (!requireNamespace("VineCopula", quietly=TRUE)) {
    stop("Package 'VineCopula' is required.")
  }

  family <- as.integer(family)
  d <- as.integer(d)

  if (family == 0L) {
    return(matrix(1, nrow=d, ncol=d))
  }

  g <- seq(0, 1, length.out=d + 1L)

  # Build C(u,v) without evaluating VineCopula exactly at 0 or 1.
  # Boundary values follow exactly from copula identities.
  Cmat <- matrix(0, nrow=d + 1L, ncol=d + 1L)

  # Rows correspond to u, columns to v.
  Cmat[d + 1L, ] <- g       # C(1,v) = v
  Cmat[, d + 1L] <- g       # C(u,1) = u

  if (d >= 2L) {

    gi <- g[2:d]

    uu <- rep(gi, times=length(gi))
    vv <- rep(gi, each=length(gi))

    CC <- VineCopula::BiCopCDF(
      u1=uu,
      u2=vv,
      family=family,
      par=par,
      par2=par2,
      check.pars=TRUE
    )

    if (any(!is.finite(CC))) {
      stop(
        "BiCopCDF returned non-finite values for family ",
        family, ", par=", par, ", par2=", par2
      )
    }

    Cmat[2:d, 2:d] <- matrix(
      CC,
      nrow=d - 1L,
      ncol=d - 1L
    )
  }

  mass <- Cmat[-1L,-1L,drop=FALSE] -
    Cmat[-(d+1L),-1L,drop=FALSE] -
    Cmat[-1L,-(d+1L),drop=FALSE] +
    Cmat[-(d+1L),-(d+1L),drop=FALSE]

  mass[abs(mass) < 1e-13] <- 0

  if (any(!is.finite(mass))) {
    stop("Non-finite checkerboard mass for family ", family)
  }

  if (any(mass < -1e-8)) {
    stop(
      "Negative checkerboard mass for family ", family,
      "; min=", signif(min(mass),6),
      "; par=", signif(par,6),
      "; par2=", signif(par2,6)
    )
  }

  mass[mass < 0] <- 0

  sm <- sum(mass)
  if (!is.finite(sm) || sm <= 0) {
    stop("Invalid total checkerboard mass for family ", family)
  }

  # Numerical renormalization.
  mass <- mass / sm

  d^2 * mass
}


.vine_family_name <- function(family) {
  if (!requireNamespace('VineCopula', quietly=TRUE)) return(as.character(family))
  z <- try(VineCopula::BiCopName(as.integer(family)), silent=TRUE)
  if (inherits(z,'try-error')) as.character(family) else as.character(z)
}

# ============================================================
# Main function
# ============================================================

.fit_bizimp_dependence_engine <- function(
    margin = NULL,
    fit_margin_1 = NULL,
    fit_margin_2 = NULL,
    n_gh = 20,
    gh = NULL,
    copula.family.lambda.list = c(1, 2, 3, 4, 5),
    exact.copula.family = NULL,
    B.depen.option = c("normal", "comonotone", "countermonotone"),
    counting.depen.option = c("bicopula"),
    d0 = 10,
    method = "IFM",
    gaussian.sign.tol = 1e-6,
    t.df.grid = c(3,5,10,30),
    progress = TRUE,
    trace = 1,
    eval.max = 200,
    iter.max = 100
) {
  
  
  B.depen.option <- match.arg(B.depen.option)
  counting.depen.option <- match.arg(counting.depen.option)

  .progress_msg <- function(...) {
    if (isTRUE(progress)) {
      cat(..., "\n", sep = "")
      flush.console()
    }
  }
  
  # ------------------------------------------------------------
  # 1. Obtain already-fitted margins
  #
  # Marginal estimation is performed outside this dependence engine.
  # Public fit_bizimp() fits both margins with fit_zimp() and passes the
  # resulting objects through fit_margin_1 and fit_margin_2.  Alternatively,
  # a selected bivariate margin object may be supplied through margin.
  # ------------------------------------------------------------

  if (!is.null(margin)) {

    if (!inherits(margin, "biv_margin_fit")) {
      stop("`margin` must be NULL or an object returned by select_biv_margin().")
    }

    if (!is.null(fit_margin_1) || !is.null(fit_margin_2)) {
      stop(
        "When `margin` is supplied, do not also supply ",
        "`fit_margin_1` or `fit_margin_2`."
      )
    }

    mm <- .extract_selected_biv_margin(margin)
    margin1 <- mm$margin1
    margin2 <- mm$margin2

  } else {

    if (is.null(fit_margin_1) || is.null(fit_margin_2)) {
      stop(
        "When `margin` is NULL, both `fit_margin_1` and `fit_margin_2` ",
        "must be supplied as fitted objects returned by fit_zimp()."
      )
    }

    margin1 <- .extract_zimp_margin(fit_margin_1)
    margin2 <- .extract_zimp_margin(fit_margin_2)
  }

  y1 <- margin1$y
  y2 <- margin2$y

  if (length(y1) != length(y2)) {
    stop("The two marginal fits must use the same number of observations.")
  }

  n <- length(y1)
  Y <- cbind(y1, y2)

  zi1 <- isTRUE(margin1$zero_inflation)
  zi2 <- isTRUE(margin2$zero_inflation)

  marginal_pattern <- if (zi1 && zi2) {
    "zimp_zimp"
  } else if (zi1) {
    "zimp_pm"
  } else if (zi2) {
    "pm_zimp"
  } else {
    "pm_pm"
  }

  estimate_B_dependence <- zi1 && zi2

  X1 <- margin1$X
  X2 <- margin2$X
  Z1 <- if (zi1) margin1$Z else NULL
  Z2 <- if (zi2) margin2$Z else NULL

  be1 <- margin1$be
  be2 <- margin2$be
  al1 <- if (zi1) margin1$al else NULL
  al2 <- if (zi2) margin2$al else NULL
  phi1 <- margin1$phi
  phi2 <- margin2$phi

  # ------------------------------------------------------------
  # 2. Utilities
  # ------------------------------------------------------------
  
  expit <- function(eta) {
    eta <- pmax(pmin(eta, 35), -35)
    1 / (1 + exp(-eta))
  }
  
  p11_from_B_option <- function(pb1, pb2, option, rho = NULL) {
    rho_use <- if (option == "normal") {
      if (is.null(rho) || length(rho) != 1L || !is.finite(rho))
        stop("rho must be supplied for Gaussian zero dependence.")
      as.numeric(rho)
    } else if (option == "comonotone") {
      0.999
    } else if (option == "countermonotone") {
      -0.999
    } else {
      stop("Unknown B.depen.option.")
    }

    p11_from_rho(rho = rho_use, p1 = pb1, p2 = pb2)
  }
  
  # ------------------------------------------------------------
  # 3. Compute marginal mixed-Poisson probabilities
  # ------------------------------------------------------------
  
  Pn1 <- .compute_I_margin(
    y = y1, X = X1, be = be1, phi = phi1,
    dist = margin1$dist, d = 1,
    use_custom_var_power = margin1$use_custom_var_power,
    var_power = margin1$var_power,
    n_gh = n_gh, gh = gh
  )
  Pn1 <- as.numeric(Pn1)
  
  I1.vec <- .compute_I_margin(
    y = y1, X = X1, be = be1, phi = phi1,
    dist = margin1$dist, d = d0,
    use_custom_var_power = margin1$use_custom_var_power,
    var_power = margin1$var_power,
    n_gh = n_gh, gh = gh
  )
  I1.vec <- I1.vec / rowSums(I1.vec) * Pn1
  
  Pn2 <- .compute_I_margin(
    y = y2, X = X2, be = be2, phi = phi2,
    dist = margin2$dist, d = 1,
    use_custom_var_power = margin2$use_custom_var_power,
    var_power = margin2$var_power,
    n_gh = n_gh, gh = gh
  )
  Pn2 <- as.numeric(Pn2)
  
  I2.vec <- .compute_I_margin(
    y = y2, X = X2, be = be2, phi = phi2,
    dist = margin2$dist, d = d0,
    use_custom_var_power = margin2$use_custom_var_power,
    var_power = margin2$var_power,
    n_gh = n_gh, gh = gh
  )
  I2.vec <- I2.vec / rowSums(I2.vec) * Pn2
  
  # ------------------------------------------------------------
  # 4. Zero-inflation probabilities
  #
  # pb_j = P(B_j = 1), where B_j = 1 means that margin j
  # enters the mixed-Poisson count process.
  #
  # For a non-zero-inflated PM margin, B_j is identically 1.
  # ------------------------------------------------------------

  pb1 <- if (zi1) {
    as.numeric(expit(Z1 %*% al1))
  } else {
    rep(1, n)
  }

  pb2 <- if (zi2) {
    as.numeric(expit(Z2 %*% al2))
  } else {
    rep(1, n)
  }
  
  # ------------------------------------------------------------
  # 5. Fit dependence models: Gaussian prescreen + sign filtering
  # ------------------------------------------------------------

  exact_mode <- !is.null(exact.copula.family)
  requested_families <- if (exact_mode) as.integer(exact.copula.family) else
    unique(as.integer(copula.family.lambda.list))
  bad <- setdiff(requested_families, .vine_all_biv_families)
  if (length(bad)) stop('Unsupported VineCopula family code(s): ', paste(bad, collapse=', '))
  if (!exact_mode) requested_families <- unique(c(1L, requested_families))

  t.df.grid <- sort(unique(as.integer(t.df.grid)))
  t.df.grid <- t.df.grid[is.finite(t.df.grid) & t.df.grid > 2L]
  if (2L %in% requested_families && !length(t.df.grid)) stop('`t.df.grid` must contain integers > 2.')

  .biv_observed_nll <- function(density_mat, theta.B=NA_real_) {
    big_penalty <- 1e12; tiny <- 1e-300
    if (any(!is.finite(density_mat)) || any(density_mat < -1e-10)) return(big_penalty)
    density_mat[density_mat < 0] <- 0
    Pn_y1_y2 <- rowSums((I1.vec %*% density_mat) * I2.vec)
    if (any(!is.finite(Pn_y1_y2)) || any(Pn_y1_y2 < 0)) return(big_penalty)

    if (!zi1 && !zi2) {
      prob <- Pn_y1_y2
    } else if (zi1 && !zi2) {
      prob <- ifelse(y1>0, pb1*Pn_y1_y2, (1-pb1)*Pn2 + pb1*Pn_y1_y2)
    } else if (!zi1 && zi2) {
      prob <- ifelse(y2>0, pb2*Pn_y1_y2, (1-pb2)*Pn1 + pb2*Pn_y1_y2)
    } else {
      Cb <- try(p11_from_B_option(pb1,pb2,B.depen.option,
        rho=if (B.depen.option=='normal') theta.B else NULL),silent=TRUE)
      if (inherits(Cb,'try-error')) return(big_penalty)
      lo <- pmax(pb1+pb2-1,0); up <- pmin(pb1,pb2)
      if (any(!is.finite(Cb)) || any(Cb<lo-1e-10) || any(Cb>up+1e-10)) return(big_penalty)
      Cb <- pmin(pmax(Cb,lo),up)
      is11 <- y1>0 & y2>0; is10 <- y1>0 & y2==0; is01 <- y1==0 & y2>0; is00 <- y1==0 & y2==0
      prob <- numeric(n)
      prob[is11] <- Cb[is11]*Pn_y1_y2[is11]
      prob[is10] <- pb1[is10]*Pn1[is10] + Cb[is10]*(Pn_y1_y2[is10]-Pn1[is10])
      prob[is01] <- pb2[is01]*Pn2[is01] + Cb[is01]*(Pn_y1_y2[is01]-Pn2[is01])
      prob[is00] <- 1-pb1[is00]-pb2[is00]+pb1[is00]*Pn1[is00]+pb2[is00]*Pn2[is00]+Cb[is00]*(1-Pn1[is00]-Pn2[is00]+Pn_y1_y2[is00])
    }
    if (any(!is.finite(prob)) || any(prob<=tiny)) return(big_penalty)
    v <- -sum(log(prob)); if (is.finite(v)) v else big_penalty
  }

  .fit_one_vine_family <- function(family, gaussian.rho=0) {
    family <- as.integer(family)
    hasB <- estimate_B_dependence && B.depen.option=='normal'
    fixedB <- if (estimate_B_dependence && !hasB) {
      if (B.depen.option == 'comonotone') 0.999 else -0.999
    } else NA_real_

    if (family==0L) {
      W <- matrix(1,d0,d0)
      if (hasB) {
        re <- nlminb(0,function(rhoB).biv_observed_nll(W,rhoB),lower=-0.999,upper=0.999,
                     control=list(trace=trace,eval.max=eval.max,iter.max=iter.max))
        return(list(family=0L,family_name=.vine_family_name(0L),par=0,par2=0,theta.B=re$par,
                    neg_log_lik_value=re$objective,log_lik_value=-re$objective,convergence=re$convergence,message=re$message,par_raw=re$par))
      }
      nll <- .biv_observed_nll(W,fixedB)
      return(list(family=0L,family_name=.vine_family_name(0L),par=0,par2=0,theta.B=fixedB,
                  neg_log_lik_value=nll,log_lik_value=-nll,convergence=0L,message='No latent parameter',par_raw=numeric(0)))
    }

    if (family==2L) {
      prof <- lapply(t.df.grid,function(df0){
        if (gaussian.rho>gaussian.sign.tol) {lr<-1e-5; ur<-0.999; sr<-max(min(gaussian.rho,0.8),0.05)}
        else if (gaussian.rho < -gaussian.sign.tol) {lr<--0.999; ur<--1e-5; sr<-min(max(gaussian.rho,-0.8),-0.05)}
        else {lr<--0.999; ur<-0.999; sr<-0}
        st <- if(hasB)c(0,sr) else sr; lo<-if(hasB)c(-0.999,lr) else lr; up<-if(hasB)c(0.999,ur) else ur
        ff<-function(z){rhoB<-if(hasB)z[1] else fixedB; rhoG<-if(hasB)z[2] else z[1]
          W<-try(.vine_checkerboard_density(2L,rhoG,df0,d0),silent=TRUE); if(inherits(W,'try-error'))return(1e12)
          .biv_observed_nll(W,rhoB)}
        re<-nlminb(st,ff,lower=lo,upper=up,control=list(trace=trace,eval.max=eval.max,iter.max=iter.max))
        list(df=df0,re=re)
      })
      ib<-which.min(vapply(prof,function(z)z$re$objective,numeric(1))); re<-prof[[ib]]$re
      rhoB<-if(hasB)re$par[1] else fixedB; rhoG<-if(hasB)re$par[2] else re$par[1]
      return(list(family=2L,family_name=.vine_family_name(2L),par=rhoG,par2=prof[[ib]]$df,theta.B=rhoB,
                  neg_log_lik_value=re$objective,log_lik_value=-re$objective,convergence=re$convergence,message=re$message,par_raw=re$par,t_profile=prof))
    }

    sp <- .vine_family_spec(family,gaussian.rho)
    st <- if(hasB)c(0,sp$start) else sp$start; lo<-if(hasB)c(-0.999,sp$lower) else sp$lower; up<-if(hasB)c(0.999,sp$upper) else sp$upper
    ff<-function(z){rhoB<-if(hasB)z[1] else fixedB; cp<-if(hasB)z[-1] else z; p1<-cp[1]; p2<-if(length(cp)>1)cp[2] else 0
      W<-try(.vine_checkerboard_density(family,p1,p2,d0),silent=TRUE); if(inherits(W,'try-error'))return(1e12)
      .biv_observed_nll(W,rhoB)}
    re<-nlminb(st,ff,lower=lo,upper=up,control=list(trace=trace,eval.max=eval.max,iter.max=iter.max))
    rhoB<-if(hasB)re$par[1] else fixedB; cp<-if(hasB)re$par[-1] else re$par
    list(family=family,family_name=.vine_family_name(family),par=cp[1],par2=if(length(cp)>1)cp[2] else 0,theta.B=rhoB,
         neg_log_lik_value=re$objective,log_lik_value=-re$objective,convergence=re$convergence,message=re$message,par_raw=re$par)
  }

  if (exact_mode) {
    # Exactly one requested latent-intensity copula is fitted.
    fam <- requested_families[1L]

    if (fam == 1L) {
      gaussian_fit <- .fit_one_vine_family(1L, 0)
      gaussian_rho <- gaussian_fit$par
      models <- list(gaussian_fit)
    } else {
      # Gaussian is only a pilot for starting value/sign; it is not a candidate.
      gaussian_fit <- .fit_one_vine_family(1L, 0)
      gaussian_rho <- gaussian_fit$par
      models <- list(.fit_one_vine_family(fam, gaussian_rho))
    }

    fitted_families <- fam
  } else {
    # Legacy selection mode.
    gaussian_fit <- .fit_one_vine_family(1L, 0)
    gaussian_rho <- gaussian_fit$par
    fitted_families <- .vine_filter_by_gaussian_sign(
      requested_families, gaussian_rho, gaussian.sign.tol
    )

    models <- list(gaussian_fit)
    other_families <- setdiff(fitted_families, 1L)
    if (length(other_families)) {
      for (fam in other_families)
        models[[length(models) + 1L]] <- .fit_one_vine_family(fam, gaussian_rho)
    }
  }

  two_par <- c(2L,7L,8L,9L,10L,17L,18L,19L,20L,27L,28L,29L,30L,37L,38L,39L,40L,104L,114L,124L,134L,204L,214L,224L,234L)
  for (j in seq_along(models)) {
    models[[j]]$model_id<-j; models[[j]]$marginal_pattern<-marginal_pattern; models[[j]]$zero_dependence_estimated<-estimate_B_dependence
    models[[j]]$B.depen.option<-if(estimate_B_dependence)B.depen.option else NA_character_; models[[j]]$counting.depen.option<-counting.depen.option
    models[[j]]$copula.family.lambda<-models[[j]]$family; models[[j]]$theta.G<-if(models[[j]]$family %in% two_par)c(models[[j]]$par,models[[j]]$par2) else models[[j]]$par
  }

  summary <- do.call(rbind,lapply(models,function(x)data.frame(
    model_id=x$model_id,marginal_pattern=x$marginal_pattern,zero_dependence_estimated=x$zero_dependence_estimated,
    B.depen.option=x$B.depen.option,counting.depen.option=x$counting.depen.option,theta.B=x$theta.B,
    copula.family.lambda=x$family,copula.family.name=x$family_name,par=x$par,par2=x$par2,
    theta.G=paste(signif(x$theta.G,7),collapse=', '),logLik=x$log_lik_value,negLogLik=x$neg_log_lik_value,
    convergence=x$convergence,stringsAsFactors=FALSE)))
  summary<-summary[order(-summary$logLik),,drop=FALSE]; rownames(summary)<-NULL
  best_id<-summary$model_id[1L]

  .progress_msg(
    "\n[Best copula] ",
    summary$copula.family.name[1L],
    " (family ", summary$copula.family.lambda[1L], ")",
    " | logLik = ", signif(summary$logLik[1L], 8)
  )

  out<-list(method=method,d0=d0,marginal_pattern=marginal_pattern,zero_dependence_estimated=estimate_B_dependence,
    B.depen.option=if(estimate_B_dependence)B.depen.option else NA_character_,counting.depen.option=counting.depen.option,
    gaussian_prescreen=gaussian_fit,gaussian_rho=gaussian_rho,
    gaussian_dependence_sign=if(gaussian_rho>gaussian.sign.tol)'positive' else if(gaussian_rho < -gaussian.sign.tol)'negative' else 'near_zero',
    requested_copula_families=requested_families,fitted_copula_families=fitted_families,
    marginal_1=margin1,marginal_2=margin2,models=models,summary=summary,best_model=models[[best_id]])

  class(out) <- c("biv_zimp_fit", class(out))
  out
}


# =============================================================================
# PUBLIC BIVARIATE API
# =============================================================================

.bizimp_resolve_input <- function(formula1 = NULL, formula2 = NULL,
                                  data = NULL,
                                  Y = NULL, X = NULL, Z = NULL) {
  using_formula <- !is.null(formula1) || !is.null(formula2) || !is.null(data)
  using_direct  <- !is.null(Y) || !is.null(X) || !is.null(Z)

  if (using_formula && using_direct)
    stop("Use either 'formula1' + 'formula2' + 'data' OR 'Y' + 'X' (+ optional 'Z'), not both.")

  if (using_formula) {
    if (is.null(formula1) || is.null(formula2) || is.null(data))
      stop("'formula1', 'formula2', and 'data' must all be supplied.")
    if (!inherits(formula1, "formula"))
      stop("'formula1' must be a formula.")
    if (!inherits(formula2, "formula"))
      stop("'formula2' must be a formula.")
    if (!is.data.frame(data)) stop("'data' must be a data frame.")
    return(list(
      interface = "formula",
      formula1 = formula1,
      formula2 = formula2,
      formula = list(formula1, formula2),
      data = data
    ))
  }

  if (is.null(Y) || is.null(X))
    stop("Direct input requires both 'Y' and 'X'.")

  if (is.list(Y) && !is.data.frame(Y)) {
    if (length(Y) != 2L) stop("'Y' as a list must have exactly two elements.")
    yy <- Y
  } else {
    if (!(is.matrix(Y) || is.data.frame(Y)) || ncol(Y) != 2L)
      stop("'Y' must be an n x 2 matrix/data frame or a list of two vectors.")
    yy <- list(Y[,1L], Y[,2L])
  }

  check_y <- function(y,j) {
    if (!is.numeric(y) || !is.null(dim(y)))
      stop("Y[[",j,"]] must be a one-dimensional numeric vector.")
    if (!length(y) || any(!is.finite(y)) || any(y < 0) ||
        any(abs(y-round(y)) > 1e-8))
      stop("Y[[",j,"]] must contain finite non-negative integer counts.")
    as.numeric(y)
  }
  yy <- list(check_y(yy[[1]],1), check_y(yy[[2]],2))
  if (length(yy[[1]]) != length(yy[[2]]))
    stop("The two response margins must have the same number of observations.")
  n <- length(yy[[1]])

  split_design <- function(A, name, allow_null=FALSE) {
    if (is.null(A)) {
      if (allow_null) return(list(NULL,NULL))
      stop("'",name,"' cannot be NULL.")
    }
    ans <- if (is.list(A) && !is.data.frame(A)) {
      if (length(A) != 2L) stop("'",name,"' as a list must have exactly two elements.")
      A
    } else list(A,A)

    for (j in 1:2) {
      if (is.null(ans[[j]])) {
        if (allow_null) next
        stop(name,"[[",j,"]] cannot be NULL.")
      }
      if (!(is.matrix(ans[[j]]) || is.data.frame(ans[[j]])))
        stop(name,"[[",j,"]] must be a matrix or data frame.")
      if (nrow(ans[[j]]) != n) stop(name,"[[",j,"]] must have ",n," rows.")
      if (ncol(ans[[j]]) < 1L) stop(name,"[[",j,"]] must have at least one column.")
      if (!all(vapply(as.data.frame(ans[[j]]), is.numeric, logical(1))))
        stop(name,"[[",j,"]] must contain numeric columns only.")
      if (any(!is.finite(as.matrix(ans[[j]]))))
        stop(name,"[[",j,"]] contains NA/NaN/Inf.")
    }
    ans
  }

  list(interface="direct", Y=yy,
       X=split_design(X,"X",FALSE),
       Z=split_design(Z,"Z",TRUE),
       formula1=NULL, formula2=NULL, formula=NULL, data=NULL)
}


.bizimp_fmt_num <- function(x, digits = 6L) {
  if (length(x) == 0L || is.null(x)) return("")
  paste(format(signif(as.numeric(x), digits), trim = TRUE), collapse = ", ")
}

.bizimp_print_margin_progress <- function(j, fit) {
  m <- .extract_zimp_margin(fit)

  cat("
------------------------------------------------------------
")
  cat("[Margin ", j, "] ",
      if (m$zero_inflation) "ZIMP" else "Mixed Poisson",
      " | mixing = ", m$distribution,
      paste0(" | var.power = ",
             if (!is.null(m$var_power) && length(m$var_power) == 1L &&
                 is.finite(m$var_power)) {
               m$var_power
             } else {
               defaults <- c(
                 gamma = 2, lnorm = 2, ig = 3,
                 weibull = 2, loglogistic = 2, betaprime = 2
               )
               unname(defaults[[m$distribution]])
             },
             if (is.null(m$var_power) || !length(m$var_power) ||
                 !is.finite(m$var_power)) " (default)" else ""),
      "
", sep = "")

  .print_horizontal <- function(title, nm, val) {
    if (is.null(nm) || length(nm) != length(val))
      nm <- paste0("par", seq_along(val))

    vals <- format(signif(as.numeric(val), 6), trim = TRUE)
    widths <- pmax(nchar(nm), nchar(vals), 8L)

    cat(title, "
", sep = "")
    cat("  ")
    for (i in seq_along(nm))
      cat(format(nm[i], width = widths[i], justify = "right"),
          if (i < length(nm)) "  " else "
", sep = "")
    cat("  ")
    for (i in seq_along(vals))
      cat(format(vals[i], width = widths[i], justify = "right"),
          if (i < length(vals)) "  " else "
", sep = "")
  }

  bn <- colnames(m$X)
  if (is.null(bn) || length(bn) != length(m$be))
    bn <- paste0("beta", seq_along(m$be))
  .print_horizontal("Count coefficients:", c(bn, "phi"), c(m$be, m$phi))

  if (m$zero_inflation) {
    an <- colnames(m$Z)
    if (is.null(an) || length(an) != length(m$al))
      an <- paste0("alpha", seq_along(m$al))
    .print_horizontal("Zero coefficients:", an, m$al)
  }

  flush.console()
}

.bizimp_public_copula_name <- function(family, fallback = NULL) {
  family <- as.integer(family)

  nm <- c(
    `0` = "Independence",
    `1` = "Gaussian",
    `2` = "Student-t",
    `3` = "Clayton",
    `4` = "Gumbel",
    `5` = "Frank",
    `6` = "Joe",
    `7` = "BB1",
    `8` = "BB6",
    `9` = "BB7",
    `10` = "BB8"
  )

  key <- as.character(family)
  if (key %in% names(nm))
    return(unname(nm[key]))

  # Standard VineCopula rotated-family codes.
  rot <- list(`10` = "180", `20` = "90", `30` = "270")
  for (shift in c(10L, 20L, 30L)) {
    base_family <- family - shift
    base_key <- as.character(base_family)
    if (base_family >= 3L && base_key %in% names(nm)) {
      return(paste0(unname(nm[base_key]), " ",
                    rot[[as.character(shift)]], "\u00b0"))
    }
  }

  if (!is.null(fallback) && length(fallback) == 1L &&
      !is.na(fallback) && nzchar(fallback))
    return(as.character(fallback))

  paste0("Family ", family)
}

.bizimp_display_copula_name <- function(model) {
  fam <- as.integer(model$family)
  if (fam == 2L) {
    df <- suppressWarnings(as.numeric(model$par2)[1L])
    if (is.finite(df)) {
      return(paste0(
        "Student-t (",
        format(round(df, 3L), nsmall = 0, trim = TRUE),
        ")"
      ))
    }
    return("Student-t")
  }
  .bizimp_public_copula_name(fam, fallback = model$family_name)
}

.bizimp_display_var_power <- function(fit) {
  vp <- fit$var.power
  if (!is.null(vp) && length(vp) && is.finite(vp[1L])) return(as.numeric(vp[1L]))
  defaults <- c(
    gamma = 2, lnorm = 2, ig = 3,
    weibull = 2, loglogistic = 2, betaprime = 2
  )
  dist <- as.character(fit$mixing.distribution)[1L]
  if (!is.na(dist) && dist %in% names(defaults)) return(unname(defaults[[dist]]))
  NA_real_
}



.bizimp_margin_npar <- function(fit) {
  if (!is.null(fit$par) && is.numeric(fit$par)) {
    return(length(fit$par))
  }
  p <- if (!is.null(fit$X)) ncol(as.matrix(fit$X)) else 0L
  q <- if (isTRUE(fit$zero_inflation) && !is.null(fit$Z)) {
    ncol(as.matrix(fit$Z))
  } else {
    0L
  }
  as.integer(p + 1L + q)  # beta + phi + gamma
}

.bizimp_mixing_npar <- function(model) {
  fam <- as.integer(model$family)
  if (fam == 0L) return(0L)
  if (fam == 2L) return(2L)  # rho + df

  th <- model$theta.G
  if (!is.null(th) && length(th)) {
    return(length(as.numeric(th)))
  }

  # Fallback to the VineCopula family specification.
  sp <- try(.vine_family_spec(fam, gaussian.rho = 0), silent = TRUE)
  if (!inherits(sp, "try-error") && !is.null(sp$npar)) {
    return(as.integer(sp$npar))
  }
  1L
}

.bizimp_total_npar <- function(margin1, margin2, model,
                               zero_option = NULL) {
  k <- .bizimp_margin_npar(margin1) +
       .bizimp_margin_npar(margin2) +
       .bizimp_mixing_npar(model)

  both_zero <- isTRUE(margin1$zero_inflation) &&
               isTRUE(margin2$zero_inflation)

  # Gaussian zero copula estimates rho. Frechet upper/lower are fixed.
  if (both_zero && !is.null(zero_option) &&
      identical(zero_option, "normal")) {
    k <- k + 1L
  }
  as.integer(k)
}

.bizimp_information_criteria <- function(logLik, n, k) {
  if (!is.finite(logLik) || !is.finite(n) || n <= 0L ||
      !is.finite(k) || k < 0L) {
    return(c(AIC = Inf, BIC = Inf))
  }
  c(
    AIC = -2 * logLik + 2 * k,
    BIC = -2 * logLik + log(n) * k
  )
}


.bizimp_print_dependence_progress <- function(out, d0) {
  bm <- out$best_model

  cat("\n------------------------------------------------------------\n")
  cat("[Dependence]\n")

  # The continuous copula is discretized to the d x d checkerboard
  # density used by the likelihood.
  cat("Mixing copula:\n")
  copula_name <- if (bm$family == 2L) {
    paste0("Student-t (", .bizimp_fmt_num(bm$par2), ")")
  } else {
    .bizimp_public_copula_name(bm$family, fallback = bm$family_name)
  }
  cat("  ", copula_name, " checkerboard copula",
      " | d = ", d0, "\n", sep = "")

  if (bm$family == 0L) {
    cat("  parameter = none (independence)\n")
  } else if (bm$family == 2L) {
    cat("  rho = ", .bizimp_fmt_num(bm$par), "\n", sep = "")
  } else {
    th <- bm$theta.G
    if (is.null(th))
      th <- if (bm$family %in% c(
        7L,8L,9L,10L,17L,18L,19L,20L,27L,28L,29L,30L,
        37L,38L,39L,40L,104L,114L,124L,134L,204L,214L,224L,234L
      )) c(bm$par, bm$par2) else bm$par
    cat("  parameter = ", .bizimp_fmt_num(th), "\n", sep = "")
  }

  cat("Zero copula:\n")
  if (!isTRUE(out$zero_dependence_estimated)) {
    cat("  not active (both margins must be zero-inflated)\n")
  } else if (identical(out$B.depen.option, "normal")) {
    cat("  Gaussian copula | rho = ",
        .bizimp_fmt_num(bm$theta.B), "\n", sep = "")
  } else if (identical(out$B.depen.option, "comonotone")) {
    cat("  Frechet upper\n")
  } else if (identical(out$B.depen.option, "countermonotone")) {
    cat("  Frechet lower\n")
  }

  ll_display <- if (!is.finite(bm$log_lik_value) ||
                    bm$log_lik_value <= -1e12 + 1) {
    "-Inf"
  } else {
    format(signif(bm$log_lik_value, 8), trim = TRUE)
  }

  zero_option <- if (isTRUE(out$zero_dependence_estimated)) {
    out$B.depen.option
  } else {
    NA_character_
  }
  k_total <- .bizimp_total_npar(
    out$marginal_1, out$marginal_2, bm, zero_option
  )
  ic <- .bizimp_information_criteria(
    bm$log_lik_value, length(out$marginal_1$y), k_total
  )

  cat("Log-likelihood = ", ll_display, "\n", sep = "")
  cat("AIC = ", format(round(ic[["AIC"]], 2), nsmall = 2, trim = TRUE),
      " | BIC = ", format(round(ic[["BIC"]], 2), nsmall = 2, trim = TRUE),
      "\n", sep = "")
  cat("============================================================\n")
  flush.console()
}



#' Fit a Bivariate Mixed Poisson or Zero-Inflated Mixed Poisson Model
#'
#' @description
#' Fits a bivariate Mixed Poisson (MP) or Zero-Inflated Mixed Poisson
#' (ZIMP) regression model using inference functions for margins (IFM).
#'
#' The two marginal models are first fitted separately using [fit_zimp()].
#' Dependence between their latent mixing variables is then modeled by a
#' specified bivariate copula. When both margins are zero-inflated,
#' dependence between the Bernoulli activation indicators may also be
#' modeled.
#'
#' @details
#' `fit_bizimp()` fits one fully specified bivariate model. It does not
#' perform model selection over marginal distributions or dependence
#' structures.
#'
#' Let the two observed counts be
#'
#' \deqn{Y_j = B_j N_j,\qquad j=1,2,}
#'
#' where \eqn{N_j} follows a Mixed Poisson model and \eqn{B_j} is an
#' optional Bernoulli activation indicator. For an MP margin,
#' \eqn{B_j = 1} identically.
#'
#' Conditional on the latent intensity,
#'
#' \deqn{N_j \mid \Lambda_j \sim
#'       \mathrm{Poisson}(\Lambda_j),}
#'
#' with marginal mean
#'
#' \deqn{\mu_{ij} = \exp(x_{ij}^{\top}\beta_j).}
#'
#' For a zero-inflated margin,
#'
#' \deqn{\Pr(B_{ij}=1)=\pi_{ij},\qquad
#'       \mathrm{logit}(\pi_{ij})=z_{ij}^{\top}\gamma_j.}
#'
#' Thus, \eqn{B_{ij}=1} indicates that observation \eqn{i} belongs to the
#' count-generating component of margin \eqn{j}.
#'
#' Dependence between \eqn{\Lambda_1} and \eqn{\Lambda_2} is represented by
#' the copula specified by `mixing.copula`. The copula is evaluated through
#' a checkerboard approximation with grid size `d0`.
#'
#' When both margins are zero-inflated, dependence between
#' \eqn{B_1} and \eqn{B_2} is specified separately through
#' `zero.dependence`.
#'
#' @section Estimation:
#'
#' Estimation follows the IFM principle:
#'
#' \enumerate{
#'   \item Fit margin 1 using [fit_zimp()].
#'   \item Fit margin 2 independently using [fit_zimp()].
#'   \item Conditional on the fitted marginal parameters, estimate the
#'     dependence parameter(s) from the bivariate likelihood.
#' }
#'
#' Consequently, the marginal estimates returned in `fit_margin_1` and
#' `fit_margin_2` are the corresponding independent [fit_zimp()] estimates
#' and are not re-estimated during dependence fitting.
#'
#' `fit_bizimp()` fits exactly one requested mixing-copula family. When the
#' requested family is not Gaussian, a Gaussian dependence model may be
#' fitted internally as a pilot for initialization and dependence-sign
#' information; it is not an additional candidate model.
#'
#' @section Marginal configurations:
#'
#' Each margin may independently be an MP or ZIMP model. The four supported
#' configurations are therefore:
#'
#' \itemize{
#'   \item MP--MP;
#'   \item ZIMP--MP;
#'   \item MP--ZIMP;
#'   \item ZIMP--ZIMP.
#' }
#'
#' For direct input, `Z = NULL` specifies MP--MP. A list such as
#' `Z = list(Z1, NULL)` specifies ZIMP--MP, while
#' `Z = list(NULL, Z2)` specifies MP--ZIMP. Supplying a shared matrix
#' `Z` applies the zero-component design to both margins.
#'
#' Dependence between the Bernoulli activation indicators is estimated only
#' for the ZIMP--ZIMP configuration.
#'
#' @section Model specification:
#'
#' **Formula interface**
#'
#' Two formulas are supplied, one for each margin. For example,
#'
#' ```
#' formula1 = y1 ~ 1 + x1 + x2
#' formula2 = y2 ~ 1 + x1
#' ```
#'
#' specifies two MP margins. Zero-inflated margins are specified using
#' two-part formulas:
#'
#' ```
#' formula1 = y1 ~ 1 + x1 + x2 | z1 + z2
#' formula2 = y2 ~ 1 + x1      | z1
#' ```
#'
#' The expression before `|` specifies the count component and the
#' expression after `|` specifies the zero-inflation component.
#'
#' As in [fit_zimp()], intercepts are included only when explicitly
#' requested with `1`. This convention applies separately to the count and
#' zero-inflation components.
#'
#' **Matrix interface**
#'
#' `Y` may be an \eqn{n \times 2} matrix/data frame or a list containing
#' two response vectors.
#'
#' `X` may be a common design matrix for both margins or a list
#' `list(X1, X2)` containing margin-specific count designs. Similarly, `Z`
#' may be a common zero-component design or a list containing
#' margin-specific designs. Intercept columns must be included explicitly
#' when required.
#'
#' @section Mixing distributions:
#'
#' The mixing distribution is specified separately for each margin through
#' `mixing.distribution.1` and `mixing.distribution.2`. Supported
#' distributions are:
#'
#' \describe{
#'   \item{`"gamma"`}{Gamma mixing distribution.}
#'   \item{`"lnorm"`}{Lognormal mixing distribution.}
#'   \item{`"ig"`}{Inverse Gaussian mixing distribution.}
#'   \item{`"weibull"`}{Weibull mixing distribution.}
#'   \item{`"loglogistic"`}{Log-logistic mixing distribution.}
#'   \item{`"betaprime"`}{Beta-prime mixing distribution.}
#' }
#'
#' The two margins may use different mixing distributions.
#'
#' @section Dependence structure:
#'
#' `mixing.copula` is a single integer family code supported by
#' `VineCopula`. It specifies the copula connecting the two latent mixing
#' variables. The default `mixing.copula = 1` corresponds to the Gaussian
#' copula.
#'
#' For `mixing.copula = 2`, corresponding to the Student-t copula, the
#' degrees of freedom are selected over the values supplied in `t.df.grid`.
#'
#' When both margins are zero-inflated, `zero.dependence` specifies the
#' dependence structure between their Bernoulli activation indicators:
#'
#' \describe{
#'   \item{`"gaussian"`}{
#'     Gaussian copula dependence with an estimated correlation parameter.
#'   }
#'   \item{`"comonotone"`}{
#'     Upper Fréchet-bound dependence.
#'   }
#'   \item{`"countermonotone"`}{
#'     Lower Fréchet-bound dependence.
#'   }
#' }
#'
#' @param formula1 Optional formula specifying the marginal model for the
#'   first response. Use a one-part formula for an MP margin and a two-part
#'   formula of the form `y1 ~ ... | ...` for a ZIMP margin.
#'
#' @param formula2 Optional formula specifying the marginal model for the
#'   second response. It follows the same syntax as `formula1`.
#'
#' @param data Data frame containing the variables referenced by `formula1`
#'   and `formula2`. When the formula interface is used, `formula1`,
#'   `formula2`, and `data` must all be supplied.
#'
#' @param Y Response for the direct matrix interface. Either an
#'   \eqn{n \times 2} numeric matrix/data frame or a list of exactly two
#'   numeric vectors. Responses must contain finite, non-negative integer
#'   counts.
#'
#' @param X Count-component design for the direct interface. May be a
#'   numeric matrix/data frame shared by both margins or a list of exactly
#'   two margin-specific design matrices. Include intercept columns
#'   explicitly when required.
#'
#' @param Z Optional zero-component design for the direct interface. May be
#'   `NULL`, a numeric matrix/data frame shared by both margins, or a list
#'   of exactly two margin-specific designs. A `NULL` list element specifies
#'   an MP margin; a non-`NULL` element specifies a ZIMP margin.
#'
#' @param mixing.distribution.1 Character string specifying the mixing
#'   distribution for margin 1. One of `"gamma"`, `"lnorm"`, `"ig"`,
#'   `"weibull"`, `"loglogistic"`, or `"betaprime"`.
#'
#' @param mixing.distribution.2 Character string specifying the mixing
#'   distribution for margin 2. One of `"gamma"`, `"lnorm"`, `"ig"`,
#'   `"weibull"`, `"loglogistic"`, or `"betaprime"`.
#'
#' @param var.power Optional variance-power specification for the two
#'   marginal models. If `NULL`, each margin uses the native/default
#'   parameterization of its mixing distribution. A scalar value from
#'   `1`, `2`, or `3` is applied to both margins. A vector of length two
#'   specifies the variance power separately for margins 1 and 2.
#'
#' @param mixing.copula Integer `VineCopula` family code specifying exactly
#'   one copula family for dependence between the two latent mixing
#'   variables. Supported family codes are:
#'
#'   \describe{
#'     \item{`0L`}{Independence copula.}
#'     \item{`1L`}{Gaussian copula (default).}
#'     \item{`2L`}{Student-t copula.}
#'     \item{`3L`}{Clayton copula.}
#'     \item{`4L`}{Gumbel copula.}
#'     \item{`5L`}{Frank copula.}
#'     \item{`6L`}{Joe copula.}
#'     \item{`7L`}{BB1 copula.}
#'     \item{`8L`}{BB6 copula.}
#'     \item{`9L`}{BB7 copula.}
#'     \item{`10L`}{BB8 copula.}
#'     \item{`13L`, `14L`, `16L`--`20L`}{
#'       180-degree rotations of Clayton, Gumbel, Joe, and BB copulas.
#'     }
#'     \item{`23L`, `24L`, `26L`--`30L`}{
#'       90-degree rotations of Clayton, Gumbel, Joe, and BB copulas.
#'     }
#'     \item{`33L`, `34L`, `36L`--`40L`}{
#'       270-degree rotations of Clayton, Gumbel, Joe, and BB copulas.
#'     }
#'     \item{`104L`, `114L`, `124L`, `134L`,
#'           `204L`, `214L`, `224L`, `234L`}{
#'       Tawn copula families and their rotations.
#'     }
#'   }
#'
#'   Family codes follow the convention used by the `VineCopula` package.
#'   The default is `1L`, corresponding to the Gaussian copula.
#'   `fit_bizimp()` fits exactly one requested family; use
#'   [select_bizimp()] to compare multiple copula families.
#'
#' @param zero.dependence Character string specifying dependence between
#'   the Bernoulli activation indicators when both margins are
#'   zero-inflated. One of `"gaussian"`, `"comonotone"`, or
#'   `"countermonotone"`. This argument is inactive unless both margins
#'   are zero-inflated.
#'
#' @param d0 Positive integer specifying the checkerboard grid size used to
#'   approximate the mixing copula. The default is `10L`.
#'
#' @param n_gh Positive integer specifying the number of Gauss--Hermite
#'   quadrature nodes used for lognormal marginal mixing. The default is
#'   `20`.
#'
#' @param gh Optional precomputed Gauss--Hermite quadrature object. It is
#'   passed to the lognormal marginal calculations. If `NULL`, the required
#'   quadrature rule is constructed internally.
#'
#' @param ini.par Optional list of length two containing initial parameter
#'   vectors for the two marginal [fit_zimp()] fits. `NULL` uses the default
#'   initialization procedure for both margins.
#'
#' @param progress Logical indicating whether fitting progress should be
#'   displayed. The default is `FALSE`.
#'
#' @param trace Numeric optimizer trace level. The default is `0`.
#'
#' @param eval.max Positive integer giving the maximum number of objective
#'   function evaluations. The default is `100`.
#'
#' @param rel.tol Positive numeric value giving the relative convergence
#'   tolerance used in marginal fitting. The default is `1e-5`.
#'
#' @param x.tol Positive numeric value giving the parameter convergence
#'   tolerance used in marginal fitting. The default is `1e-5`.
#'
#' @param iter.max Positive integer giving the maximum number of optimizer
#'   iterations. The default is `100`.
#'
#' @param gaussian.sign.tol Positive numeric tolerance used internally to
#'   classify the sign of the Gaussian pilot dependence estimate. The
#'   default is `1e-6`.
#'
#' @param t.df.grid Integer vector of candidate degrees of freedom greater
#'   than two used when `mixing.copula = 2`, corresponding to the Student-t
#'   copula. The default is `c(3, 5, 10, 30)`.
#'
#' @return
#' An object of class `"bizimp_fit"` containing the fitted marginal models,
#' dependence model, likelihood, and information criteria. Important
#' components include:
#'
#' \describe{
#'   \item{`fit_margin_1`, `fit_margin_2`}{
#'     Fitted univariate marginal models returned by [fit_zimp()].
#'   }
#'   \item{`marginal_1`, `marginal_2`}{
#'     Internal representations of the fitted marginal models used in the
#'     dependence likelihood.
#'   }
#'   \item{`best_model`}{
#'     Fitted dependence model. It contains the requested mixing-copula
#'     family and estimated copula parameter(s), the zero-dependence
#'     parameter when applicable, the maximized bivariate log-likelihood,
#'     convergence information, and the full-model AIC and BIC.
#'   }
#'   \item{`mixing.distribution.1`, `mixing.distribution.2`}{
#'     Mixing distributions used for margins 1 and 2.
#'   }
#'   \item{`var.power`}{
#'     User-supplied marginal variance-power specification, or `NULL` when
#'     native/default parameterizations were used.
#'   }
#'   \item{`mixing.copula`}{
#'     `VineCopula` family code used for latent mixing dependence.
#'   }
#'   \item{`zero.dependence`}{
#'     Zero-dependence specification when both margins are zero-inflated;
#'     otherwise `NA`.
#'   }
#'   \item{`zero.copula.rho`}{
#'     Estimated Gaussian zero-copula correlation when both margins are
#'     zero-inflated and `zero.dependence = "gaussian"`; otherwise `NA`.
#'   }
#'   \item{`zero_inflation`}{
#'     Logical vector of length two indicating whether each margin is
#'     zero-inflated.
#'   }
#'   \item{`npar`}{
#'     Total number of parameters used to calculate the full-model
#'     information criteria.
#'   }
#'   \item{`AIC`, `BIC`}{
#'     Akaike and Bayesian information criteria for the complete bivariate
#'     model.
#'   }
#'   \item{`interface`}{
#'     Character string identifying whether the formula or direct matrix
#'     interface was used.
#'   }
#'   \item{`formula1`, `formula2`}{
#'     Marginal formulas when the formula interface was used.
#'   }
#'   \item{`call`}{
#'     Matched call to `fit_bizimp()`.
#'   }
#' }
#'
#' @seealso
#' [fit_zimp()] for fitting a univariate marginal model and
#' [select_bizimp()] for selecting among multiple bivariate model
#' specifications.
#'
#' @examples
#' ## Generate a small bivariate Mixed Poisson data set
#' set.seed(12345)
#' n <- 50
#' x <- rnorm(n)
#' X <- cbind(1, x)
#'
#' mu1 <- exp(X %*% c(0.4, 0.3))
#' mu2 <- exp(X %*% c(0.2, -0.2))
#'
#' ## Gamma-mixed Poisson margins
#' y1 <- rnbinom(n, size = 1 / 0.2, mu = mu1)
#' y2 <- rnbinom(n, size = 1 / 0.3, mu = mu2)
#' Y <- cbind(y1, y2)
#'
#' ## MP--MP model with a Gaussian mixing copula
#' fit <- fit_bizimp(
#'   Y = Y,
#'   X = X,
#'   mixing.distribution.1 = "gamma",
#'   mixing.distribution.2 = "gamma",
#'   var.power = 2,
#'   mixing.copula = 1L,
#'   progress = FALSE
#' )
#'
#' fit$best_model
#' fit$AIC
#' fit$BIC
#'
#' ## Margin-specific designs and a ZIMP--MP configuration
#' z1 <- rnorm(n)
#' Z1 <- cbind(z1)
#'
#' fit_zm <- fit_bizimp(
#'   Y = Y,
#'   X = list(X, X),
#'   Z = list(Z1, NULL),
#'   mixing.distribution.1 = "gamma",
#'   mixing.distribution.2 = "gamma",
#'   var.power = c(2, 2),
#'   mixing.copula = 1L,
#'   progress = FALSE
#' )
#'
#' fit_zm$zero_inflation
#'
#' @export

fit_bizimp <- function(
    formula1 = NULL, formula2 = NULL, data = NULL,
    Y = NULL, X = NULL, Z = NULL,
    mixing.distribution.1 = "gamma",
    mixing.distribution.2 = "gamma",
    var.power = NULL,
    mixing.copula = 1L,
    zero.dependence = c("gaussian", "comonotone", "countermonotone"),
    d0 = 10L,
    n_gh = 20, gh = NULL,
    ini.par = NULL,
    progress = FALSE, trace = 0,
    eval.max = 100, rel.tol = 1e-5,
    x.tol = 1e-5, iter.max = 100,
    gaussian.sign.tol = 1e-6,
    t.df.grid = c(3, 5, 10, 30)) {

  inp <- .bizimp_resolve_input(formula1, formula2, data, Y, X, Z)
  zero.dependence <- match.arg(zero.dependence)

  if (isTRUE(progress)) {
    cat("\n============================================================\n")
    cat("Fitting bivariate zero-inflated mixed-Poisson model\n")
    cat("============================================================\n")
  }

  recycle2 <- function(x, name) {
    if (length(x) == 1L) x <- rep(x, 2L)
    if (length(x) != 2L)
      stop("'", name, "' must have length 1 or 2.")
    x
  }

  allowed_dist <- c(
    "gamma", "lnorm", "ig",
    "weibull", "loglogistic", "betaprime"
  )
  if (length(mixing.distribution.1) != 1L ||
      !is.character(mixing.distribution.1) ||
      !(mixing.distribution.1 %in% allowed_dist)) {
    stop("Unsupported 'mixing.distribution.1': ",
         paste(mixing.distribution.1, collapse = ", "), ".")
  }

  if (length(mixing.distribution.2) != 1L ||
      !is.character(mixing.distribution.2) ||
      !(mixing.distribution.2 %in% allowed_dist)) {
    stop("Unsupported 'mixing.distribution.2': ",
         paste(mixing.distribution.2, collapse = ", "), ".")
  }

  mixing.distribution.by.margin <- c(
    mixing.distribution.1,
    mixing.distribution.2
  )

  # var.power = NULL -> native/default parameterization for both margins.
  # Scalar p -> same p for both; length 2 -> margin-specific p.
  if (is.null(var.power)) {
    vp <- list(NULL, NULL)
  } else {
    if (!is.numeric(var.power) || any(!is.finite(var.power)) ||
        any(abs(var.power - round(var.power)) > 1e-12) ||
        any(!(as.integer(round(var.power)) %in% c(1L, 2L, 3L)))) {
      stop("'var.power' must be NULL or contain values from {1, 2, 3}.")
    }
    var.power <- as.integer(round(var.power))
    if (length(var.power) == 1L) var.power <- rep(var.power, 2L)
    if (length(var.power) != 2L)
      stop("'var.power' must be NULL, a scalar, or a vector of length 2.")
    vp <- as.list(var.power)
  }

  if (length(mixing.copula) != 1L ||
      !is.numeric(mixing.copula) || !is.finite(mixing.copula)) {
    stop("fit_bizimp() fits one configuration; 'mixing.copula' ",
         "must be one VineCopula family code.")
  }
  mixing.copula <- as.integer(mixing.copula)

  if (is.null(ini.par)) ini.par <- list(NULL, NULL)
  if (!is.list(ini.par) || length(ini.par) != 2L)
    stop("'ini.par' must be NULL or a list of length 2.")

  fit_one <- function(j) {
    args <- list(
      mixing.distribution = mixing.distribution.by.margin[j],
      var.power = vp[[j]],
      n_gh = n_gh,
      gh = gh,
      ini.par = ini.par[[j]],
      progress = FALSE,
      trace = trace,
      eval.max = eval.max,
      rel.tol = rel.tol,
      x.tol = x.tol,
      iter.max = iter.max
    )

    if (inp$interface == "formula") {
      do.call(
        fit_zimp,
        c(list(formula = if (j == 1L) inp$formula1 else inp$formula2, data = inp$data), args)
      )
    } else {
      do.call(
        fit_zimp,
        c(list(
          Y = inp$Y[[j]],
          X = inp$X[[j]],
          Z = inp$Z[[j]]
        ), args)
      )
    }
  }

  if (isTRUE(progress)) cat("\n[1/3] Fitting margin 1 ...\n")
  f1 <- fit_one(1L)
  if (isTRUE(progress)) .bizimp_print_margin_progress(1L, f1)

  if (isTRUE(progress)) cat("\n[2/3] Fitting margin 2 ...\n")
  f2 <- fit_one(2L)
  if (isTRUE(progress)) .bizimp_print_margin_progress(2L, f2)

  if (length(f1$y) != length(f2$y))
    stop("The two margins must use the same number of observations.")

  if (inp$interface == "formula" &&
      !identical(f1$na_removed, f2$na_removed)) {
    stop("The two formulas removed different rows; ",
         "bivariate fitting requires identical observations.")
  }

  both_zero <- isTRUE(f1$zero_inflation) &&
               isTRUE(f2$zero_inflation)

  # Legacy names are kept inside the old dependence engine only.
  zero_internal <- switch(
    zero.dependence,
    gaussian = "normal",
    comonotone = "comonotone",
    countermonotone = "countermonotone"
  )

  if (isTRUE(progress))
    cat("\n[3/3] Fitting dependence structure ...\n")

  out <- .fit_bizimp_dependence_engine(
    fit_margin_1 = f1,
    fit_margin_2 = f2,
    n_gh = n_gh,
    gh = gh,
    exact.copula.family = mixing.copula,
    B.depen.option = zero_internal,
    d0 = d0,
    gaussian.sign.tol = gaussian.sign.tol,
    t.df.grid = t.df.grid,
    progress = FALSE,
    trace = trace,
    eval.max = eval.max,
    iter.max = iter.max
  )

  # Full-model information criteria: both marginal regressions/dispersion
  # parameters + zero-copula parameter (when estimated) + mixing-copula
  # parameter(s).
  k_total <- .bizimp_total_npar(
    out$marginal_1, out$marginal_2, out$best_model,
    if (isTRUE(out$zero_dependence_estimated)) out$B.depen.option else NA_character_
  )
  ic <- .bizimp_information_criteria(
    out$best_model$log_lik_value, length(out$marginal_1$y), k_total
  )
  out$npar <- k_total
  out$AIC <- unname(ic[["AIC"]])
  out$BIC <- unname(ic[["BIC"]])
  out$best_model$npar <- k_total
  out$best_model$AIC <- out$AIC
  out$best_model$BIC <- out$BIC

  if (isTRUE(progress))
    .bizimp_print_dependence_progress(out, d0)

  out$call <- match.call()
  out$interface <- inp$interface
  out$formula1 <- inp$formula1
  out$formula2 <- inp$formula2
  out$formula <- if (inp$interface == "formula") list(inp$formula1, inp$formula2) else NULL
  out$fit_margin_1 <- f1
  out$fit_margin_2 <- f2

  # Public-facing terminology.
  out$mixing.distribution.1 <- mixing.distribution.1
  out$mixing.distribution.2 <- mixing.distribution.2
  out$var.power <- if (is.null(var.power)) NULL else var.power
  out$mixing.copula <- mixing.copula
  out$zero.dependence <- if (both_zero) zero.dependence else NA_character_
  out$zero.copula.rho <- if (!both_zero) {
    NA_real_
  } else if (zero.dependence == "gaussian") {
    out$best_model$theta.B
  } else {
    NA_real_
  }
  out$zero_inflation <- c(f1$zero_inflation, f2$zero_inflation)

  class(out) <- unique(c("bizimp_fit", class(out)))
  out
}
