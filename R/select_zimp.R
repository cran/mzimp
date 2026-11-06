# =============================================================================
# Univariate Poisson-mixture / zero-inflated Poisson-mixture models
#
# PUBLIC API
#   fit_zimp()    : fit exactly one specified marginal model.
#   select_zimp() : fit a set of candidate marginal models and return the best.
#
# INTERNAL HELPERS
#   Every other top-level function in this file begins with "." and is not
#   exported. These helpers parse formulas, initialize parameters, run the
#   numerical optimizer, and compute post-estimation inference.
#
# CRAN / roxygen2
#   Only fit_zimp() and select_zimp() below carry @export tags. Therefore,
#   roxygen2::roxygenise() will expose only these two user-facing functions
#   from this file (apart from any S3 methods exported elsewhere).
# =============================================================================

# Internal formula helper: test whether an expression is exactly a given finite number.
.zimp_is_number <- function(x, value) {
  is.numeric(x) &&
    length(x) == 1L &&
    is.finite(x) &&
    identical(as.numeric(x), as.numeric(value))
}

.zimp_has_explicit_intercept <- function(expr) {
  if (.zimp_is_number(expr, 1)) return(TRUE)
  if (!is.call(expr)) return(FALSE)

  op <- as.character(expr[[1L]])

  if (identical(op, "+")) {
    if (length(expr) == 2L) {
      return(.zimp_has_explicit_intercept(expr[[2L]]))
    }
    return(
      .zimp_has_explicit_intercept(expr[[2L]]) ||
        .zimp_has_explicit_intercept(expr[[3L]])
    )
  }

  FALSE
}

.zimp_contains_disallowed_intercept_syntax <- function(expr) {
  if (!is.call(expr)) return(FALSE)

  op <- as.character(expr[[1L]])

  if (identical(op, "-") &&
      length(expr) == 2L &&
      .zimp_is_number(expr[[2L]], 1)) {
    return(TRUE)
  }

  if (identical(op, "-") &&
      length(expr) == 3L &&
      (.zimp_is_number(expr[[3L]], 1) ||
       .zimp_is_number(expr[[3L]], 0))) {
    return(TRUE)
  }

  if (identical(op, "+")) {
    if (length(expr) == 2L) {
      return(.zimp_contains_disallowed_intercept_syntax(expr[[2L]]))
    }

    if (.zimp_is_number(expr[[2L]], 0) ||
        .zimp_is_number(expr[[3L]], 0)) {
      return(TRUE)
    }

    return(
      .zimp_contains_disallowed_intercept_syntax(expr[[2L]]) ||
        .zimp_contains_disallowed_intercept_syntax(expr[[3L]])
    )
  }

  FALSE
}

# Internal formula helper: convert a count/zero RHS expression to a model-matrix formula
# while enforcing zimp's explicit-intercept convention.
.zimp_make_rhs_formula <- function(rhs_expr, env, part = c("count", "zero")) {
  part <- match.arg(part)

  if (.zimp_contains_disallowed_intercept_syntax(rhs_expr)) {
    stop(
      "Do not use '-1', 'x - 1', or '0 +' to control the intercept in the ",
      part,
      " formula. In the zimp formula interface, the intercept is excluded by default. ",
      "Write an explicit '1 +' only when an intercept is required."
    )
  }

  has_intercept <- .zimp_has_explicit_intercept(rhs_expr)

  rhs_chr <- paste(
    deparse(rhs_expr, width.cutoff = 500L),
    collapse = ""
  )

  rhs_model_chr <- if (has_intercept) {
    rhs_chr
  } else {
    paste("0 +", rhs_chr)
  }

  list(
    formula = stats::as.formula(
      paste("~", rhs_model_chr),
      env = env
    ),
    has_intercept = has_intercept
  )
}

# Internal initialization helper: initialize beta, phi, and (when present) gamma
# using GLM fits plus distribution-specific moment calculations.
.zimp_init_par <- function(y, X, Z = NULL, offset.count = NULL,
                                 dist = c("gamma", "lnorm", "pareto",
                                          "ig", "weibull", "loglogistic", "betaprime"),
                                 zero_inflation = FALSE,
                                 var.power = NULL,
                                 eps = 1e-6) {
  dist <- match.arg(dist)
  
  if (!is.null(var.power) &&
      (length(var.power) != 1L || !is.numeric(var.power) || !is.finite(var.power) ||
       !(var.power %in% c(1, 2, 3)))) {
    stop("var.power must be NULL, 1, 2, or 3.")
  }
  use_custom_var_power <- !is.null(var.power)
  var_power <- if (is.null(var.power)) 2 else as.integer(var.power)
  
  X <- as.matrix(X)
  y <- as.numeric(y)
  
  n <- length(y)
  p <- ncol(X)

  if (is.null(offset.count)) offset.count <- rep(0, n)
  offset.count <- as.numeric(offset.count)
  if (length(offset.count) != n || any(!is.finite(offset.count)))
    stop("offset.count must be NULL or a finite numeric vector of length(y).")
  
  if (length(y) != nrow(X)) stop("length(y) must equal nrow(X).")
  if (is.null(offset.count)) offset.count <- rep(0, length(y))
  offset.count <- as.numeric(offset.count)
  if (length(offset.count) != length(y) || any(!is.finite(offset.count)))
    stop("offset.count must be NULL or a finite numeric vector of length(y).")
  
  if (zero_inflation) {
    Z <- as.matrix(Z)
    if (length(y) != nrow(Z)) stop("length(y) must equal nrow(Z).")
    p1 <- ncol(Z)
  } else {
    Z  <- NULL
    p1 <- 0
  }
  
  ## --------------------------------------------------
  ## 1. Khởi tạo gamma cho zero layer
  ## --------------------------------------------------
  if (zero_inflation) {
    y_pos <- as.numeric(y > 0)
    
    gam_init <- rep(0, p1)
    
    fit_gam <- tryCatch(
      suppressWarnings(glm(y_pos ~ Z - 1, family = binomial())),
      error = function(e) NULL
    )
    
    if (!is.null(fit_gam)) {
      cc <- coef(fit_gam)
      cc[!is.finite(cc)] <- 0
      if (length(cc) == p1) gam_init <- cc
    } else {
      p_pos <- min(1 - eps, max(eps, mean(y_pos)))
      gam_init <- rep(0, p1)
      if (p1 >= 1) gam_init[1] <- qlogis(p_pos)
    }
    
    pr_init <- safe_plogis(as.numeric(Z %*% gam_init))
    pr_init <- pmin(1 - eps, pmax(eps, pr_init))
    
  } else {
    gam_init <- NULL
    pr_init <- rep(1, n)
  }
  
  ## --------------------------------------------------
  ## 2. Khởi tạo beta bằng Poisson GLM
  ## --------------------------------------------------
  be_init <- rep(0, p)
  
  fit_be <- tryCatch(
    suppressWarnings(glm(y ~ X - 1, family = poisson(), offset = offset.count)),
    error = function(e) NULL
  )
  
  if (!is.null(fit_be)) {
    cc <- coef(fit_be)
    cc[!is.finite(cc)] <- 0
    if (length(cc) == p) be_init <- cc
  } else {
    mu_bar <- max(mean(y) / max(mean(pr_init), eps), eps)
    be_init <- rep(0, p)
    if (p >= 1) be_init[1] <- log(mu_bar)
  }
  
  mu_tmp <- safe_mu_from_Xbe(cbind(X, offset.count), c(be_init, 1))
  Ey_tmp <- mean(pr_init * mu_tmp)
  ybar   <- mean(y)
  
  if (p >= 1 &&
      all(abs(X[, 1] - 1) < 1e-12) &&
      is.finite(Ey_tmp) && Ey_tmp > eps) {
    be_init[1] <- be_init[1] + log(max(ybar, eps) / Ey_tmp)
  }
  
  mu_init <- safe_mu_from_Xbe(cbind(X, offset.count), c(be_init, 1))
  mu_init <- pmax(mu_init, eps)
  
  ## --------------------------------------------------
  ## 3. Moment quantities
  ## --------------------------------------------------
  Ey <- mean(y)
  Vy <- var(y)
  
  m_i   <- pr_init * mu_init
  extra <- if (zero_inflation) mean(pr_init * (1 - pr_init) * mu_init^2) else 0
  heter <- var(m_i)
  
  D <- Vy - Ey - extra - heter
  
  q95 <- as.numeric(stats::quantile(mu_init, 0.95, na.rm = TRUE))
  mu_trunc <- pmin(mu_init, q95)
  
  phi_init <- 1
  
  ## --------------------------------------------------
  ## 4. Khởi tạo phi theo dist và var.power
  ## --------------------------------------------------
  
  if (dist == "gamma") {
    
    p_gamma <- if (use_custom_var_power) {
      var_power
    } else {
      2
    }
    
    A <- mean(
      pr_init * mu_trunc^p_gamma
    )
    
    if (is.finite(D) && D > eps &&
        is.finite(A) && A > eps) {
      
      phi_init <- max(D / A, eps)
      
    } else {
      
      phi_init <- 1
    }
    
  } else if (dist == "lnorm") {
    
    if (use_custom_var_power && var_power == 1) {
      ## tương ứng lnorm1 cũ
      A <- mean(pr_init * mu_trunc)
      
      if (is.finite(D) && D > eps &&
          is.finite(A) && A > eps) {
        phi_init <- D / A
        phi_init <- min(max(phi_init, 0.3), 8)
      } else {
        phi_init <- 1
      }
      
    } else {
      ## lnorm gốc
      od <- if (is.finite(Ey) && Ey > eps) (Vy - Ey) / Ey else 1
      od <- max(od, eps)
      phi_init <- min(max(sqrt(od), 0.3), 3)
    }
  } else if (dist == "pareto") {

  A <- mean(
    pr_init * mu_trunc^2
  )

  if (is.finite(D) &&
      is.finite(A) &&
      D > eps &&
      A > eps) {

    ratio <- D / A

    if (is.finite(ratio) && ratio > 1) {
      phi_init <- (ratio - 1) / 2
    } else {
      phi_init <- 0.01
    }

    phi_init <- min(
      max(phi_init, 0.01),
      20
    )

  } else {

    phi_init <- 0.5
  }
 } else if (dist == "ig") {
    if (use_custom_var_power && var_power == 2) {
      ## tương ứng ig1 cũ
      A <- mean(pr_init * mu_trunc^2)
      
      if (is.finite(D) && is.finite(A) && A > eps && D > eps) {
        phi_init <- D / A
        phi_init <- min(max(phi_init, 0.3), 8)
      } else {
        phi_init <- 1
      }
      
    } else {
      ## ig gốc
      A <- mean(pr_init * mu_trunc^3)
      
      if (is.finite(D) && is.finite(A) && A > eps && D > eps) {
        phi_init <- D / A
        phi_init <- min(max(phi_init, 0.3), 8)
      } else {
        phi_init <- 1
      }
    }
    
  } else if (dist == "weibull") {
    
    phi_init <- 0.7
    
  } else if (dist == "loglogistic") {
    
    od <- if (is.finite(Ey) && Ey > eps) (Vy - Ey) / Ey else 1
    od <- max(od, eps)
    phi_init <- 1 / sqrt(od + eps)
    phi_init <- min(max(phi_init, 0.3), 3)
  

  } else if (dist == "betaprime") {

    # Symmetric beta-prime:
    # Var(Lambda_i) = phi * mu_i^p.

    p_betaprime <- if (use_custom_var_power) {
      var_power
    } else {
      2
    }

    A <- mean(
      pr_init * mu_trunc^p_betaprime
    )

    if (is.finite(D) &&
        is.finite(A) &&
        D > eps &&
        A > eps) {

      phi_init <- D / A

    } else {

      phi_init <- 0.5
    }

    phi_init <- min(
      max(phi_init, 0.01),
      20
    )
  }
  
  phi_init <- max(phi_init, eps)
  
  if (zero_inflation) {
    list(
      be  = be_init,
      phi = phi_init,
      gam = gam_init,
      par = c(be_init, log(phi_init), gam_init)
    )
  } else {
    list(
      be  = be_init,
      phi = phi_init,
      gam = NULL,
      par = c(be_init, log(phi_init))
    )
  }
}


# Internal inference helper: obtain the observed Hessian from the analytic gradient,
# invert it for vcov, and construct standard errors / z tests / p-values.
.zimp_compute_inference <- function(par_hat, fn.gr, p, p1 = 0,
                                   zero_inflation = FALSE,
                                   eps_se = 1e-12) {
  
  if (!requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Package 'numDeriv' is required. Please install.packages('numDeriv').")
  }
  
  format_p_lm <- function(p) {
    ifelse(
      is.na(p), NA_character_,
      ifelse(p < 2e-16, "< 2e-16",
             format.pval(p, digits = 3, eps = 2e-16))
    )
  }
  
  signif_code_lm <- function(p) {
    ifelse(
      is.na(p), "",
      ifelse(p < 0.001, "***",
             ifelse(p < 0.01, "**",
                    ifelse(p < 0.05, "*", "")))
    )
  }
  
  make_coef_table <- function(est, se, par_names) {
    zval <- est / se
    pval <- 2 * stats::pnorm(-abs(zval))
    
    data.frame(
      Estimate = est,
      `Std. Error` = se,
      `z value` = zval,
      `Pr(>|z|)` = format_p_lm(pval),
      ` ` = signif_code_lm(pval),
      row.names = par_names,
      check.names = FALSE
    )
  }
  
  H <- tryCatch(
    numDeriv::jacobian(fn.gr, par_hat),
    error = function(e) NULL
  )
  
  if (is.null(H)) {
    return(list(
      hessian = NULL,
      vcov = NULL,
      coef_table = NULL,
      coef_table_internal = NULL,
      hessian_ok = FALSE
    ))
  }
  
  H <- 0.5 * (H + t(H))
  
  V <- tryCatch(
    solve(H),
    error = function(e) {
      tryCatch(MASS::ginv(H), error = function(e2) NULL)
    }
  )
  
  if (is.null(V)) {
    return(list(
      hessian = H,
      vcov = NULL,
      coef_table = NULL,
      coef_table_internal = NULL,
      hessian_ok = FALSE
    ))
  }
  
  se <- sqrt(pmax(diag(V), eps_se))
  
  par_names_internal <- c(
    paste0("beta_", seq_len(p)),
    "log_phi",
    if (zero_inflation) paste0("gamma_", seq_len(p1)) else NULL
  )
  
  coef_table_internal <- make_coef_table(
    est = par_hat,
    se = se,
    par_names = par_names_internal
  )
  
  phi_hat <- exp(par_hat[p + 1])
  se_phi <- phi_hat * se[p + 1]
  
  est_report <- c(
    par_hat[seq_len(p)],
    phi_hat,
    if (zero_inflation) par_hat[(p + 2):(p + 1 + p1)] else NULL
  )
  
  se_report <- c(
    se[seq_len(p)],
    se_phi,
    if (zero_inflation) se[(p + 2):(p + 1 + p1)] else NULL
  )
  
  par_names_report <- c(
    paste0("beta_", seq_len(p)),
    "phi",
    if (zero_inflation) paste0("gamma_", seq_len(p1)) else NULL
  )
  
  coef_table <- make_coef_table(
    est = est_report,
    se = se_report,
    par_names = par_names_report
  )
  
  list(
    hessian = H,
    vcov = V,
    coef_table = coef_table,
    coef_table_internal = coef_table_internal,
    hessian_ok = TRUE
  )
}
control = list(
  trace    = trace,
  iter.max = 200L,
  eval.max = 400L,
  rel.tol  = 1e-8,
  x.tol    = 1e-8
)


# =============================================================================
# LOW-LEVEL FITTER: FIT EXACTLY ONE MARGINAL MODEL
# =============================================================================
# Responsibility:
#   - Fit one fully specified PM/zimp marginal model.
#   - No model search and no comparison across distributions.
#   - This is the computational engine reused by fit_zimp() and select_zimp().
#
# Parameterization used by nlminb:
#   par = (beta, log(phi), gamma)
# where gamma is present only for a zero-inflated model.
# Internal numerical engine: fit exactly ONE already-prepared PM/zimp model.
# This function performs initialization, likelihood/gradient construction,
# nlminb optimization, and post-estimation inference. It does no model selection.
# `optimize = FALSE` is an internal-only mode used by select_zimp() to compute
# inference at an already optimized parameter vector without re-optimization.
.fit_zimp_raw <- function(y, X, Z = NULL, offset.count = NULL,
                          dist.family = c("gamma", "lnorm",
                                          "ig", "weibull", "loglogistic", "betaprime"),
                          zero_inflation = FALSE,
                          var.power = NULL,
                          ini.par = NULL,
                          progress = TRUE,
                          trace = 1,
                          eval.max = 200,
                          iter.max = 100,
                          rel.tol = 1e-8,
                          x.tol = 1e-8,
                          n_gh = 20,
                          gh = NULL,
                          compute.inference = TRUE,
                          optimize = TRUE) {

  dist.family <- match.arg(dist.family)  # exactly one distribution
  X <- as.matrix(X)
  y <- as.numeric(y)

  if (length(y) != nrow(X)) stop("length(y) must equal nrow(X).")
  if (!is.null(var.power) &&
      (length(var.power) != 1L || !is.numeric(var.power) || !is.finite(var.power) ||
       !(var.power %in% c(1, 2, 3)))) {
    stop("For fit_zimp(), var.power must be NULL, 1, 2, or 3.")
  }
  use_custom_var_power <- !is.null(var.power)
  var_power <- if (is.null(var.power)) 2 else as.integer(var.power)

  p <- ncol(X)
  if (zero_inflation) {
    Z <- as.matrix(Z)
    if (length(y) != nrow(Z)) stop("length(y) must equal nrow(Z).")
    p1 <- ncol(Z)
  } else {
    Z <- NULL
    p1 <- 0L
  }

  n_par <- p + 1L + p1
  if (!is.null(ini.par) && length(ini.par) != n_par) {
    stop("length(ini.par) is not correct for the chosen zero_inflation setting.")
  }

  .progress_msg <- function(...) {
    if (isTRUE(progress)) {
      cat(..., "\n", sep = "")
      flush.console()
    }
  }
 

 model_type <- if (zero_inflation) "Zero-Inflated Mixed Poisson" else "Mixed Poisson"
 .progress_msg(
  "\nStarting model fitting\n",
  "Model = ", model_type,
  " | mixing distribution = ", dist.family,
  if (!is.null(var.power)) paste0(" | var.power = ", var.power) else ""
)

  # Lognormal mixing requires Gauss-Hermite nodes/weights.  For other
  # distributions gh0 remains NULL because their marginal likelihood code
  # does not use this quadrature object.
  gh0 <- NULL
  if (dist.family == "lnorm") {
    gh0 <- if (is.null(gh)) statmod::gauss.quad(n_gh, kind = "hermite") else gh
  }

  # Objective passed to nlminb: negative marginal log-likelihood.
  ff <- function(par) {
    if (any(!is.finite(par))) return(1e100)
    be1  <- par[1:p]
    phi1 <- exp(par[p + 1L])
    gam1  <- if (zero_inflation) par[(p + 2L):(p + 1L + p1)] else NULL

    val <- tryCatch(
      neg_loglik_marginal(
        y = y, X = cbind(X, offset.count), Z = Z,
        dist = dist.family,
        be = c(be1, 1), phi = phi1, al = gam1,
        zero_inflation = zero_inflation,
        use_custom_var_power = use_custom_var_power,
        var_power = var_power,
        n_gh = n_gh, gh = gh0
      ),
      error = function(e) {
        message("Objective error [", dist.family, "]: ", e$message)
        1e100
      }
    )
    if (!is.finite(val)) 1e100 else val
  }

  # Analytic gradient of the negative log-likelihood.  Since optimization is
  # performed in log(phi), dL/dlog(phi) = dL/dphi * phi.
  fn.gr <- function(par) {
    if (any(!is.finite(par))) return(rep(1e6, n_par))
    be1  <- par[1:p]
    phi1 <- exp(par[p + 1L])
    gam1  <- if (zero_inflation) par[(p + 2L):(p + 1L + p1)] else NULL

    re1 <- tryCatch(
      neg_loglik_marginal_grad(
        y = y, X = cbind(X, offset.count), Z = Z,
        dist = dist.family,
        be = c(be1, 1), phi = phi1, al = gam1,
        zero_inflation = zero_inflation,
        use_custom_var_power = use_custom_var_power,
        var_power = var_power,
        n_gh = n_gh, gh = gh0
      ),
      error = function(e) {
        message("Gradient error [", dist.family, "]: ", e$message)
        NULL
      }
    )

    if (is.null(re1)) return(rep(1e6, n_par))
    g <- if (zero_inflation) {
      c(re1$d_be[seq_len(p)], re1$d_phi * phi1, re1$d_al)
    } else {
      c(re1$d_be[seq_len(p)], re1$d_phi * phi1)
    }
    g[!is.finite(g)] <- 1e6
    g
  }

  # If the user does not supply starting values, construct stable starts from
  # simple binomial/Poisson regressions and moment-based dispersion estimates.
  ini.par.0 <- if (is.null(ini.par)) {
    .zimp_init_par(
      y = y, X = X, Z = Z, offset.count = offset.count,
      dist = dist.family,
      zero_inflation = zero_inflation,
      var.power = var.power
    )$par
  } else {
    as.numeric(ini.par)
  }

  # Optimize once unless this internal call is used only to compute inference
  # at an already estimated parameter vector.
  if (isTRUE(optimize)) {
    re <- stats::nlminb(
      start = ini.par.0,
      objective = ff,
      gradient = fn.gr,
      control = list(trace = trace, eval.max = eval.max, iter.max = iter.max,
                     rel.tol = rel.tol, x.tol = x.tol)
    )
  } else {
    if (is.null(ini.par)) stop("'ini.par' is required when optimize = FALSE.")
    re <- list(
      par = ini.par.0,
      objective = ff(ini.par.0),
      convergence = NA_integer_,
      message = "Inference evaluated at supplied parameter vector; no optimization performed."
    )
  }

  # Numerical Hessian/Jacobian of the analytic score is used only after
  # optimization to obtain vcov, standard errors, z statistics and p-values.
  infer <- if (isTRUE(compute.inference)) {
    .zimp_compute_inference(
      par_hat = re$par, fn.gr = fn.gr, p = p, p1 = p1,
      zero_inflation = zero_inflation
    )
  } else {
    list(hessian = NULL, vcov = NULL, coef_table = NULL,
         coef_table_internal = NULL, hessian_ok = NA)
  }

  be1  <- re$par[1:p]
  phi1 <- exp(re$par[p + 1L])
  gam1  <- if (zero_inflation) re$par[(p + 2L):(p + 1L + p1)] else NULL

  list(
    distribution = dist.family,
    var.power = var.power,
    neg_log_lik_value = re$objective,
    log_lik_value = -re$objective,
    be = be1, phi = phi1, gam = gam1, par = re$par,
    convergence = re$convergence, message = re$message,
    hessian = infer$hessian, vcov = infer$vcov,
    coef_table = infer$coef_table,
    coef_table_internal = infer$coef_table_internal,
    hessian_ok = infer$hessian_ok
  )
}


# =============================================================================
# FORMULA/DATA PREPARATION SHARED BY fit_zimp() AND select_zimp()
# =============================================================================
# Responsibility:
#   - Parse y ~ X | Z using the package's explicit-intercept convention.
#   - Validate variables and response counts.
#   - Remove incomplete rows once, consistently for both count and zero layers.
#   - Build y, X and (when present) Z.
#
# Keeping this logic in one helper guarantees that fit_zimp() and
# select_zimp() interpret an identical formula in exactly the same way.
# Internal interface helper: parse `y ~ X | Z`, validate variables/count response,
# remove incomplete rows, and build y, X, Z plus formula/model-frame metadata.
.zimp_prepare_formula_data <- function(formula, data) {
  if (missing(formula)) stop("formula is required.")
  if (missing(data)) stop("data is required.")
  if (!inherits(formula, "formula")) stop("formula must be a formula.")
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (length(formula) != 3L) {
    stop("formula must have the form y ~ x1 + x2 or y ~ x1 + x2 | z1 + z2.")
  }

  lhs <- formula[[2L]]
  rhs <- formula[[3L]]
  has_zero_part <- is.call(rhs) && identical(rhs[[1L]], as.name("|"))

  if (has_zero_part) {
    rhs_count <- rhs[[2L]]
    rhs_zero  <- rhs[[3L]]
    zero_inflation <- TRUE
  } else {
    rhs_count <- rhs
    rhs_zero <- NULL
    zero_inflation <- FALSE
  }

  formula_env <- environment(formula)
  count_spec <- .zimp_make_rhs_formula(rhs_count, formula_env, part = "count")

  deparse1_safe <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = "")
  lhs_chr <- deparse1_safe(lhs)
  f_count <- stats::as.formula(
    paste(lhs_chr, "~", deparse1_safe(count_spec$formula[[2L]])),
    env = formula_env
  )

  if (zero_inflation) {
    zero_spec <- .zimp_make_rhs_formula(rhs_zero, formula_env, part = "zero")
    f_zero <- zero_spec$formula
    zero_has_intercept <- zero_spec$has_intercept
  } else {
    f_zero <- NULL
    zero_has_intercept <- FALSE
  }

  vars_needed <- all.vars(formula)
  vars_missing <- setdiff(vars_needed, names(data))
  if (length(vars_missing) > 0L) {
    stop("The following variables are not found in data: ", paste(vars_missing, collapse = ", "))
  }

  cc <- stats::complete.cases(data[, vars_needed, drop = FALSE])
  if (!all(cc)) warning(sum(!cc), " rows removed because of missing values.")
  data_use <- data[cc, , drop = FALSE]
  if (nrow(data_use) == 0L) stop("No complete observations remain after removing missing values.")

  mf_count <- stats::model.frame(f_count, data = data_use, na.action = stats::na.pass)
  y <- stats::model.response(mf_count)
  if (is.factor(y) || is.character(y)) stop("Response variable must be numeric/integer count data.")
  y <- as.numeric(y)
  if (any(!is.finite(y))) stop("Response variable contains non-finite values.")
  if (any(y < 0)) stop("Response variable must be non-negative.")
  if (any(abs(y - round(y)) > 1e-8)) stop("Response variable must contain integer counts.")
  y <- as.integer(round(y))

  offset.count <- stats::model.offset(mf_count)
  if (is.null(offset.count)) offset.count <- rep(0, length(y))
  offset.count <- as.numeric(offset.count)

  X <- stats::model.matrix(f_count, data = mf_count)
  if (ncol(X) == 0L) stop("Count model matrix X has zero columns.")

  if (zero_inflation) {
    mf_zero <- stats::model.frame(f_zero, data = data_use, na.action = stats::na.pass)
    Z <- stats::model.matrix(f_zero, data = mf_zero)
    if (ncol(Z) == 0L) stop("Zero-inflation formula produces a Z matrix with zero columns.")
  } else {
    mf_zero <- NULL
    Z <- NULL
  }

  list(
    y = y, X = X, Z = Z, offset.count = offset.count,
    zero_inflation = zero_inflation,
    count_formula = f_count, zero_formula = f_zero,
    count_intercept = count_spec$has_intercept,
    zero_intercept = if (zero_inflation) zero_has_intercept else NA,
    model_frame_count = mf_count, model_frame_zero = mf_zero,
    na_removed = which(!cc)
  )
}



# =============================================================================
# INTERNAL INPUT ROUTER
# =============================================================================
# Converts either supported user interface into the same validated internal
# representation (y, X, Z).
#
# Interface 1: formula + data
# Interface 2: Y + X + optional Z
#
# The two interfaces are mutually exclusive.
.zimp_resolve_input <- function(formula = NULL, data = NULL,
                                Y = NULL, X = NULL, Z = NULL, offset.count = NULL) {

  formula_used <- !is.null(formula) || !is.null(data)
  direct_used  <- !is.null(Y) || !is.null(X) || !is.null(Z) || !is.null(offset.count)

  if (formula_used && direct_used) {
    stop("Use either 'formula' + 'data' OR 'Y' + 'X' (+ optional 'Z'), not both.")
  }
  if (!formula_used && !direct_used) {
    stop("Supply either 'formula' + 'data' or 'Y' + 'X' (+ optional 'Z').")
  }

  # ---- Formula interface -----------------------------------------------------
  if (formula_used) {
    if (!is.null(offset.count)) stop("With the formula interface, specify the count offset using offset() in the formula, not via 'offset.count'.")
    if (is.null(formula) || is.null(data)) {
      stop("The formula interface requires both 'formula' and 'data'.")
    }
    if (!inherits(formula, "formula")) stop("'formula' must be a formula.")
    if (!is.data.frame(data)) stop("'data' must be a data frame.")

    out <- .zimp_prepare_formula_data(formula, data)
    out$interface <- "formula"
    return(out)
  }

  # ---- Direct Y/X/Z interface -----------------------------------------------
  if (is.null(Y) || is.null(X)) {
    stop("The direct-data interface requires both 'Y' and 'X'.")
  }

  # Y must really be one-dimensional; n x 1 matrices/data frames are rejected.
  if (is.matrix(Y) || is.data.frame(Y) || !is.null(dim(Y))) {
    stop("'Y' must be a one-dimensional numeric vector, not a matrix/data frame.")
  }
  if (!is.numeric(Y) || !is.atomic(Y)) stop("'Y' must be a numeric vector.")
  if (length(Y) < 1L) stop("'Y' must contain at least one observation.")
  if (any(!is.finite(Y))) stop("'Y' contains NA/NaN/Inf values.")
  if (any(Y < 0)) stop("'Y' must contain non-negative counts.")
  if (any(abs(Y - round(Y)) > 1e-8)) stop("'Y' must contain integer counts.")

  y <- as.integer(round(Y))
  n <- length(y)

  if (is.null(offset.count)) offset.count <- rep(0, n)
  if (!is.numeric(offset.count) || is.matrix(offset.count) || is.data.frame(offset.count) || !is.null(dim(offset.count)))
    stop("'offset.count' must be NULL or a one-dimensional numeric vector.")
  offset.count <- as.numeric(offset.count)
  if (length(offset.count) != n) stop("'offset.count' must have length(Y) elements.")
  if (any(!is.finite(offset.count))) stop("'offset.count' contains NA/NaN/Inf values.")

  validate_design <- function(M, nm) {
    if (!(is.matrix(M) || is.data.frame(M))) {
      stop("'", nm, "' must be a matrix or data frame.")
    }
    if (length(dim(M)) != 2L) stop("'", nm, "' must be two-dimensional.")
    if (nrow(M) != n) {
      stop("'", nm, "' has ", nrow(M), " rows but length(Y) = ", n, ".")
    }
    if (ncol(M) < 1L) stop("'", nm, "' must contain at least one column.")

    if (is.data.frame(M)) {
      ok <- vapply(M, is.numeric, logical(1))
      if (!all(ok)) {
        bad <- names(M)[!ok]
        stop("'", nm, "' must contain only numeric columns. Non-numeric: ",
             paste(bad, collapse = ", "), ".")
      }
    } else if (!is.numeric(M)) {
      stop("'", nm, "' must be a numeric matrix.")
    }

    M <- as.matrix(M)
    storage.mode(M) <- "double"
    if (any(!is.finite(M))) stop("'", nm, "' contains NA/NaN/Inf values.")
    M
  }

  # Preserve user-supplied column names before coercion. If no names were
  # supplied, coefficients will be displayed as beta1, beta2, ... / gamma1, ...
  X_names <- colnames(X)
  X_has_names <- !is.null(X_names) && length(X_names) == ncol(X) &&
    all(!is.na(X_names)) && all(nzchar(X_names))

  Z_names <- if (!is.null(Z)) colnames(Z) else NULL
  Z_has_names <- !is.null(Z) && !is.null(Z_names) &&
    length(Z_names) == ncol(Z) && all(!is.na(Z_names)) && all(nzchar(Z_names))

  X <- validate_design(X, "X")
  if (X_has_names) colnames(X) <- X_names

  if (is.null(Z)) {
    Z <- NULL
    zero_inflation <- FALSE
  } else {
    Z <- validate_design(Z, "Z")
    if (Z_has_names) colnames(Z) <- Z_names
    zero_inflation <- TRUE
  }

  list(
    y = y, X = X, Z = Z, offset.count = offset.count,
    zero_inflation = zero_inflation,
    interface = "direct",
    formula = NULL,
    count_formula = NULL,
    zero_formula = NULL,
    count_intercept = NA,
    zero_intercept = NA,
    model_frame_count = NULL,
    model_frame_zero = NULL,
    na_removed = integer(0),
    X_has_names = X_has_names,
    Z_has_names = Z_has_names
  )
}

# =============================================================================
# INTERNAL OUTPUT HELPERS
# =============================================================================

# Choose readable coefficient labels.
# Formula input uses model.matrix() column names, including "(Intercept)" and
# expanded factor/one-hot columns. Direct matrix/data-frame input uses supplied
# column names when available; otherwise beta1,... and gamma1,... are used.
.zimp_coef_names <- function(d) {
  p <- ncol(d$X)
  q <- if (d$zero_inflation) ncol(d$Z) else 0L

  if (identical(d$interface, "formula")) {
    beta_names <- colnames(d$X)
    gamma_names <- if (q > 0L) colnames(d$Z) else character(0)
  } else {
    beta_names <- if (isTRUE(d$X_has_names)) colnames(d$X) else paste0("beta", seq_len(p))
    gamma_names <- if (q > 0L) {
      if (isTRUE(d$Z_has_names)) colnames(d$Z) else paste0("gamma", seq_len(q))
    } else character(0)
  }

  if (is.null(beta_names) || length(beta_names) != p)
    beta_names <- paste0("beta", seq_len(p))
  if (q > 0L && (is.null(gamma_names) || length(gamma_names) != q))
    gamma_names <- paste0("gamma", seq_len(q))

  list(beta = beta_names, gamma = gamma_names)
}

# Convert the engine's combined coefficient table into two user-facing tables.
.zimp_split_coef_table <- function(fit, d) {
  nm <- .zimp_coef_names(d)
  p <- ncol(d$X)
  q <- if (d$zero_inflation) ncol(d$Z) else 0L
  tab <- fit$coef_table

  if (is.null(tab) || nrow(tab) < p + 1L + q) {
    return(list(count = NULL, zero = NULL, names = nm))
  }

  count <- tab[seq_len(p), , drop = FALSE]
  rownames(count) <- nm$beta

  zero <- NULL
  if (q > 0L) {
    zero_idx <- p + 1L + seq_len(q)
    zero <- tab[zero_idx, , drop = FALSE]
    rownames(zero) <- nm$gamma
  }

  list(count = count, zero = zero, names = nm)
}

# Pretty console printer used at the end of fit_zimp().
.zimp_print_fit <- function(fit, show.call = TRUE) {

  # Print coefficient tables in an lm-like layout. The last column deliberately
  # has a blank header and contains only significance stars.
  print_coef_block <- function(tab) {
    if (is.null(tab)) {
      cat("Inference unavailable.\n")
      return(invisible(NULL))
    }
    print(tab, quote = FALSE)
    invisible(NULL)
  }

  if (isTRUE(show.call)) {
   cat("\nCall:\n")
   print(fit$call)
  }

  cat("\nCount model coefficients:\n")
  print_coef_block(fit$count_coef_table)

  if (isTRUE(fit$zero_inflation)) {
    cat("\nZero-inflation model coefficients:\n")
    print_coef_block(fit$zero_coef_table)
  }

  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 ' ' 1\n")

  # phi is reported on its natural scale. Its standard error is obtained by
  # the delta method from log(phi), already computed in the inference engine.
  phi_se <- fit$phi_se
  if (is.null(phi_se) || length(phi_se) != 1L || !is.finite(phi_se)) {
    phi_se_txt <- "NA"
  } else {
    phi_se_txt <- formatC(phi_se, digits = 5, format = "fg")
  }

  cat("\nDispersion parameter:\n")
  cat("phi = ", formatC(fit$phi, digits = 6, format = "fg"),
      "  (Std. Error = ", phi_se_txt, ")\n", sep = "")

  cat("\nLog-likelihood: ",
      formatC(fit$log_lik_value, digits = 6, format = "f"), sep = "")
  cat("\nAIC: ",
      formatC(-2*fit$log_lik_value + 2*sum(c(length(fit$be),length(fit$gam),1),na.rm=TRUE), digits = 2, format = "f"), sep = "")
  cat("\nBIC: ",
      formatC(-2*fit$log_lik_value + log(length(fit$y))*sum(c(length(fit$be),length(fit$gam),1),na.rm=TRUE), digits = 2, format = "f"), "\n", sep = "")

  cat("\nConvergence: ", ifelse(fit$convergence == 0, "Yes", "No") , "\n", sep = "")

  invisible(fit)
}

# =============================================================================
# PUBLIC API 2: select_zimp() -- MULTIPLE CANDIDATE MODELS
# =============================================================================
#' Select a Univariate Mixed Poisson or Zero-Inflated Mixed Poisson Model
#'
#' @description
#' Fits and compares multiple univariate Mixed Poisson (MP) or Zero-Inflated
#' Mixed Poisson (ZIMP) regression models across candidate mixing
#' distributions and variance-power specifications. Candidate models are
#' ranked by maximized log-likelihood, and the model with the largest
#' log-likelihood is selected.
#'
#' @details
#' `select_zimp()` is the model-selection counterpart of [fit_zimp()].
#' While [fit_zimp()] fits one fully specified model, `select_zimp()` fits a
#' collection of candidate models and returns the best-fitting specification.
#'
#' The function supports the same formula and matrix interfaces as
#' [fit_zimp()].
#'
#' @section Model specification:
#'
#' **Formula interface**
#'
#' A Mixed Poisson model can be specified as
#'
#' ```
#' y ~ x1 + x2
#' ```
#'
#' and a Zero-Inflated Mixed Poisson model as
#'
#' ```
#' y ~ x1 + x2 | z1 + z2
#' ```
#'
#' The expression before `|` specifies the count component and the expression
#' after `|` specifies the zero-inflation component.
#'
#' In the `mzimp` formula interface, intercepts are included only when
#' explicitly requested with `1`. Thus,
#'
#' ```
#' y ~ x1 + x2
#' ```
#'
#' specifies a model without an intercept, whereas
#'
#' ```
#' y ~ 1 + x1 + x2
#' ```
#'
#' includes an intercept. The same convention applies independently to the
#' zero-inflation component.
#'
#' **Matrix interface**
#'
#' The response and design matrices may instead be supplied through `Y`, `X`,
#' and optionally `Z`. If `Z = NULL`, MP candidate models are fitted. If `Z`
#' is supplied, ZIMP candidate models are fitted. Intercept columns, when
#' required, must be included explicitly in `X` and `Z`.
#'
#' @section Candidate models:
#'
#' Candidate models are generated from the requested `mixing.distribution`
#' and `var.power` values.
#'
#' The supported mixing distributions are `"gamma"`, `"lnorm"`, `"ig"`,
#' `"weibull"`, `"loglogistic"`, and `"betaprime"`.
#'
#' If `var.power = NULL`, one model is fitted for each requested mixing
#' distribution using its native/default variance-power specification:
#'
#' \describe{
#'   \item{`"gamma"`}{`var.power = 2`}
#'   \item{`"lnorm"`}{`var.power = 2`}
#'   \item{`"ig"`}{`var.power = 3`}
#'   \item{`"weibull"`}{`var.power = 2`}
#'   \item{`"loglogistic"`}{`var.power = 2`}
#'   \item{`"betaprime"`}{`var.power = 2`}
#' }
#'
#' If one or more values are supplied to `var.power`, every requested mixing
#' distribution is fitted for each requested variance power. For example,
#'
#' ```
#' mixing.distribution = c("gamma", "lnorm")
#' var.power = c(1, 2, 3)
#' ```
#'
#' defines six candidate models.
#'
#' Duplicate distribution and variance-power specifications are removed
#' before fitting.
#'
#' @section Model fitting and selection:
#'
#' Candidate models are fitted by maximum likelihood using the same
#' distribution-specific likelihood and gradient routines as [fit_zimp()].
#'
#' Gamma, lognormal, and inverse-Gaussian candidates are fitted directly.
#' Weibull, log-logistic, and beta-prime candidates use matched Gamma models
#' as initialization anchors. For these candidates, the count coefficients
#' and, for ZIMP models, the zero-component coefficients are initialized from
#' a Gamma model with the same effective variance power. The dispersion
#' parameter is initialized using the target distribution's own
#' initialization procedure.
#'
#' When a required Gamma anchor is not among the requested candidate models,
#' it is fitted internally for initialization only and is not included in the
#' returned candidate table.
#'
#' Candidate models are initially fitted without post-estimation inference.
#' They are then ordered by decreasing maximized log-likelihood. The first row
#' of `summary_table` therefore corresponds to the selected model.
#'
#' Full post-estimation inference is subsequently computed only for the
#' selected model, at its already optimized parameter estimates; the selected
#' model is not re-optimized. The resulting `best_model` is an object of class
#' `"zimp_fit"`.
#'
#' Selection is based exclusively on maximized log-likelihood. AIC and BIC
#' are available for the selected model through `best_model`, but are not
#' used to rank candidate models.
#'
#' @section Parallel computation:
#'
#' Candidate models within the same fitting stage may be evaluated in
#' parallel by setting `parallel` to an integer greater than one.
#'
#' Parallelization is performed across candidate models; optimization of an
#' individual candidate is not itself parallelized by `select_zimp()`.
#' Dependencies between fitting stages, including construction of Gamma
#' initialization anchors, are preserved.
#'
#' On Unix-like systems, fork-based parallel processing may be used. On
#' Windows, PSOCK workers are used.
#'
#' @inheritParams fit_zimp
#'
#' @param mixing.distribution Character vector specifying one or more
#'   candidate mixing distributions. Supported values are `"gamma"`,
#'   `"lnorm"`, `"ig"`, `"weibull"`, `"loglogistic"`, and `"betaprime"`.
#'   By default, all supported distributions are considered.
#'
#' @param var.power Optional numeric vector specifying one or more
#'   variance-power parameterizations. Values must be selected from `1`, `2`,
#'   and `3`. If `NULL`, each requested mixing distribution is fitted using
#'   its native/default variance-power specification.
#'
#' @param parallel `NULL` or a positive integer specifying the number of
#'   workers requested for parallel candidate-model fitting. `NULL` or `1`
#'   uses serial execution. The number of workers actually used in a fitting
#'   stage cannot exceed the number of candidate jobs in that stage.
#'
#' @param show.result Logical compatibility argument. Selection is performed
#'   silently and results are stored in the returned object. The selection
#'   report can be displayed by evaluating or printing the `summary`
#'   component.
#'
#' @param .return.all Logical. Internal-use option. If `FALSE`, the default,
#'   only the compact public selection object is returned. If `TRUE`,
#'   candidate fit objects are additionally retained in `all_results`.
#'   Users should normally leave this argument at its default value.
#'
#' @return
#' An object of class `"zimp_selection"`. By default, the object contains:
#'
#' \describe{
#'   \item{`summary_table`}{
#'     A data frame containing one row for each requested candidate model,
#'     ordered by decreasing maximized log-likelihood. Columns include
#'     `distribution`, `var.power`, `log_lik`, `neg_log_lik`, and `phi`.
#'   }
#'   \item{`best_model`}{
#'     The selected model, returned as an object of class `"zimp_fit"`.
#'     It corresponds to the first row of `summary_table` and contains full
#'     post-estimation inference.
#'   }
#'   \item{`summary`}{
#'     An object of class `"zimp_selection_summary"` containing the candidate
#'     table and selected model for formatted display.
#'   }
#'   \item{`call`}{
#'     The matched call to `select_zimp()`.
#'   }
#' }
#'
#' If `.return.all = TRUE`, an additional internal component `all_results`
#' containing the candidate fit objects is returned.
#'
#' @seealso
#' [fit_zimp()] for fitting one specified univariate MP or ZIMP model.
#'
#' @examples
#' ## Generate a small Mixed Poisson data set
#' set.seed(12345)
#' n <- 50
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' X <- cbind(1, x1, x2)
#' mu <- exp(X %*% c(0.5, 0.3, 0.7))
#'
#' ## Gamma-mixed Poisson data
#' phi <- 0.2
#' y <- rnbinom(n, size = 1 / phi, mu = mu)
#'
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' ## Compare candidate mixing distributions
#' sel <- select_zimp(
#'   y ~ 1 + x1 + x2,
#'   data = dat,
#'   mixing.distribution = c("gamma", "lnorm", "ig"),
#'   progress = FALSE
#' )
#'
#' ## Ranked candidate models
#' sel$summary_table
#'
#' ## Selected model
#' sel$best_model
#'
#' ## Compare variance-power specifications
#' sel_power <- select_zimp(
#'   Y = y,
#'   X = X,
#'   mixing.distribution = c("gamma", "lnorm"),
#'   var.power = c(1, 2, 3),
#'   progress = FALSE
#' )
#'
#' sel_power$summary_table
#'
#' @export

select_zimp <- function(formula = NULL, data = NULL,
                         Y = NULL, X = NULL, Z = NULL, offset.count = NULL,
                         mixing.distribution = c("gamma","lnorm","ig",
                                                 "weibull","loglogistic","betaprime"),
                         var.power = NULL, n_gh = 20, gh = NULL, ini.par = NULL,
                         parallel = NULL, progress = TRUE, show.result = TRUE, trace = 0,
                         eval.max = 100, rel.tol = 1e-5,
                         x.tol = 1e-5, iter.max = 100,
                         .return.all = FALSE) {

  d <- .zimp_resolve_input(formula, data, Y, X, Z, offset.count)

  if (length(.return.all) != 1L || !is.logical(.return.all) || is.na(.return.all)) stop("'.return.all' must be TRUE or FALSE.")

  if (length(show.result) != 1L || !is.logical(show.result) || is.na(show.result)) {
    stop("'show.result' must be TRUE or FALSE.")
  }

  choices <- c("gamma","lnorm","ig",
               "weibull","loglogistic","betaprime")
  mixing.distribution <- unique(
    match.arg(mixing.distribution, choices, several.ok = TRUE)
  )

  if (!is.null(parallel)) {
    if (length(parallel) != 1L || !is.numeric(parallel) ||
        !is.finite(parallel) || parallel < 1 ||
        parallel != as.integer(parallel)) {
      stop("'parallel' must be NULL or a positive integer.")
    }
    parallel <- as.integer(parallel)
  }

  default_p <- c(
    gamma = 2L, lnorm = 2L, pareto = 2L, ig = 3L,
    weibull = 2L, loglogistic = 2L, betaprime = 2L
  )

  if (!is.null(var.power)) {
    if (!is.numeric(var.power) || any(!is.finite(var.power)) ||
        any(!var.power %in% 1:3)) {
      stop("'var.power' must be NULL or contain only 1, 2, or 3.")
    }
    var.power <- unique(as.integer(var.power))
  }

  has_pareto <- "pareto" %in% mixing.distribution
  other_dist <- setdiff(mixing.distribution, "pareto")

  if (is.null(var.power)) {
    grid <- data.frame(
      distribution = mixing.distribution,
      fit.var.power = rep(NA_integer_, length(mixing.distribution)),
      var.power = unname(default_p[mixing.distribution]),
      stringsAsFactors = FALSE
    )
  } else {
    parts <- list()
    if (length(other_dist)) {
      g <- expand.grid(
        distribution = other_dist,
        fit.var.power = var.power,
        stringsAsFactors = FALSE
      )
      g$var.power <- g$fit.var.power
      parts[[length(parts) + 1L]] <- g
    }
    if (has_pareto) {
      parts[[length(parts) + 1L]] <- data.frame(
        distribution = "pareto",
        fit.var.power = NA_integer_,
        var.power = 2L,
        stringsAsFactors = FALSE
      )
    }
    grid <- do.call(rbind, parts)
    rownames(grid) <- NULL
  }

  n_candidates <- nrow(grid)

  # select_zimp() is intentionally silent. Candidate fits also run silently.

  fits <- vector("list", n_candidates)

  # Buffer serial progress so it can be printed in candidate-grid order.
  serial_progress <- vector("list", n_candidates)

  # ---------------------------------------------------------------------------
  # Console helpers
  # ---------------------------------------------------------------------------
  vp_text <- function(p) if (is.null(p)) "default" else as.character(p)

  print_start <- function(...) invisible(NULL)

  print_done <- function(...) invisible(NULL)

  # ---------------------------------------------------------------------------
  # Self-contained job descriptor.
  # This is deliberately data-only so it can be sent safely to PSOCK workers.
  # ---------------------------------------------------------------------------
  make_job <- function(id, dist, p, start, kind, anchor_label = NULL, display_id = NULL) {
    list(
      id = id,
      display_id = display_id,
      dist = dist,
      p = p,
      start = start,
      kind = kind,
      anchor_label = anchor_label
    )
  }

  # ---------------------------------------------------------------------------
  # One job executor.
  #
  # IMPORTANT for Windows PSOCK:
  # workers receive the numerical functions explicitly below, so they do not
  # depend on closures such as raw_fit()/fit_idx() defined only in the master.
  # ---------------------------------------------------------------------------
  worker_fit <- function(job, dat, ctrl) {
    t0 <- proc.time()[["elapsed"]]

    fit <- .fit_zimp_raw(
      y = dat$y,
      X = dat$X,
      Z = dat$Z,
      offset.count = dat$offset.count,
      dist.family = job$dist,
      zero_inflation = dat$zero_inflation,
      var.power = job$p,
      ini.par = job$start,
      progress = FALSE,
      trace = ctrl$trace,
      eval.max = ctrl$eval.max,
      iter.max = ctrl$iter.max,
      rel.tol = ctrl$rel.tol,
      x.tol = ctrl$x.tol,
      n_gh = ctrl$n_gh,
      gh = ctrl$gh,
      compute.inference = FALSE
    )

    if (!is.finite(fit$log_lik_value) ||
        fit$neg_log_lik_value >= 1e90) {
      stop(
        "Invalid likelihood on worker for ", job$dist,
        " (var.power = ", if (is.null(job$p)) "default" else job$p, "). ",
        "The objective returned its numerical-error sentinel (>= 1e90). ",
        "Check that the package DLL/Rcpp likelihood code is loaded on the worker."
      )
    }

    list(
      job = job,
      fit = fit,
      elapsed = proc.time()[["elapsed"]] - t0
    )
  }

  ctrl <- list(
    trace = trace, eval.max = eval.max, iter.max = iter.max,
    rel.tol = rel.tol, x.tol = x.tol, n_gh = n_gh, gh = gh
  )

  dat_worker <- list(
    y = d$y, X = d$X, Z = d$Z,
    offset.count = d$offset.count,
    zero_inflation = d$zero_inflation
  )

  # ---------------------------------------------------------------------------
  # Find all functions needed recursively by .fit_zimp_raw().
  # We use codetools to avoid the old PSOCK error:
  # "could not find function '.fit_zimp_raw'".
  # ---------------------------------------------------------------------------
  collect_dependencies <- function(root_names, env = .GlobalEnv) {
    if (!requireNamespace("codetools", quietly = TRUE)) {
      stop("Package 'codetools' is required for Windows parallel execution.")
    }

    found <- character(0)
    queue <- unique(root_names)

    while (length(queue)) {
      nm <- queue[1L]
      queue <- queue[-1L]
      if (nm %in% found) next
      if (!exists(nm, envir = env, mode = "function", inherits = TRUE)) next

      fn <- get(nm, envir = env, mode = "function", inherits = TRUE)
      found <- c(found, nm)

      gl <- tryCatch(
        codetools::findGlobals(fn, merge = FALSE)$functions,
        error = function(e) character(0)
      )

      gl <- setdiff(gl, c(
        "{","(","[","[[","$","<-","=","::",":::",
        "+","-","*","/","^","%%","%/%","%*%",
        "<","<=",">",">=","==","!=","!","&","&&","|","||",
        ":","~","if","for","while","repeat","return","function"
      ))

      local_gl <- gl[vapply(
        gl,
        function(z) exists(z, envir = env, mode = "function", inherits = TRUE),
        logical(1)
      )]

      queue <- unique(c(queue, setdiff(local_gl, found)))
    }

    unique(found)
  }

  # ---------------------------------------------------------------------------
  # Run a batch either serially or in parallel.
  # Master prints START messages before dispatch and DONE messages as results
  # are returned. In PSOCK, worker stdout is also connected to the console.
  # ---------------------------------------------------------------------------
  run_jobs <- function(jobs, stage_label) {
    if (!length(jobs)) return(list())

    # Serial
    if (is.null(parallel) || parallel == 1L || length(jobs) == 1L) {
      out <- vector("list", length(jobs))
      for (j in seq_along(jobs)) {
        out[[j]] <- worker_fit(jobs[[j]], dat_worker, ctrl)
        did <- jobs[[j]]$display_id
        if (!is.null(did) && length(did) == 1L &&
            is.finite(did) && did >= 1L && did <= n_candidates) {
          serial_progress[[as.integer(did)]] <<- list(
            dist = jobs[[j]]$dist,
            p = jobs[[j]]$p
          )
        }
      }
      return(out)
    }

    nc <- min(parallel, length(jobs))

    # Unix/macOS fork
    if (.Platform$OS.type != "windows") {
      out <- parallel::mclapply(
        jobs,
        worker_fit,
        dat = dat_worker,
        ctrl = ctrl,
        mc.cores = nc,
        mc.preschedule = FALSE
      )
      if (FALSE) {
        for (j in seq_along(out)) {
        }
      }
      return(out)
    }

    # Windows PSOCK
    #
    # IMPORTANT:
    # Do NOT export only the R wrappers of the likelihood to workers. The
    # likelihood uses compiled Rcpp/C++ routines; exporting an R function does
    # not load its DLL/native symbols in a fresh PSOCK R session. In that case
    # .fit_zimp_raw() catches the objective error and returns 1e100, which looks
    # like a "converged" fit with logLik = -1e100.
    #
    # Therefore every worker loads the package itself (installed package first;
    # development tree through pkgload::load_all() otherwise). This loads both
    # the R namespace and the compiled shared library.

    cl <- parallel::makeCluster(nc, outfile = "")
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # Locate the development package root from the master session.
    find_pkg_root <- function(path = getwd()) {
      path <- normalizePath(path, winslash = "/", mustWork = TRUE)
      repeat {
        if (file.exists(file.path(path, "DESCRIPTION"))) return(path)
        parent <- dirname(path)
        if (identical(parent, path)) break
        path <- parent
      }
      NULL
    }

    pkg_root <- find_pkg_root()
    pkg_name <- NULL
    if (!is.null(pkg_root)) {
      desc <- tryCatch(
        read.dcf(file.path(pkg_root, "DESCRIPTION")),
        error = function(e) NULL
      )
      if (!is.null(desc) && "Package" %in% colnames(desc))
        pkg_name <- unname(desc[1L, "Package"])
    }

    # Initialize each PSOCK worker. Prefer the development tree because this
    # guarantees workers use the same current source/C++ code as the master.
    init_worker <- function(pkg_root, pkg_name) {
      ok <- FALSE
      how <- NULL

      if (!is.null(pkg_root) &&
          requireNamespace("pkgload", quietly = TRUE)) {
        tryCatch({
          pkgload::load_all(
            pkg_root,
            quiet = TRUE,
            recompile = FALSE,
            export_all = TRUE
          )
          ok <- TRUE
          how <- paste0("pkgload::load_all(", pkg_root, ")")
        }, error = function(e) {
          how <<- paste0("load_all failed: ", conditionMessage(e))
        })
      }

      if (!ok && !is.null(pkg_name) &&
          requireNamespace(pkg_name, quietly = TRUE)) {
        suppressPackageStartupMessages(
          library(pkg_name, character.only = TRUE)
        )
        ok <- TRUE
        how <- paste0("library(", pkg_name, ")")
      }

      list(ok = ok, how = how)
    }

    init_status <- parallel::clusterCall(
      cl, init_worker, pkg_root = pkg_root, pkg_name = pkg_name
    )

    bad <- which(!vapply(init_status, `[[`, logical(1), "ok"))
    if (length(bad)) {
      details <- vapply(
        init_status[bad],
        function(z) if (is.null(z$how)) "package could not be loaded" else z$how,
        character(1)
      )
      stop(
        "Parallel workers could not load the package/C++ likelihood code. ",
        "Run select_zimp() from the package project root and make sure pkgload ",
        "is installed. Worker message(s): ",
        paste(unique(details), collapse = " | ")
      )
    }

    # worker_fit is local to select_zimp(), so export only this small wrapper.
    # All numerical/package functions are resolved from the namespace loaded
    # above, including the compiled Rcpp routines.
    parallel::clusterExport(
      cl,
      varlist = "worker_fit",
      envir = environment()
    )

    # Fail early if the core fitter is not visible on workers. This avoids
    # silently converting a setup error into logLik = -1e100.
    fitter_ok <- parallel::clusterCall(
      cl,
      function() exists(".fit_zimp_raw", mode = "function", inherits = TRUE)
    )
    if (!all(unlist(fitter_ok))) {
      stop("'.fit_zimp_raw' is not available on one or more parallel workers.")
    }

    out <- parallel::parLapplyLB(
      cl,
      jobs,
      function(job, dat, ctrl) worker_fit(job, dat, ctrl),
      dat = dat_worker,
      ctrl = ctrl
    )

    if (FALSE) {
      for (j in seq_along(out)) {
      }
    }

    out
  }

  # ---------------------------------------------------------------------------
  # Target-specific starting vector:
  # beta and gamma from anchor; phi from target distribution initialization.
  # ---------------------------------------------------------------------------
  make_warm_start <- function(i, anchor) {
    p_fit <- if (is.na(grid$fit.var.power[i])) NULL else grid$fit.var.power[i]

    b <- .zimp_init_par(
      y = d$y, X = d$X, Z = d$Z, offset.count = d$offset.count,
      dist = grid$distribution[i],
      zero_inflation = d$zero_inflation,
      var.power = p_fit
    )

    b$be <- anchor$be
    if (isTRUE(d$zero_inflation)) b$gam <- anchor$gam

    c(
      b$be,
      log(b$phi),
      if (isTRUE(d$zero_inflation)) b$gam else NULL
    )
  }

  if (FALSE) {
    cat("\nModel selection: ", n_candidates, " candidates | ",
        if (is.null(parallel) || parallel == 1L) "serial"
        else paste0("parallel (", parallel, " cores)"),
        "\n", sep = "")
  }

  custom <- !is.null(var.power)
  anchors <- list()

  # ===========================================================================
  # STAGE 1: direct fits + matched Gamma anchors
  #
  # Warm-start rule:
  #   Weibull / Loglogistic / Betaprime with p=1 <- Gamma p=1
  #   Weibull / Loglogistic / Betaprime with p=2 <- Gamma p=2
  #   Weibull / Loglogistic / Betaprime with p=3 <- Gamma p=3
  #
  # beta and gamma (zero layer, when present) come from the matched Gamma fit.
  # phi remains initialized from the target distribution's own initializer.
  # Gamma anchors are internal unless Gamma itself is a requested candidate.
  # ===========================================================================

  warm_dist <- c("weibull", "loglogistic", "betaprime")

  warm_idx <- which(grid$distribution %in% warm_dist)
  gamma_anchor_p <- if (length(warm_idx)) {
    sort(unique(as.integer(grid$var.power[warm_idx])))
  } else {
    integer(0)
  }

  anchor_jobs <- list()

  # Requested Gamma / Lognormal / IG candidates are fitted independently.
  # Requested Gamma fits can simultaneously serve as matched anchors.
  direct_dist <- intersect(c("gamma", "lnorm", "ig"), mixing.distribution)

  if (length(direct_dist)) {
    for (dd in direct_dist) {
      idx_dd <- which(grid$distribution == dd)

      for (ii in idx_dd) {
        pp_eff <- as.integer(grid$var.power[ii])
        pp_fit <- if (is.na(grid$fit.var.power[ii])) NULL else
          as.integer(grid$fit.var.power[ii])

        anchor_jobs[[length(anchor_jobs) + 1L]] <- make_job(
          id = paste0("direct_", dd, "_p", pp_eff),
          dist = dd,
          p = pp_fit,
          start = if (identical(dd, "gamma")) ini.par else NULL,
          kind = "direct",
          display_id = ii
        )
      }
    }
  }

  # Determine which Gamma variance powers are already being fitted above.
  gamma_p_already <- integer(0)
  if (length(anchor_jobs)) {
    gamma_p_already <- vapply(
      anchor_jobs,
      function(jb) {
        if (!identical(jb$dist, "gamma")) return(NA_integer_)
        if (is.null(jb$p)) return(unname(default_p["gamma"]))
        as.integer(jb$p)
      },
      integer(1)
    )
    gamma_p_already <- unique(gamma_p_already[!is.na(gamma_p_already)])
  }

  # Add missing Gamma fits solely as internal anchors.
  gamma_p_missing <- setdiff(gamma_anchor_p, gamma_p_already)

  if (length(gamma_p_missing)) {
    for (pp in gamma_p_missing) {
      # Gamma's native/default variance power is 2.
      pp_fit <- if (pp == unname(default_p["gamma"])) NULL else pp

      anchor_jobs[[length(anchor_jobs) + 1L]] <- make_job(
        id = paste0("internal_gamma_anchor_p", pp),
        dist = "gamma",
        p = pp_fit,
        start = ini.par,
        kind = "anchor",
        anchor_label = paste0("gamma, var.power = ", pp)
      )
    }
  }

  anchor_results <- run_jobs(anchor_jobs, "[Anchor]")

  # Store Gamma fits by effective variance power and reuse requested direct fits.
  for (j in seq_along(anchor_results)) {
    aj <- anchor_jobs[[j]]

    pp <- if (is.null(aj$p)) {
      unname(default_p[aj$dist])
    } else {
      as.integer(aj$p)
    }

    if (identical(aj$dist, "gamma")) {
      anchors[[as.character(pp)]] <- anchor_results[[j]]$fit
    }

    if (!is.null(aj$display_id) &&
        length(aj$display_id) == 1L &&
        is.finite(aj$display_id)) {
      fits[[as.integer(aj$display_id)]] <- anchor_results[[j]]$fit
    }
  }

  # ===========================================================================
  # STAGE 2: remaining candidates, using matched anchor beta + gamma
  # ===========================================================================
  rem <- which(vapply(fits, is.null, logical(1)))

  if (length(rem)) {
    rem_jobs <- lapply(seq_along(rem), function(j) {
      i <- rem[j]
      pp <- grid$var.power[i]
      anchor <- anchors[[as.character(pp)]]

      if (grid$distribution[i] %in% warm_dist && is.null(anchor)) {
        stop(
          "Internal error: matched Gamma anchor for var.power = ",
          pp, " was not fitted."
        )
      }

      anchor_name <- if (grid$distribution[i] %in% warm_dist) {
        paste0("gamma, var.power = ", pp)
      } else {
        NULL
      }

      make_job(
        id = i,
        dist = grid$distribution[i],
        p = if (is.na(grid$fit.var.power[i])) NULL else grid$fit.var.power[i],
        start = if (grid$distribution[i] %in% warm_dist) {
          make_warm_start(i, anchor)
        } else {
          NULL
        },
        kind = "candidate",
        anchor_label = anchor_name
      )
    })

    for (jj in seq_along(rem_jobs)) rem_jobs[[jj]]$display_id <- rem_jobs[[jj]]$id
    rem_results <- run_jobs(rem_jobs, "[Candidate]")

    for (j in seq_along(rem)) {
      fits[[rem[j]]] <- rem_results[[j]]$fit
    }
  }

  # ===========================================================================
  # Rank candidates
  # ===========================================================================
  ll <- vapply(fits, function(z) z$log_lik_value, numeric(1))
  ord <- order(ll, decreasing = TRUE)
  fits <- fits[ord]
  grid <- grid[ord, , drop = FALSE]
  rownames(grid) <- NULL

  # ===========================================================================
  # STAGE 3: inference for best model only -- NO second optimization.
  # The selected candidate parameters are kept exactly as estimated in Stage 1/2.
  # .fit_zimp_raw(..., optimize = FALSE) only rebuilds the analytic gradient and
  # evaluates Hessian/vcov/SE at best0$par.
  # ===========================================================================
  best0 <- fits[[1L]]
  pbest <- if (is.na(grid$fit.var.power[1L])) NULL else grid$fit.var.power[1L]

  infer_fit <- .fit_zimp_raw(
    y = d$y, X = d$X, Z = d$Z, offset.count = d$offset.count,
    dist.family = grid$distribution[1L],
    zero_inflation = d$zero_inflation,
    var.power = pbest,
    ini.par = best0$par,
    progress = FALSE, trace = 0,
    eval.max = eval.max, iter.max = iter.max,
    rel.tol = rel.tol, x.tol = x.tol,
    n_gh = n_gh, gh = gh,
    compute.inference = TRUE,
    optimize = FALSE
  )

  best <- best0
  best$hessian <- infer_fit$hessian
  best$vcov <- infer_fit$vcov
  best$coef_table <- infer_fit$coef_table
  best$coef_table_internal <- infer_fit$coef_table_internal
  best$hessian_ok <- infer_fit$hessian_ok

  # User-facing metadata
  best$call <- match.call()
  best$interface <- d$interface
  best$formula <- d$formula
  best$count_formula <- d$count_formula
  best$zero_formula <- d$zero_formula
  best$count_intercept <- d$count_intercept
  best$zero_intercept <- d$zero_intercept
  best$zero_inflation <- d$zero_inflation
  best$mixing.distribution <- grid$distribution[1L]
  best$var.power <- pbest
  best$y <- d$y
  best$X <- d$X
  best$Z <- d$Z
  best$offset.count <- d$offset.count
  best$model_frame_count <- d$model_frame_count
  best$model_frame_zero <- d$model_frame_zero
  best$na_removed <- d$na_removed

  split_tab <- .zimp_split_coef_table(best, d)
  best$count_coef_table <- split_tab$count
  best$zero_coef_table <- split_tab$zero
  best$coefficient_names <- split_tab$names

  p_count <- ncol(d$X)
  best$phi_se <- if (!is.null(best$coef_table) &&
                     nrow(best$coef_table) >= p_count + 1L) {
    as.numeric(best$coef_table[p_count + 1L, "Std. Error"])
  } else {
    NA_real_
  }

  class(best) <- unique(c("zimp_fit", class(best)))

  # Internal ranked candidates for downstream IFM selection. No candidate is
  # re-optimized here; only metadata required by the dependence engine is added.
  for (ii in seq_along(fits)) {
    fits[[ii]]$call <- match.call()
    fits[[ii]]$interface <- d$interface
    fits[[ii]]$formula <- d$formula
    fits[[ii]]$count_formula <- d$count_formula
    fits[[ii]]$zero_formula <- d$zero_formula
    fits[[ii]]$count_intercept <- d$count_intercept
    fits[[ii]]$zero_intercept <- d$zero_intercept
    fits[[ii]]$zero_inflation <- d$zero_inflation
    fits[[ii]]$mixing.distribution <- grid$distribution[ii]
    fits[[ii]]$var.power <- if (is.na(grid$fit.var.power[ii])) NULL else grid$fit.var.power[ii]
    fits[[ii]]$y <- d$y
    fits[[ii]]$X <- d$X
    fits[[ii]]$Z <- d$Z
    fits[[ii]]$offset.count <- d$offset.count
    fits[[ii]]$model_frame_count <- d$model_frame_count
    fits[[ii]]$model_frame_zero <- d$model_frame_zero
    fits[[ii]]$na_removed <- d$na_removed
    class(fits[[ii]]) <- unique(c("zimp_fit", class(fits[[ii]])))
  }
  fits[[1L]] <- best

  summary_table <- data.frame(
    distribution = grid$distribution,
    var.power = grid$var.power,
    log_lik = vapply(fits, function(z) z$log_lik_value, numeric(1)),
    neg_log_lik = vapply(fits, function(z) z$neg_log_lik_value, numeric(1)),
    phi = vapply(fits, function(z) z$phi, numeric(1))
  )

  # Store selection output in a dedicated summary object.  Because it has its
  # own print method, `sel$summary` displays the selection report, while
  # `sel <- select_zimp(...)` remains completely silent.
  selection_summary <- list(
    summary_table = summary_table,
    best_model = best,
    call = match.call()
  )
  class(selection_summary) <- "zimp_selection_summary"

  # Compact public selection object. Candidate fits and the duplicated
  # candidate_grid are intentionally not retained after model selection.
  out <- list(
    summary_table = summary_table,
    best_model = best,
    summary = selection_summary,
    call = match.call()
  )
  if (isTRUE(.return.all)) out$all_results <- fits
  class(out) <- "zimp_selection"

  out
}

#' Print a select_zimp summary
#'
#' @param x An object of class `zimp_selection_summary`.
#' @param ... Additional arguments (currently unused).
#' @return The summary object, invisibly.
#' @export
print.zimp_selection_summary <- function(x, ...) {
  cat("\nModel selection completed\n")
  cat(strrep("=", 72), "\n", sep = "")
  print(x$summary_table, row.names = FALSE)

  best <- x$best_model
  cat("\nBest model\n")
  cat(strrep("=", 72), "\n", sep = "")
  cat("Model = ", if (isTRUE(best$zero_inflation))
    "Zero-Inflated Mixed Poisson" else "Mixed Poisson", "\n", sep = "")
  cat("Mixing distribution = ", best$mixing.distribution,
      " | var.power = ", best$var.power, "\n", sep = "")
  .zimp_print_fit(best, show.call = FALSE)

  invisible(x)
}
