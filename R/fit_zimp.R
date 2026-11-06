# =============================================================================
# PUBLIC API 1: fit_zimp() -- ONE SPECIFIED MODEL
# =============================================================================
#' Fit a Univariate Mixed Poisson or Zero-Inflated Mixed Poisson Model
#'
#' @description
#' Fits a single univariate Mixed Poisson (MP) or Zero-Inflated Mixed
#' Poisson (ZIMP) regression model for a specified mixing distribution and
#' variance-power parameterization.
#'
#' Both formula and matrix interfaces are supported. The count component
#' models the mean of the Poisson mixing distribution, while an optional
#' zero-inflation component models the probability that the count-generating
#' component is active.
#'
#' @details
#' `fit_zimp()` fits one specified model. It does not perform model selection
#' over mixing distributions or variance-power specifications. For model
#' selection, see [select_zimp()].
#'
#' The conditional count model is
#'
#' \deqn{Y \mid \Lambda \sim \mathrm{Poisson}(\Lambda),}
#'
#' where the distribution of the latent intensity \eqn{\Lambda} is determined
#' by `mixing.distribution`. The marginal mean is linked to the count
#' covariates through
#'
#' \deqn{\mu_i = \exp(x_i^\top \beta + o_i),}
#'
#' where \eqn{o_i} is an optional count offset.
#'
#' For a zero-inflated model, let \eqn{B_i} denote the indicator that the
#' count-generating component is active. The model uses
#'
#' \deqn{Y_i = B_i N_i,}
#'
#' with
#'
#' \deqn{\Pr(B_i = 1) = \pi_i,\qquad
#'       \mathrm{logit}(\pi_i) = z_i^\top \gamma.}
#'
#' Thus, `Z` models the probability of belonging to the count-generating
#' component rather than the structural-zero component.
#'
#' The dispersion parameter `phi` is constrained to be positive and is
#' optimized internally on the log scale.
#'
#' @section Model specification:
#'
#' Two input interfaces are available.
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
#' explicitly requested with `1`. For example,
#'
#' ```
#' y ~ x1 + x2
#' ```
#'
#' specifies a count model without an intercept, whereas
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
#' The response and design matrices can instead be supplied directly through
#' `Y`, `X`, and optionally `Z`. If `Z = NULL`, an MP model is fitted. If `Z`
#' is supplied, a ZIMP model is fitted. Intercept columns, when required, must
#' be included explicitly in `X` and `Z`.
#'
#' @section Mixing distributions:
#'
#' The following mixing distributions are supported:
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
#' The optional `var.power` argument controls the variance-power
#' parameterization. If `var.power = NULL`, the native/default
#' parameterization of the selected mixing distribution is used. Otherwise,
#' `var.power` must be one of `1`, `2`, or `3`.
#'
#' @section Initial values:
#'
#' Initial parameter values may be supplied through `ini.par`. Parameters are
#' ordered internally as
#'
#' ```
#' c(beta, log(phi), gamma)
#' ```
#'
#' where `beta` contains the count-component coefficients, `log(phi)` is the
#' log-dispersion parameter, and `gamma` contains the zero-component
#' coefficients. The `gamma` block is omitted for an MP model.
#'
#' @param formula Optional model formula. For an MP model, use a formula such
#'   as `y ~ 1 + x1 + x2`. For a ZIMP model, use a two-part formula such as
#'   `y ~ 1 + x1 + x2 | 1 + z1 + z2`. Intercepts must be included explicitly
#'   using `1`. When `formula` is supplied, `data` must also be supplied.
#'
#' @param data Optional data frame containing the variables referenced in
#'   `formula`.
#'
#' @param Y Numeric vector containing the count response for the matrix
#'   interface. Values must be finite, non-negative integers.
#'
#' @param X Numeric matrix or data frame containing the design variables for
#'   the count component when using the matrix interface. It must have
#'   `length(Y)` rows. Include a column of ones explicitly if an intercept is
#'   required.
#'
#' @param Z Optional numeric matrix or data frame containing the design
#'   variables for the zero-inflation component when using the matrix
#'   interface. It must have `length(Y)` rows. Include a column of ones
#'   explicitly if an intercept is required. If `NULL`, an MP model is fitted.
#'
#' @param offset.count Optional numeric offset for the count linear predictor.
#'   For the matrix interface, it must be a finite numeric vector of length
#'   `length(Y)`. If `NULL`, no count offset is used. With the formula
#'   interface, an offset may be specified in the count component using
#'   `offset()`.
#'
#' @param mixing.distribution Character string specifying the mixing
#'   distribution. One of `"gamma"`, `"lnorm"`, `"ig"`, `"weibull"`,
#'   `"loglogistic"`, or `"betaprime"`. The default is `"gamma"`.
#'
#' @param var.power Optional variance-power specification. If `NULL`, the
#'   native/default parameterization of the selected mixing distribution is
#'   used. Otherwise, it must be one of `1`, `2`, or `3`.
#'
#' @param n_gh Positive integer specifying the number of Gauss--Hermite
#'   quadrature nodes when Gauss--Hermite integration is required. The default
#'   is `20`.
#'
#' @param gh Optional precomputed Gauss--Hermite quadrature object. If `NULL`,
#'   the required quadrature rule is constructed internally.
#'
#' @param ini.par Optional numeric vector of initial parameter values. Its
#'   length must be `ncol(X) + 1` for an MP model and
#'   `ncol(X) + 1 + ncol(Z)` for a ZIMP model. The parameter order is
#'   `c(beta, log(phi), gamma)`.
#'
#' @param progress Logical indicating whether fitting progress should be
#'   displayed. The default is `TRUE`.
#'
#' @param trace Numeric optimizer trace level. The default is `0`.
#'
#' @param eval.max Positive integer giving the maximum number of objective
#'   function evaluations. The default is `100`.
#'
#' @param rel.tol Positive numeric value giving the relative convergence
#'   tolerance. The default is `1e-5`.
#'
#' @param x.tol Positive numeric value giving the parameter convergence
#'   tolerance. The default is `1e-5`.
#'
#' @param iter.max Positive integer giving the maximum number of optimizer
#'   iterations. The default is `100`.
#'
#' @return
#' An object of class `"zimp_fit"` containing the fitted model. Important
#' components include:
#'
#' \describe{
#'   \item{`be`}{Estimated regression coefficients \eqn{\hat\beta} for the
#'     count component.}
#'   \item{`phi`}{Estimated dispersion parameter on its natural scale.}
#'   \item{`gam`}{Estimated zero-component coefficients \eqn{\hat\gamma} for
#'     a ZIMP model.}
#'   \item{`log_lik_value`}{Maximized log-likelihood.}
#'   \item{`neg_log_lik_value`}{Negative log-likelihood evaluated at the
#'     fitted parameters.}
#'   \item{`AIC`}{Akaike information criterion.}
#'   \item{`BIC`}{Bayesian information criterion.}
#'   \item{`par`}{Estimated parameter vector on the internal optimization
#'     scale.}
#'   \item{`convergence`}{Optimizer convergence code.}
#'   \item{`message`}{Optimizer convergence message, when available.}
#'   \item{`hessian`}{Estimated Hessian matrix, when available.}
#'   \item{`vcov`}{Estimated covariance matrix, when available.}
#'   \item{`coef_table`}{Coefficient and inference table produced by the
#'     fitting engine.}
#'   \item{`count_coef_table`}{Coefficient table for the count component.}
#'   \item{`zero_coef_table`}{Coefficient table for the zero component for
#'     a ZIMP model.}
#'   \item{`phi_se`}{Standard error of the reported dispersion parameter,
#'     when available.}
#'   \item{`y`, `X`, `Z`}{Response and design matrices used for fitting.}
#'   \item{`offset.count`}{Count offset used in the fitted model.}
#'   \item{`mixing.distribution`}{Mixing distribution used in the model.}
#'   \item{`var.power`}{Variance-power specification used in the model.}
#'   \item{`zero_inflation`}{Logical indicating whether the fitted model
#'     contains a zero-inflation component.}
#'   \item{`nobs`}{Number of observations used in the fit.}
#'   \item{`npar`}{Number of estimated parameters.}
#'   \item{`call`}{Matched function call.}
#'   \item{`summary`}{Printable summary containing coefficient tables,
#'     dispersion, log-likelihood, AIC, BIC, and convergence information.}
#' }
#'
#' Additional components describing the formula, model frames, coefficient
#' names, interface, and missing-data handling may be included when relevant.
#'
#' @seealso
#' [select_zimp()] for selecting among multiple univariate MP or ZIMP model
#' specifications.
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
#' ## Gamma-mixed Poisson data are Negative Binomial marginally.
#' ## This provides a short reproducible example without relying on
#' ## non-exported simulation functions.
#' phi <- 0.2
#' y <- rnbinom(n, size = 1 / phi, mu = mu)
#'
#' dat <- data.frame(y = y, x1 = x1, x2 = x2)
#'
#' ## Formula interface
#' fit_formula <- fit_zimp(
#'   y ~ 1 + x1 + x2,
#'   data = dat,
#'   mixing.distribution = "gamma",
#'   var.power = 2,
#'   progress = FALSE
#' )
#'
#' fit_formula
#'
#' ## Equivalent matrix interface
#' fit_matrix <- fit_zimp(
#'   Y = y,
#'   X = X,
#'   mixing.distribution = "gamma",
#'   var.power = 2,
#'   progress = FALSE
#' )
#'
#' @export
fit_zimp <- function(formula = NULL, data = NULL,
                     Y = NULL, X = NULL, Z = NULL,
                     offset.count = NULL,
                     mixing.distribution = "gamma",
                     var.power = NULL,
                     n_gh = 20, gh = NULL, ini.par = NULL,
                     progress = TRUE, trace = 0,
                     eval.max = 100, rel.tol = 1e-5,
                     x.tol = 1e-5, iter.max = 100) {

  d <- .zimp_resolve_input(formula, data, Y, X, Z, offset.count)

  choices <- c("gamma","lnorm","ig","weibull","loglogistic","betaprime")
  mixing.distribution <- match.arg(mixing.distribution, choices)

  if (!is.null(var.power) &&
      (length(var.power) != 1L || !is.numeric(var.power) || !is.finite(var.power) ||
       !(var.power %in% c(1,2,3))))
    stop("fit_zimp() requires 'var.power' to be NULL, 1, 2, or 3.")

  # Pareto has only its native/default variance parameterization.
  # Ignore a user-supplied variance power rather than passing it to the
  # Pareto likelihood, which does not accept var_power.
  if (identical(mixing.distribution, "pareto") && !is.null(var.power)) {
    if (isTRUE(progress)) {
      cat(
        "Note: Pareto always uses its default variance power (var.power = 2); ",
        "the user-specified var.power is ignored.\n",
        sep = ""
      )
    }
    var.power <- NULL
  }

  p <- ncol(d$X)
  pz <- if (d$zero_inflation) ncol(d$Z) else 0L
  if (!is.null(ini.par) && length(ini.par) != p + 1L + pz)
    stop("Incorrect length of 'ini.par'.")

  fit <- .fit_zimp_raw(
    y = d$y, X = d$X, Z = d$Z, offset.count = d$offset.count,
    dist.family = mixing.distribution,
    zero_inflation = d$zero_inflation,
    var.power = var.power,
    ini.par = ini.par, progress = FALSE, trace = trace,
    eval.max = eval.max, iter.max = iter.max,
    rel.tol = rel.tol, x.tol = x.tol, n_gh = n_gh, gh = gh
  )

  fit$call <- match.call()
  fit$interface <- d$interface
  fit$formula <- d$formula
  fit$count_formula <- d$count_formula
  fit$zero_formula <- d$zero_formula
  fit$count_intercept <- d$count_intercept
  fit$zero_intercept <- d$zero_intercept
  fit$zero_inflation <- d$zero_inflation
  fit$mixing.distribution <- mixing.distribution
  fit$var.power <- var.power
  fit$y <- d$y; fit$X <- d$X; fit$Z <- d$Z; fit$offset.count <- d$offset.count
  fit$model_frame_count <- d$model_frame_count
  fit$model_frame_zero <- d$model_frame_zero
  fit$na_removed <- d$na_removed

  # Information criteria.
  fit$nobs <- length(d$y)
  fit$npar <- length(fit$par)
  fit$AIC <- -2 * fit$log_lik_value + 2 * fit$npar
  fit$BIC <- -2 * fit$log_lik_value + log(fit$nobs) * fit$npar

  # User-facing coefficient tables with meaningful variable names.
  split_tab <- .zimp_split_coef_table(fit, d)
  fit$count_coef_table <- split_tab$count
  fit$zero_coef_table <- split_tab$zero
  fit$coefficient_names <- split_tab$names

  # phi is the (p + 1)-th row of the engine's report-scale coefficient table.
  # Keep its delta-method standard error separately for clean console output.
  p_count <- ncol(d$X)
  fit$phi_se <- if (!is.null(fit$coef_table) &&
                    nrow(fit$coef_table) >= p_count + 1L) {
    as.numeric(fit$coef_table[p_count + 1L, "Std. Error"])
  } else {
    NA_real_
  }

  # Store the user-facing display in a dedicated summary object.
  fit$summary <- .zimp_make_summary(fit)

  class(fit) <- unique(c("zimp_fit", class(fit)))
  fit
}


# Construct the printable summary stored in fit$summary.
.zimp_make_summary <- function(fit) {
  out <- list(
    call = fit$call,
    count_coef_table = fit$count_coef_table,
    zero_coef_table = fit$zero_coef_table,
    phi = fit$phi,
    phi_se = fit$phi_se,
    log_lik_value = fit$log_lik_value,
    AIC = fit$AIC,
    BIC = fit$BIC,
    convergence = fit$convergence,
    zero_inflation = fit$zero_inflation
  )
  class(out) <- "zimp_summary"
  out
}


# Printing a fitted object (including an unassigned direct call) displays only
# its stored summary. Assignment itself remains silent under normal R semantics.
#' @export
print.zimp_fit <- function(x, ...) {
  print(x$summary, ...)
  invisible(x)
}


#' @export
print.zimp_summary <- function(x, ...) {

  # Match the original fit_zimp() console display exactly, except that the
  # former "Starting model fitting / Model = ..." progress banner is omitted.
  print_coef_block <- function(tab) {
    if (is.null(tab)) {
      cat("Inference unavailable.\n")
      return(invisible(NULL))
    }
    print(tab, quote = FALSE)
    invisible(NULL)
  }

  cat("Call:\n")
  print(x$call)

  cat("\nCount model coefficients:\n")
  print_coef_block(x$count_coef_table)

  if (isTRUE(x$zero_inflation)) {
    cat("\nZero-inflation model coefficients:\n")
    print_coef_block(x$zero_coef_table)
  }

  cat("---\n")
  cat("Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 ' ' 1\n")

  phi_se <- x$phi_se
  if (is.null(phi_se) || length(phi_se) != 1L || !is.finite(phi_se)) {
    phi_se_txt <- "NA"
  } else {
    phi_se_txt <- formatC(phi_se, digits = 5, format = "fg")
  }

  cat("\nDispersion parameter:\n")
  cat("phi = ", formatC(x$phi, digits = 6, format = "fg"),
      "  (Std. Error = ", phi_se_txt, ")\n", sep = "")

  cat("\nLog-likelihood: ",
      formatC(x$log_lik_value, digits = 6, format = "f"), sep = "")
  cat("\nAIC: ", formatC(x$AIC, digits = 2, format = "f"), sep = "")
  cat("\nBIC: ", formatC(x$BIC, digits = 2, format = "f"), "\n", sep = "")

  cat("\nConvergence: ", ifelse(x$convergence == 0, "Yes", "No"), "\n", sep = "")

  invisible(x)
}
