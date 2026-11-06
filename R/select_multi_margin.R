# ============================================================
# select_multi_margin
# Fast staged parallel marginal selection.
#
# Stage 1: fit direct models, benchmark, and required Gamma anchors in parallel.
# Stage 2: fit warm-start families in parallel.
#   * lnorm, gamma, ig: ordinary distribution-specific initialization.
#   * weibull, loglogistic, betaprime: target-specific initialization for phi,
#     but beta/gamma are replaced by fitted Gamma beta/gamma at the SAME var.power.
#
# Requires current zimp helpers:
#   .fit_zimp_raw(), .zimp_init_par(), I_gamma(), I_lnorm(), ...
# ============================================================


# -----------------------------------------------------------------------------
# Stand-alone worker helper for optional I-cache computation.
# IMPORTANT: this function is deliberately defined outside select_multi_margin()
# so a PSOCK worker does not receive the enclosing estimation environment.
# Each call receives only one margin's compact task.
# -----------------------------------------------------------------------------
.mzimp_compute_I_margin_worker <- function(task) {
  mj <- task$margin_fit
  dist <- mj$distribution

  get_fun <- function(nm) {
    if (exists(nm, mode="function", inherits=TRUE))
      get(nm, mode="function", inherits=TRUE)
    else
      getFromNamespace(nm, "mzimp")
  }

  one <- function(z) {
    nm <- switch(dist,
      gamma="I_gamma", lnorm="I_lnorm", ig="I_ig", weibull="I_weibull", loglogistic="I_loglogistic", betaprime="I_betaprime")
    fun <- get_fun(nm)
    args <- list(z=z, X=task$X, be=mj$be, phi=mj$phi, d=task$I.d0)
    fml <- names(formals(fun))
    if ("var_power" %in% fml) args$var_power <- mj$var.power
    if ("var.power" %in% fml) args$var.power <- mj$var.power
    if ("n_gh" %in% fml) args$n_gh <- task$n_gh
    if ("gh" %in% fml) args$gh <- task$gh
    ans <- as.matrix(do.call(fun,args))
    dimnames(ans) <- NULL
    ans
  }

  y <- task$y
  Iy <- one(y)
  I0 <- one(rep.int(0L, length(y)))

  if (isTRUE(mj$zero_inflation)) {
    if (is.null(task$Z) || is.null(mj$gam))
      stop("Missing Z/gamma for zero-inflated margin ", mj$margin, ".")
    pb <- unname(plogis(as.vector(task$Z %*% mj$gam)))
    y0 <- unname(y == 0L)
  } else {
    pb <- unname(rep.int(1, length(y)))
    y0 <- NULL
  }
  list(I_y=Iy, I_0=I0, pb=pb, y0=y0)
}


#' Select Marginal Models for Multivariate Count Data
#'
#' @description
#' Fits and selects a univariate Mixed Poisson (MP) or Zero-Inflated
#' Mixed Poisson (ZIMP) model independently for each margin of a
#' multivariate count response.
#'
#' Candidate models may differ in mixing distribution, variance-power
#' specification, and, when requested, the presence or absence of zero
#' inflation. The best model for each margin is selected using AIC or BIC.
#'
#' `select_multi_margin()` is designed primarily as the marginal-selection
#' stage preceding [fit_multizimp()]. It supports common or
#' margin-specific count and zero-component designs, serial or parallel
#' computation, and optional precomputation of the checkerboard
#' interval-probability matrices required by the multivariate dependence
#' model.
#'
#' @details
#' Let
#'
#' \deqn{\mathbf{Y}=(Y_1,\ldots,Y_k)}
#'
#' denote an \eqn{n \times k} multivariate count response. The function
#' fits a collection of candidate univariate models separately to each
#' margin \eqn{Y_j} and selects the best marginal model according to
#' `selection.criterion`.
#'
#' Marginal selection does not estimate dependence between response
#' margins. Dependence estimation is performed subsequently by
#' [fit_multizimp()].
#'
#' @section Response and design specification:
#'
#' `Y` must be an \eqn{n \times k} numeric matrix or data frame containing
#' finite, non-negative integer counts. Each column represents one response
#' margin.
#'
#' The count-component design `X` can be supplied either as:
#'
#' \itemize{
#'   \item one matrix or data frame shared by all margins; or
#'   \item a list of length `k` containing a separate design for each margin.
#' }
#'
#' When a common `X` is supplied, `X.columns` may be used to specify a
#' different subset of columns for each margin.
#'
#' The zero-component design `Z` follows the same convention. It may be
#' `NULL`, a common matrix/data frame, or a list of length `k`.
#'
#' No intercept is added automatically to either `X` or `Z`. Include an
#' intercept column explicitly when required.
#'
#' Data-frame designs may contain numeric, integer, logical, factor, and
#' character variables. Character variables are converted to factors and
#' categorical variables are expanded internally into numeric design
#' columns.
#'
#' @section MP and ZIMP margins:
#'
#' The interpretation of `Z` determines whether a margin is fitted as an
#' ordinary MP model or as a ZIMP model.
#'
#' \describe{
#'   \item{`Z = NULL`, `check.zero = FALSE`}{
#'     All margins are fitted as ordinary MP models.
#'   }
#'
#'   \item{Common matrix/data frame `Z`}{
#'     All margins are fitted as ZIMP models using the supplied
#'     zero-component design.
#'   }
#'
#'   \item{`Z` supplied as a list}{
#'     Each margin is specified independently. `Z[[j]] = NULL` gives an MP
#'     margin, whereas a non-`NULL` `Z[[j]]` gives a ZIMP margin.
#'   }
#'
#'   \item{`Z = NULL`, `check.zero = TRUE`}{
#'     For every marginal distribution and variance-power specification,
#'     both an MP model and an intercept-only ZIMP model are considered.
#'     The presence or absence of zero inflation is therefore selected
#'     independently for each margin.
#'   }
#' }
#'
#' If `Z` is explicitly supplied, `check.zero` is ignored.
#'
#' For a ZIMP margin, the zero-component regression models the Bernoulli
#' activation probability. If \eqn{B_i} denotes the activation indicator,
#'
#' \deqn{
#' \Pr(B_i=1)=\operatorname{logit}^{-1}(z_i^\top\gamma).
#' }
#'
#' Thus, \eqn{B_i=1} corresponds to the count-generating component and
#' \eqn{B_i=0} to a structural zero.
#'
#' @section Candidate marginal models:
#'
#' The following mixing distributions are supported:
#'
#' \describe{
#'   \item{`"gamma"`}{Gamma mixing distribution.}
#'   \item{`"lnorm"`}{Lognormal mixing distribution.}
#'   \item{`"ig"`}{Inverse-Gaussian mixing distribution.}
#'   \item{`"weibull"`}{Weibull mixing distribution.}
#'   \item{`"loglogistic"`}{Log-logistic mixing distribution.}
#'   \item{`"betaprime"`}{Beta-prime mixing distribution.}
#' }
#'
#' If `var.power` is supplied, every requested mixing distribution is
#' crossed with every requested variance power. For example,
#'
#' ```
#' mixing.distribution = c("gamma", "lnorm")
#' var.power = c(1, 2, 3)
#' ```
#'
#' defines six distribution--variance-power specifications before any
#' MP/ZIMP expansion.
#'
#' If `var.power = NULL`, each distribution uses its native/default
#' variance-power parameterization:
#'
#' \describe{
#'   \item{Gamma}{2}
#'   \item{Lognormal}{2}
#'   \item{Inverse Gaussian}{3}
#'   \item{Weibull}{2}
#'   \item{Log-logistic}{2}
#'   \item{Beta-prime}{2}
#' }
#'
#' The complete candidate specifications used for marginal selection are
#' returned in `candidate_grid`.
#'
#' @section Marginal model selection:
#'
#' Candidate models are fitted independently for each response margin.
#' The selected model minimizes the information criterion specified by
#' `selection.criterion`.
#'
#' For `selection.criterion = "BIC"`,
#'
#' \deqn{
#' \mathrm{BIC}=-2\ell+\log(n)q,
#' }
#'
#' whereas for `selection.criterion = "AIC"`,
#'
#' \deqn{
#' \mathrm{AIC}=-2\ell+2q,
#' }
#'
#' where \eqn{\ell} is the maximized marginal log-likelihood and \eqn{q}
#' is the number of estimated marginal parameters.
#'
#' Selection is performed separately for every margin. Consequently,
#' different response margins may select different mixing distributions,
#' variance powers, and, when `check.zero = TRUE`, different MP/ZIMP
#' specifications.
#'
#' Candidate fitting uses the same initialization strategy as
#' [select_zimp()]. Gamma, lognormal, and inverse-Gaussian models use
#' distribution-specific initialization. Weibull, log-logistic, and
#' beta-prime models use fitted Gamma models with the same effective
#' variance power to initialize their regression coefficients, while
#' retaining distribution-specific initialization for dispersion.
#'
#' Gamma models required only for initialization are internal fits and
#' do not participate in final model selection unless the corresponding
#' Gamma specification was explicitly requested.
#'
#' @section Benchmark model:
#'
#' `benchmark.distribution` specifies a reference marginal distribution
#' used to compare the selected marginal models. Supported benchmark
#' distributions are `"lnorm"`, `"gamma"`, and `"ig"`.
#'
#' The benchmark is fitted using its native/default variance-power
#' parameterization under the corresponding zero-inflation specification.
#' It is used only for reporting and does not affect marginal model
#' selection.
#'
#' The returned object contains the benchmark log-likelihood for each
#' margin and the log-likelihood gain of the selected model relative to
#' the benchmark.
#'
#' @section Parallel computation:
#'
#' Candidate marginal models can be fitted in parallel.
#'
#' If `parallel = NULL` or `parallel = 1`, fitting is serial. If
#' `parallel` is an integer greater than one, a PSOCK cluster is created
#' and marginal-model fitting jobs are distributed across workers.
#'
#' Parallelization is across independent margin--candidate fitting jobs.
#' The number of workers actually used is bounded by the number of
#' available fitting tasks.
#'
#' `block_size` optionally controls how fitting jobs are grouped for
#' load-balanced parallel execution. If `block_size = NULL`, an
#' appropriate block size is selected automatically.
#'
#' Serial and parallel execution use the same model-selection procedure
#' and return the same public result structure.
#'
#' @section Checkerboard probability cache:
#'
#' If `compute.I = TRUE`, checkerboard interval probabilities are
#' precomputed for each selected marginal model for subsequent use by
#' [fit_multizimp()].
#'
#' For each margin \eqn{j}, the cache contains:
#'
#' \describe{
#'   \item{`I_y_list[[j]]`}{
#'     An \eqn{n \times d} matrix of interval probabilities evaluated at
#'     the observed counts, where \eqn{d} is `I.d0`.
#'   }
#'
#'   \item{`I_0_list[[j]]`}{
#'     An \eqn{n \times d} matrix of interval probabilities evaluated at
#'     zero.
#'   }
#'
#'   \item{`pb_list[[j]]`}{
#'     Bernoulli activation probabilities. For an MP margin this is a
#'     vector of ones. For a ZIMP margin it is
#'     \eqn{\operatorname{logit}^{-1}(Z_j\gamma_j)}.
#'   }
#'
#'   \item{`y0_list[[j]]`}{
#'     Logical indicator of observed zeros for a ZIMP margin. It is `NULL`
#'     for an MP margin.
#'   }
#' }
#'
#' The checkerboard resolution \eqn{d} is controlled by `I.d0`.
#'
#' If `compute.I = FALSE`, no cache is created and `I_d0`, `I_y_list`,
#' `I_0_list`, `pb_list`, and `y0_list` are returned as `NULL`.
#'
#' @param Y An \eqn{n \times k} numeric matrix or data frame containing
#'   multivariate count responses. All entries must be finite,
#'   non-negative, and integer-valued. Columns represent response margins.
#'
#' @param X Count-component design. Either a numeric matrix/data frame
#'   shared by all margins or a list of length `k` containing
#'   margin-specific design matrices/data frames. Intercepts must be
#'   supplied explicitly.
#'
#' @param Z Optional zero-component design. May be `NULL`, a matrix/data
#'   frame shared by all margins, or a list of length `k`. For list input,
#'   `Z[[j]] = NULL` specifies an MP margin and a non-`NULL` element
#'   specifies a ZIMP margin.
#'
#' @param X.columns Optional list of length `k` specifying which columns
#'   of a common `X` are used for each margin. Elements may contain column
#'   names or column indices. A `NULL` element uses all columns. Ignored
#'   when `X` is supplied as a list.
#'
#' @param Z.columns Optional list of length `k` specifying which columns
#'   of a common `Z` are used for each margin. Elements may contain column
#'   names or column indices. A `NULL` element uses all columns. Ignored
#'   when `Z` is supplied as a list.
#'
#' @param check.zero Logical. If `TRUE` and `Z = NULL`, both MP and
#'   intercept-only ZIMP specifications are considered for every candidate
#'   marginal model. If `Z` is explicitly supplied, `check.zero` is
#'   ignored. The default is `FALSE`.
#'
#' @param selection.criterion Character string specifying the information
#'   criterion used to select the best model independently for each margin.
#'   Either `"BIC"` or `"AIC"`. The default is `"BIC"`.
#'
#' @param mixing.distribution Character vector specifying candidate mixing
#'   distributions. Supported values are `"gamma"`, `"lnorm"`, `"ig"`,
#'   `"weibull"`, `"loglogistic"`, and `"betaprime"`. By default, all
#'   supported distributions are considered.
#'
#' @param var.power `NULL` or a numeric vector containing candidate
#'   variance powers. Values must be selected from `1`, `2`, and `3`.
#'   If `NULL`, each distribution uses its native/default variance-power
#'   parameterization.
#'
#' @param benchmark.distribution Character string specifying the benchmark
#'   marginal distribution used for reporting. Supported values are
#'   `"lnorm"`, `"gamma"`, and `"ig"`. The default is `"lnorm"`.
#'
#' @param ini.par Optional initial parameter specification. It may be a
#'   list of length `k` containing one initial parameter vector per margin,
#'   or a named list indexed by mixing distribution whose elements are
#'   themselves lists of length `k`.
#'
#' @param parallel `NULL` or a positive integer specifying the requested
#'   number of PSOCK workers. `NULL` or `1` uses serial execution.
#'
#' @param block_size `NULL` or a positive integer controlling the number
#'   of fitting jobs grouped into each load-balanced parallel block. If
#'   `NULL`, the block size is chosen automatically.
#'
#' @param compute.I Logical. If `TRUE`, precompute checkerboard
#'   interval-probability matrices for the selected marginal models.
#'   The default is `FALSE`.
#'
#' @param I.d0 Integer greater than or equal to `2` specifying the
#'   checkerboard grid size when `compute.I = TRUE`. The default is `10L`.
#'
#' @param progress Deprecated and ignored. Retained for API compatibility.
#'   `select_multi_margin()` is silent.
#'
#' @param show.result Deprecated and ignored. Retained for API
#'   compatibility. Selection results are stored in the returned object.
#'
#' @param trace Deprecated and ignored. Retained for API compatibility.
#'   `select_multi_margin()` is silent.
#'
#' @param eval.max Positive integer giving the maximum number of
#'   objective-function evaluations for each marginal fit.
#'
#' @param iter.max Positive integer giving the maximum number of optimizer
#'   iterations for each marginal fit.
#'
#' @param rel.tol Positive numeric value giving the relative convergence
#'   tolerance for marginal fitting.
#'
#' @param x.tol Positive numeric value giving the parameter convergence
#'   tolerance for marginal fitting.
#'
#' @param n_gh Positive integer giving the number of Gauss--Hermite
#'   quadrature nodes used when required by marginal calculations.
#'   The default is `20`.
#'
#' @param gh Optional precomputed Gauss--Hermite quadrature object. If
#'   `NULL`, the required quadrature rule is constructed internally.
#'
#' @param seed Optional integer seed used to initialize reproducible
#'   random-number streams in parallel execution.
#'
#' @return
#' An object of class `"multi_margin_fit"` containing the marginal-selection
#' results and the information required by [fit_multizimp()]. The public
#' components are:
#'
#' \describe{
#'   \item{`Y`}{
#'     Validated \eqn{n \times k} response matrix.
#'   }
#'
#'   \item{`X`, `Z`}{
#'     Count- and zero-component designs after internal conversion.
#'   }
#'
#'   \item{`X_is_list`, `Z_is_list`}{
#'     Logical indicators identifying whether the corresponding designs are
#'     stored as margin-specific lists.
#'   }
#'
#'   \item{`X_columns`, `Z_columns`}{
#'     Margin-specific column selections when common designs were supplied.
#'   }
#'
#'   \item{`selected_margins`}{
#'     List of length `k` containing the selected marginal models. Each
#'     element has class `"mzimp_margin_fit"` and contains, among other
#'     quantities, the selected mixing distribution, variance power,
#'     zero-inflation status, count coefficients (`be`), dispersion
#'     parameter (`phi`), zero-component coefficients (`gam`) when
#'     applicable, and maximized marginal log-likelihood.
#'   }
#'
#'   \item{`benchmark.distribution`}{
#'     Mixing distribution used as the benchmark.
#'   }
#'
#'   \item{`benchmark_by_margin`}{
#'     Benchmark-fit information for individual response margins.
#'   }
#'
#'   \item{`benchmark`}{
#'     Summary of the comparison between the selected models and the
#'     benchmark, including marginal and aggregate log-likelihood gains.
#'   }
#'
#'   \item{`selected_summary`}{
#'     Data frame with one row per response margin describing the selected
#'     distribution, variance power, zero-inflation status, log-likelihood,
#'     AIC, BIC, benchmark log-likelihood, and gain relative to the
#'     benchmark.
#'   }
#'
#'   \item{`I_precomputed`}{
#'     Logical indicating whether the checkerboard probability cache was
#'     computed.
#'   }
#'
#'   \item{`I_d0`}{
#'     Checkerboard grid size used for the cache, or `NULL` when
#'     `compute.I = FALSE`.
#'   }
#'
#'   \item{`I_y_list`, `I_0_list`}{
#'     Lists of checkerboard interval-probability matrices evaluated at
#'     observed counts and at zero, respectively, when `compute.I = TRUE`;
#'     otherwise `NULL`.
#'   }
#'
#'   \item{`pb_list`}{
#'     List of Bernoulli activation-probability vectors when
#'     `compute.I = TRUE`; otherwise `NULL`.
#'   }
#'
#'   \item{`y0_list`}{
#'     List of observed-zero indicators for ZIMP margins when
#'     `compute.I = TRUE`. Elements corresponding to MP margins are `NULL`.
#'   }
#'
#'   \item{`mixing.distribution`, `var.power`}{
#'     Candidate marginal specifications requested by the user.
#'   }
#'
#'   \item{`candidate_grid`}{
#'     Data frame describing the candidate marginal specifications
#'     participating in model selection.
#'   }
#'
#'   \item{`check.zero`}{
#'     Logical indicating whether MP versus intercept-only ZIMP comparison
#'     was requested when `Z = NULL`.
#'   }
#'
#'   \item{`selection.criterion`}{
#'     Information criterion used for marginal model selection.
#'   }
#'
#'   \item{`n`, `k`}{
#'     Number of observations and response margins.
#'   }
#'
#'   \item{`response_names`}{
#'     Names of the response margins.
#'   }
#'
#'   \item{`parallel`, `n_workers`}{
#'     Requested parallel setting and the actual number of workers used.
#'   }
#'
#'   \item{`performance`}{
#'     Timing and, when available, memory-use information for marginal
#'     fitting and optional checkerboard-cache computation.
#'   }
#'
#'   \item{`summary`}{
#'     A regular list containing `selected_summary`, benchmark information,
#'     selection criterion, response dimensions and names, and performance
#'     information.
#'   }
#' }
#'
#' @seealso
#' [fit_multizimp()] for fitting the multivariate dependence model after
#' marginal selection, [select_zimp()] for selecting a univariate marginal
#' model, and [fit_zimp()] for fitting one specified univariate model.
#'
#' @examples
#' set.seed(12345)
#'
#' n <- 30
#' k <- 3
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#'
#' X <- cbind(Intercept = 1, x1 = x1, x2 = x2)
#'
#' Y <- cbind(
#'   y1 = rpois(n, exp(0.5 + 0.3 * x1 + 0.2 * x2)),
#'   y2 = rpois(n, exp(0.3 - 0.2 * x1 + 0.3 * x2)),
#'   y3 = rpois(n, exp(0.4 + 0.1 * x1 - 0.2 * x2))
#' )
#'
#' ## Select among Gamma and lognormal marginal models
#' fit <- select_multi_margin(
#'   Y = Y,
#'   X = X,
#'   Z = NULL,
#'   mixing.distribution = c("gamma", "lnorm"),
#'   var.power = 2,
#'   selection.criterion = "BIC",
#'   parallel = 1
#' )
#'
#' ## One selected model per response margin
#' fit$selected_summary
#'
#' ## Selected marginal model objects
#' fit$selected_margins
#'
#' @export

select_multi_margin <- function(
    Y,
    X,
    Z = NULL,
    X.columns = NULL,
    Z.columns = NULL,
    check.zero = FALSE,
    selection.criterion = c("BIC", "AIC"),
    mixing.distribution = c(
      "gamma", "lnorm", "ig",
      "weibull", "loglogistic", "betaprime"
    ),
    var.power = NULL,
    benchmark.distribution = "lnorm",
    ini.par = NULL,
    parallel = NULL,
    block_size = NULL,
    compute.I = FALSE,
    I.d0 = 10L,
    progress = FALSE,
    show.result = FALSE,
    trace = 0,
    eval.max = 100,
    iter.max = 100,
    rel.tol = 1e-5,
    x.tol = 1e-5,
    n_gh = 20,
    gh = NULL,
    seed = NULL
) {
  # Selection is intentionally silent; progress/show.result are retained only for API compatibility.
  progress <- FALSE
  show.result <- FALSE
  trace <- 0
  selection.criterion <- match.arg(selection.criterion)
  allowed_dist <- c("gamma","lnorm","ig","weibull","loglogistic","betaprime")
  mixing.distribution <- unique(match.arg(mixing.distribution, allowed_dist, several.ok=TRUE))
  benchmark.distribution <- match.arg(benchmark.distribution, allowed_dist)

  # Y: matrix/data.frame of count responses.
  if (!(is.matrix(Y) || is.data.frame(Y)))
    stop("Y must be an n x k matrix/data.frame of nonnegative integer counts.")
  yn <- colnames(Y)
  if (is.data.frame(Y) &&
      any(!vapply(Y, function(z) is.numeric(z) || is.integer(z), logical(1))))
    stop("Every column of Y must be a numeric/integer count variable.")
  Y <- as.matrix(Y)
  if (!is.numeric(Y) || any(!is.finite(Y)) || any(Y < 0) ||
      any(abs(Y-round(Y)) > 1e-8))
    stop("Every column of Y must contain finite nonnegative integer counts.")
  storage.mode(Y) <- "integer"
  n <- nrow(Y); k <- ncol(Y)
  if (n < 1L || k < 1L) stop("Y must have positive dimensions.")
  if (is.null(yn)) yn <- paste0("Y",seq_len(k))
  colnames(Y) <- yn

  # Convert one X/Z data.frame to a numeric design matrix.
  # No intercept is ever added automatically.
  # Numeric-only data.frame -> its columns unchanged.
  # Factor/character columns -> one-hot encoding.  If there is no explicit
  # intercept column, all factor levels are retained.  Thus the coefficient
  # of the chosen/reference level is represented by its own dummy column.
  # If an explicit all-ones column is supplied, treatment coding is used.
  design_df <- function(A,name) {
    if (nrow(A)!=n) stop(name," must have nrow(Y) rows.")
    if (!ncol(A)) stop(name," must contain at least one column.")
    A[] <- lapply(A,function(z) if(is.character(z)) factor(z) else z)
    ok <- vapply(A,function(z) is.numeric(z)||is.integer(z)||is.logical(z)||is.factor(z),logical(1))
    if(any(!ok)) stop(name," contains unsupported column type(s): ",paste(names(A)[!ok],collapse=", "))
    if(anyNA(A)) stop(name," contains NA.")

    has_int <- any(vapply(A,function(z)
      (is.numeric(z)||is.integer(z)||is.logical(z)) &&
      all(is.finite(as.numeric(z))) && all(abs(as.numeric(z)-1)<1e-12),logical(1)))

    ans <- list()
    for(nm in names(A)) {
      z <- A[[nm]]
      if(is.factor(z)) {
        if(has_int) {
          M <- model.matrix(~z)[,-1L,drop=FALSE]
          if(ncol(M)) colnames(M) <- paste0(nm,levels(z)[-1L])
        } else {
          M <- model.matrix(~z-1)
          colnames(M) <- paste0(nm,levels(z))
        }
        if(ncol(M)) ans[[length(ans)+1L]] <- M
      } else {
        v <- as.numeric(z)
        if(any(!is.finite(v))) stop(name,"$",nm," contains NA/NaN/Inf.")
        ans[[length(ans)+1L]] <- matrix(v,ncol=1L,dimnames=list(NULL,nm))
      }
    }
    M <- do.call(cbind,ans)
    storage.mode(M) <- "double"
    M
  }

  design_matrix <- function(A,name) {
    if(is.data.frame(A)) return(design_df(A,name))
    if(!is.matrix(A)) stop(name," must be a matrix/data.frame.")
    if(nrow(A)!=n) stop(name," must have nrow(Y) rows.")
    if(!ncol(A) || !is.numeric(A) || any(!is.finite(A)))
      stop(name," matrix must be a finite numeric matrix with at least one column.")
    A <- as.matrix(A); storage.mode(A)<-"double"
    if(is.null(colnames(A))) colnames(A)<-paste0(name,seq_len(ncol(A)))
    A
  }
  if (!is.logical(compute.I) || length(compute.I)!=1L || is.na(compute.I)) stop("`compute.I` must be TRUE or FALSE.")
  I.d0 <- as.integer(round(I.d0)); if (!is.finite(I.d0) || I.d0 < 2L) stop("`I.d0` must be >= 2.")
  if (!is.logical(check.zero) || length(check.zero)!=1L || is.na(check.zero)) stop("check.zero must be TRUE or FALSE.")
  if (!is.null(Z) && isTRUE(check.zero)) warning("Z is supplied, so check.zero is ignored.")

  msg <- function(...) if (isTRUE(progress)) { cat(...,"\n",sep=""); flush.console() }

  ram_snapshot_MB <- function(cl=NULL) {
    if (!requireNamespace("ps",quietly=TRUE))
      return(list(master=NA_real_,workers=NA_real_,total=NA_real_,max_worker=NA_real_))
    master <- as.numeric(ps::ps_memory_info()[["rss"]]) / 2^20
    workers <- numeric(0)
    if (!is.null(cl)) {
      z <- try(parallel::clusterCall(cl,function() {
        if (!requireNamespace("ps",quietly=TRUE)) return(NA_real_)
        as.numeric(ps::ps_memory_info()[["rss"]]) / 2^20
      }),silent=TRUE)
      if (!inherits(z,"try-error")) workers <- unlist(z,use.names=FALSE)
    }
    list(master=master,workers=sum(workers,na.rm=TRUE),
         total=master+sum(workers,na.rm=TRUE),
         max_worker=if(length(workers)&&any(is.finite(workers))) max(workers,na.rm=TRUE) else NA_real_)
  }
  default_vp <- c(gamma=2L, lnorm=2L, ig=3L, weibull=2L, loglogistic=2L, betaprime=2L)
  mixed_Z_input <- is.list(Z) && length(Z)==k && any(vapply(Z,is.null,logical(1))) && any(!vapply(Z,is.null,logical(1)))

  if (!is.null(var.power)) {
    if (!is.numeric(var.power) || any(!is.finite(var.power)) ||
        any(abs(var.power-round(var.power))>1e-12) || any(!(as.integer(round(var.power)) %in% 1:3)))
      stop("var.power must be NULL or contain only 1, 2, or 3.")
    var.power <- unique(as.integer(round(var.power)))
  }

  # Selection candidate grid. Native/default is represented by fit.var.power=NA.
  if (is.null(var.power)) {
    candidate_grid <- data.frame(
      distribution=mixing.distribution,
      fit.var.power=NA_integer_,
      var.power=unname(default_vp[mixing.distribution]),
      stringsAsFactors=FALSE
    )
  } else {
    parts <- list()
    nonpareto <- setdiff(mixing.distribution,"pareto")
    if (length(nonpareto)) {
      g <- expand.grid(distribution=nonpareto, fit.var.power=var.power,
                       KEEP.OUT.ATTRS=FALSE, stringsAsFactors=FALSE)
      g$var.power <- g$fit.var.power; parts[[length(parts)+1L]] <- g
    }
    candidate_grid <- do.call(rbind,parts); rownames(candidate_grid) <- NULL
  }
  if (!nrow(candidate_grid)) stop("No admissible candidate remains.")

  if (is.null(Z) && isTRUE(check.zero)) {
    candidate_grid <- merge(candidate_grid,data.frame(zero_inflation=c(FALSE,TRUE)),by=NULL)
  } else if (mixed_Z_input) {
    candidate_grid <- merge(candidate_grid,data.frame(zero_inflation=c(FALSE,TRUE)),by=NULL)
  } else {
    candidate_grid$zero_inflation <- !is.null(Z)
  }
  candidate_grid$candidate_id <- seq_len(nrow(candidate_grid))
  candidate_grid$selection_candidate <- TRUE

  # Benchmark is always its native/default model.  The fitting schedule below
  # treats lnorm/gamma/ig as direct-initialization models (Run 1) and all other
  # families as warm-start models (Run 2).
  ref_zi <- if ((is.null(Z) && isTRUE(check.zero)) || mixed_Z_input) c(FALSE,TRUE) else !is.null(Z)
  benchmark_grid <- data.frame(
    distribution=rep(benchmark.distribution,length(ref_zi)),
    fit.var.power=NA_integer_,
    var.power=rep(unname(default_vp[benchmark.distribution]),length(ref_zi)),
    zero_inflation=as.logical(ref_zi), stringsAsFactors=FALSE
  )

  direct_init <- c("lnorm","gamma","ig")
  warm_dist <- c("weibull","loglogistic","betaprime")

  if (!(benchmark.distribution %in% direct_init))
    stop("benchmark.distribution must be one of: lnorm, gamma, ig.")

  group1_requested <- candidate_grid[
    candidate_grid$distribution %in% direct_init, , drop=FALSE
  ]

  group2_grid <- candidate_grid[
    candidate_grid$distribution %in% warm_dist, , drop=FALSE
  ]

  # Run 2 follows the same initialization policy as select_zimp():
  # Weibull, Log-logistic, and Beta-prime always borrow beta/gamma
  # from a fitted Gamma model with the SAME effective variance power.
  #
  # If the corresponding Gamma model was not requested by the user,
  # create it as an internal Run-1 anchor.  Internal anchors never
  # participate in final model selection.
  anchor_extra <- candidate_grid[FALSE, , drop=FALSE]

  if (nrow(group2_grid)) {
    for (zi0 in unique(group2_grid$zero_inflation)) {

      needed_p <- sort(unique(
        as.integer(
          group2_grid$var.power[
            group2_grid$zero_inflation == zi0
          ]
        )
      ))

      for (pp in needed_p) {

        has_gamma_anchor <-
          any(
            group1_requested$distribution == "gamma" &
            group1_requested$zero_inflation == zi0 &
            as.integer(group1_requested$var.power) == pp
          ) ||
          any(
            benchmark_grid$distribution == "gamma" &
            benchmark_grid$zero_inflation == zi0 &
            as.integer(benchmark_grid$var.power) == pp
          )

        if (!has_gamma_anchor) {
          anchor_extra <- rbind(
            anchor_extra,
            data.frame(
              distribution = "gamma",
              fit.var.power = as.integer(pp),
              var.power = as.integer(pp),
              zero_inflation = as.logical(zi0),
              candidate_id = NA_integer_,
              selection_candidate = FALSE,
              stringsAsFactors = FALSE
            )
          )
        }
      }
    }
  }

  # Store a common matrix/data.frame only once.  A true list input remains a
  # list because it represents genuinely margin-specific designs.
  prepare_design_store <- function(A,name,allow_null=FALSE,columns=NULL) {
    if (is.null(A)) {
      if (allow_null) return(list(data=NULL,is_list=FALSE,columns=NULL))
      stop(name," must not be NULL.")
    }
    if (is.matrix(A)||is.data.frame(A)) {
      M <- design_matrix(A,name)
      cols <- columns
      if (!is.null(cols)) {
        if (!is.list(cols) || length(cols)!=k)
          stop(name,".columns must be NULL or a list of length k.")
        cols <- lapply(seq_len(k), function(j) {
          z <- cols[[j]]
          if (is.null(z)) return(seq_len(ncol(M)))
          if (is.character(z)) {
            ii <- match(z,colnames(M))
            if (anyNA(ii)) stop(name,".columns[[",j,"]] contains unknown column name(s).")
            return(as.integer(ii))
          }
          if (!is.numeric(z) || any(!is.finite(z)) || any(z!=as.integer(z)) ||
              any(z<1L) || any(z>ncol(M)))
            stop(name,".columns[[",j,"]] must contain valid column names or indices.")
          as.integer(z)
        })
      }
      return(list(data=M,is_list=FALSE,columns=cols))
    }
    if (is.list(A)) {
      if (!is.null(columns)) stop(name,".columns is only used when ",name," is a common matrix/data.frame.")
      if(length(A)!=k) stop(name," must be a matrix/data.frame or list of length k.")
      out <- vector("list",k)
      for(j in seq_len(k)) {
        if(is.null(A[[j]])) {
          if(allow_null) {out[j] <- list(NULL); next}
          stop(name,"[[",j,"]] is NULL.")
        }
        out[[j]] <- design_matrix(A[[j]],paste0(name,"[[",j,"]]"))
      }
      return(list(data=out,is_list=TRUE,columns=NULL))
    }
    stop(name," must be a matrix/data.frame or list of length k.")
  }

  X_store <- prepare_design_store(X,"X",columns=X.columns)
  Z_store <- prepare_design_store(Z,"Z",allow_null=TRUE,columns=Z.columns)

  get_design_j <- function(store,j) {
    if (is.null(store$data)) return(NULL)
    if (isTRUE(store$is_list)) return(store$data[[j]])
    if (is.null(store$columns)) return(store$data)
    store$data[,store$columns[[j]],drop=FALSE]
  }
  get_Xj <- function(j) get_design_j(X_store,j)
  get_Zj <- function(j) get_design_j(Z_store,j)
  margin_has_Z <- vapply(seq_len(k), function(j) !is.null(get_Zj(j)), logical(1))
  mixed_Z <- !is.null(Z) && isTRUE(Z_store$is_list)
  Z_intercept <- matrix(1,nrow=n,ncol=1L)

  get_user_ini <- function(dist,j) {
    if (is.null(ini.par)) return(NULL)
    if (!is.list(ini.par)) stop("ini.par must be NULL or a list.")
    if (!is.null(names(ini.par)) && dist %in% names(ini.par)) {
      z <- ini.par[[dist]]; if (is.null(z)) return(NULL)
      if (!is.list(z) || length(z)!=k) stop("ini.par[['",dist,"']] must be a list of length k.")
      return(z[[j]])
    }
    if (length(ini.par)==k) return(ini.par[[j]])
    stop("ini.par must be a list of length k or a named list by distribution.")
  }

  if (is.null(gh)) {
    if (!requireNamespace("statmod",quietly=TRUE)) stop("Package 'statmod' is required.")
    gh <- statmod::gauss.quad(n_gh,kind="hermite")
  }

  # parallel=NULL/1 => serial; integer >1 => PSOCK workers.
  if (is.null(parallel)) n_workers <- 1L else {
    if (length(parallel)!=1L || !is.numeric(parallel) || !is.finite(parallel) || parallel < 1 || parallel != as.integer(parallel))
      stop("`parallel` must be NULL or a positive integer.")
    n_workers <- as.integer(parallel)
  }
  run_parallel <- n_workers > 1L
  max_tasks <- max(k*nrow(benchmark_grid), k*nrow(candidate_grid), k)
  n_workers <- min(n_workers,max_tasks)

  cl <- NULL
  if (run_parallel) {
    cl <- parallel::makeCluster(n_workers,outfile="")
    on.exit(if (!is.null(cl)) try(parallel::stopCluster(cl),silent=TRUE),add=TRUE)
    dev_pkg_path <- if (file.exists("DESCRIPTION") && dir.exists("R")) normalizePath(".",winslash="/",mustWork=TRUE) else NULL
    parallel::clusterCall(cl,function(dev_pkg_path){
      # Prefer the active development tree when available; this avoids workers
      # silently using an older installed mzimp namespace.
      if (!is.null(dev_pkg_path)) {
        if (!requireNamespace("pkgload",quietly=TRUE)) stop("pkgload is required on development workers.")
        pkgload::load_all(dev_pkg_path,quiet=TRUE,recompile=FALSE,export_all=TRUE)
        return("development")
      }
      if (requireNamespace("mzimp",quietly=TRUE)) {
        suppressPackageStartupMessages(library(mzimp)); return("installed")
      }
      stop("mzimp is unavailable on workers.")
    },dev_pkg_path)
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl,iseed=as.integer(seed))
  }

  ctrl <- list(trace=trace,eval.max=eval.max,iter.max=iter.max,rel.tol=rel.tol,x.tol=x.tol,n_gh=n_gh,gh=gh)

  # Current helper names from fit_zimp are used here. Candidate fits deliberately
  # skip Hessian/inference; model selection and dependence fitting do not need it.
  fit_job <- function(job,Y,X_store,Z_store,Z_intercept,ctrl) {
    j <- job$margin; dist <- job$distribution
    zi <- isTRUE(job$zero_inflation)
    get_j <- function(store,j) {
      if (is.null(store$data)) return(NULL)
      if (isTRUE(store$is_list)) return(store$data[[j]])
      if (is.null(store$columns)) return(store$data)
      store$data[,store$columns[[j]],drop=FALSE]
    }
    Xj <- get_j(X_store,j)
    Zraw <- get_j(Z_store,j)
    Zj <- if (zi) { if (!is.null(Zraw)) Zraw else Z_intercept } else NULL
    pfit <- job$fit.var.power
    if (length(pfit)==0L || is.na(pfit)) pfit <- NULL else pfit <- as.integer(pfit)[1L]
    t0 <- proc.time()[["elapsed"]]
    # Each worker has already loaded the active package/DLL.  A scheduled
    # job represents exactly one marginal model, so call the public fit_zimp()
    # directly.  This is especially important in Run 2, where ini.par is the
    # target-specific warm start constructed from a same-var.power Run-1 anchor.
    ans <- tryCatch({
      # Each scheduled job is already exactly one model.  Use fit_zimp()
      # rather than select_zimp(): select_zimp() may create its own internal
      # anchor models and reinterpret ini.par, which is not appropriate for
      # Run-2 target-specific warm starts.
      fit_zimp(
        Y=Y[,j], X=Xj, Z=Zj, offset.count=NULL,
        mixing.distribution=dist, var.power=pfit,
        ini.par=job$start, progress=FALSE,
        trace=ctrl$trace, eval.max=ctrl$eval.max, iter.max=ctrl$iter.max,
        rel.tol=ctrl$rel.tol, x.tol=ctrl$x.tol,
        n_gh=ctrl$n_gh, gh=if(dist=="lnorm") ctrl$gh else NULL
      )
    },error=function(e) list(distribution=dist,var.power=job$var.power,
      neg_log_lik_value=Inf,log_lik_value=-Inf,be=NULL,phi=NA_real_,gam=NULL,
      par=NULL,convergence=999L,message=conditionMessage(e),hessian_ok=NA))

    # .fit_zimp_raw() deliberately catches low-level likelihood/Rcpp errors and
    # returns objective=1e100.  That value is finite, so is.finite() alone is
    # NOT a validity check.  Convert the sentinel to an explicit failed fit.
    if (!is.null(ans$neg_log_lik_value) &&
        is.finite(ans$neg_log_lik_value) && ans$neg_log_lik_value >= 1e90) {
      ans$log_lik_value <- -Inf
      ans$neg_log_lik_value <- Inf
      ans$convergence <- 999L
      ans$message <- paste0(
        "Likelihood returned numerical-error sentinel (>=1e90) for ", dist,
        ". The worker returned the numerical-error sentinel; inspect the underlying select_zimp()/likelihood setup."
      )
    }
    # Return only what model selection / warm-start / downstream fitting needs.
    # Do NOT return the full fit_zimp object (data, design matrices, Hessian, etc.).
    npar <- length(ans$par)
    ll <- ans$log_lik_value
    list(
      margin=j, distribution=dist, fit.var.power=pfit, var.power=job$var.power,
      use_custom_var_power=!is.null(pfit), zero_inflation=zi,
      be=ans$be, phi=ans$phi, gam=ans$gam, log_lik_value=ll,
      neg_log_lik_value=ans$neg_log_lik_value, npar=npar,
      AIC=if(is.finite(ll)) -2*ll+2*npar else Inf,
      BIC=if(is.finite(ll)) -2*ll+log(nrow(Y))*npar else Inf,
      convergence=ans$convergence, message=ans$message,
      elapsed_sec=unname(proc.time()[["elapsed"]]-t0)
    )
  }

  run_jobs <- function(jobs,label) {
    if (!length(jobs)) return(list())
    msg(label," | tasks = ",length(jobs)," | ",if(run_parallel) paste0("parallel | ",n_workers," cores") else "serial")
    if (!run_parallel) return(lapply(jobs,fit_job,Y=Y,X_store=X_store,Z_store=Z_store,Z_intercept=Z_intercept,ctrl=ctrl))
    parallel::clusterExport(cl,c("Y","X_store","Z_store","Z_intercept","ctrl","fit_job"),envir=environment())
    if (is.null(block_size)) {
      bs <- max(1L,ceiling(length(jobs)/(4L*n_workers)))
    } else bs <- max(1L,as.integer(block_size))
    blocks <- split(seq_along(jobs),ceiling(seq_along(jobs)/bs))
    out <- parallel::parLapplyLB(cl,blocks,function(ii,jobs,Y,X_store,Z_store,Z_intercept,ctrl,fit_job)
      lapply(jobs[ii],fit_job,Y=Y,X_store=X_store,Z_store=Z_store,Z_intercept=Z_intercept,ctrl=ctrl),
      jobs=jobs,Y=Y,X_store=X_store,Z_store=Z_store,Z_intercept=Z_intercept,ctrl=ctrl,fit_job=fit_job)
    unlist(out,recursive=FALSE)
  }

  # ---------------- Run 1: direct models + benchmark + Gamma anchors ----------------
  # If the user requests only lnorm/gamma/ig, this is the only fitting run.
  # If Weibull/Log-logistic/Beta-prime are requested, Run 1 also creates
  # any missing same-p Gamma anchors required to initialize Run 2.
  model_key <- function(dist, p, zi) paste(dist, as.integer(p), as.integer(zi), sep="|")

  run1_grid <- rbind(
    group1_requested[, c("distribution","fit.var.power","var.power","zero_inflation",
                         "candidate_id","selection_candidate"), drop=FALSE],
    transform(benchmark_grid, candidate_id=NA_integer_, selection_candidate=FALSE),
    anchor_extra[, c("distribution","fit.var.power","var.power","zero_inflation",
                     "candidate_id","selection_candidate"), drop=FALSE]
  )
  run1_grid$.key <- mapply(model_key, run1_grid$distribution, run1_grid$var.power, run1_grid$zero_inflation)
  # Prefer a requested candidate when duplicate rows represent the same model.
  run1_grid <- run1_grid[order(!run1_grid$selection_candidate), , drop=FALSE]
  run1_grid <- run1_grid[!duplicated(run1_grid$.key), , drop=FALSE]
  rownames(run1_grid) <- NULL

  run1_jobs <- list(); q <- 0L
  for (j in seq_len(k)) for (r in seq_len(nrow(run1_grid))) {
    rr <- run1_grid[r,]
    if (mixed_Z && isTRUE(rr$zero_inflation) != margin_has_Z[j]) next
    q <- q+1L
    pfit <- if (is.na(rr$fit.var.power)) NULL else as.integer(rr$fit.var.power)
    run1_jobs[[q]] <- list(margin=j, distribution=as.character(rr$distribution), fit.var.power=rr$fit.var.power, var.power=rr$var.power, zero_inflation=if(mixed_Z) margin_has_Z[j] else isTRUE(rr$zero_inflation), start=get_user_ini(as.character(rr$distribution),j), selection_candidate=isTRUE(rr$selection_candidate), run1_key=as.character(rr$.key))
  }

  estimation_t0 <- proc.time()[["elapsed"]]
  estimation_ram_start <- ram_snapshot_MB(cl)
  msg("Multivariate marginal selection: ",k," margins | benchmark = ",benchmark.distribution,
      " | candidates per margin = ",nrow(candidate_grid))
  run1_fits <- run_jobs(run1_jobs,"Run 1: direct models + benchmark/Gamma anchors")
  bad_run1 <- which(!vapply(run1_fits, function(x) is.finite(x$log_lik_value), logical(1)))
  if (length(bad_run1)) {
    ex <- unique(vapply(run1_fits[bad_run1], function(x)
      if (is.null(x$message)) "unknown worker error" else as.character(x$message)[1L], character(1)))
    stop("Run 1 fitting failed for ", length(bad_run1), "/", length(run1_fits),
         " tasks. First error: ", ex[1L])
  }

  # Restore scheduling metadata after worker fitting.
  for (ii in seq_along(run1_fits)) {
    run1_fits[[ii]]$selection_candidate <- isTRUE(run1_jobs[[ii]]$selection_candidate)
    run1_fits[[ii]]$benchmark_candidate <-
      identical(run1_fits[[ii]]$distribution, benchmark.distribution) &&
      as.integer(run1_fits[[ii]]$var.power) == as.integer(default_vp[benchmark.distribution])
    run1_fits[[ii]]$run1_key <- run1_jobs[[ii]]$run1_key
  }

  # Benchmark fit for reporting: exact native/default benchmark model.
  benchmark_fits <- run1_fits[vapply(run1_fits, function(x) isTRUE(x$benchmark_candidate), logical(1))]
  if (!length(benchmark_fits)) stop("Internal error: benchmark model is missing from Run 1.")

  # Retrieve the Gamma anchor with the SAME variance power for a Run-2 target.
  # This matches the current select_zimp() initialization policy exactly.
  get_anchor <- function(j, zi, p) {

    rr <- run1_fits[
      vapply(
        run1_fits,
        function(x)
          x$margin == j &&
          isTRUE(x$zero_inflation) == isTRUE(zi) &&
          identical(x$distribution, "gamma") &&
          as.integer(x$var.power) == as.integer(p) &&
          is.finite(x$log_lik_value),
        logical(1)
      )
    ]

    if (!length(rr)) return(NULL)

    # Duplicates should already have been removed in run1_grid.
    # Keep this deterministic if that assumption ever changes.
    rr[[1L]]
  }

  # ---------------- Run 2: Weibull / Log-logistic / Beta-prime ----------------
  # Each target borrows beta/gamma from Gamma with the SAME var.power.
  # Target phi is still initialized from the target distribution itself.
  run2_jobs <- list(); q <- 0L
  if (nrow(group2_grid)) {
    init_fun <- if (exists(".zimp_init_par", mode="function", inherits=TRUE))
      get(".zimp_init_par", mode="function", inherits=TRUE) else
      getFromNamespace(".zimp_init_par", "mzimp")

    for (j in seq_len(k)) for (r in seq_len(nrow(group2_grid))) {
      rr <- group2_grid[r,]
      if (mixed_Z && isTRUE(rr$zero_inflation) != margin_has_Z[j]) next
      dist <- as.character(rr$distribution); zi <- if(mixed_Z) margin_has_Z[j] else isTRUE(rr$zero_inflation)
      pfit <- if (is.na(rr$fit.var.power)) NULL else as.integer(rr$fit.var.power)
      start <- get_user_ini(dist,j)

      if (is.null(start)) {
        a <- get_anchor(j, zi, rr$var.power)
        if (is.null(a)) stop("Internal error: no Gamma Run-1 anchor for ",dist,
                             " with var.power = ",rr$var.power," on margin ",j,".")
        Zj0 <- get_Zj(j); Zj <- if(zi){if(!is.null(Zj0)) Zj0 else Z_intercept} else NULL
        base <- init_fun(y=Y[,j], X=get_Xj(j), Z=Zj, offset.count=NULL,
                         dist=dist, zero_inflation=zi, var.power=pfit)
        base$be <- a$be
        if (zi && !is.null(a$gam)) base$gam <- a$gam
        base$par <- c(base$be, log(base$phi), base$gam)
        start <- base$par
      }

      q <- q+1L
      run2_jobs[[q]] <- list(
        margin=j, distribution=dist, fit.var.power=rr$fit.var.power,
        var.power=rr$var.power, zero_inflation=zi, start=start,
        selection_candidate=TRUE
      )
    }
  }

  run2_fits <- if (length(run2_jobs))
    run_jobs(run2_jobs,"Run 2: Weibull/Log-logistic/Beta-prime from Gamma anchors") else list()

  bad_run2 <- which(!vapply(run2_fits, function(x) is.finite(x$log_lik_value), logical(1)))
  if (length(bad_run2)) {
    ex <- unique(vapply(run2_fits[bad_run2], function(x)
      if (is.null(x$message)) "unknown worker error" else as.character(x$message)[1L], character(1)))
    warning("Run 2: ", length(bad_run2), "/", length(run2_fits),
            " candidate fits failed. First error: ", ex[1L])
  }
  if (length(run2_fits)) {
    run2_fits <- lapply(run2_fits, function(z) {
      z$selection_candidate <- TRUE
      z$benchmark_candidate <- FALSE
      z
    })
  }

  # all_results contains Run-1 requested candidates, Run-2 candidates, and the
  # additional benchmark/anchor fits.  Only selection_candidate=TRUE participates
  # in model selection.
  estimation_elapsed <- unname(proc.time()[["elapsed"]] - estimation_t0)
  estimation_ram_end <- ram_snapshot_MB(cl)

  fits <- c(run1_fits, run2_fits)

  # Select best candidate and benchmark per margin.
  best_by_margin <- vector("list",k); benchmark_by_margin <- vector("list",k)
  for (j in seq_len(k)) {
    rr <- fits[vapply(fits,function(x)x$margin==j && isTRUE(x$selection_candidate),logical(1))]
    ic <- vapply(rr,function(x)x[[selection.criterion]],numeric(1))
    if (length(rr) && any(is.finite(ic))) best_by_margin[[j]] <- rr[[which.min(ic)]]
    rb <- benchmark_fits[vapply(benchmark_fits,function(x)x$margin==j,logical(1))]
    ib <- vapply(rb,function(x)x[[selection.criterion]],numeric(1))
    if (length(rb) && any(is.finite(ib))) benchmark_by_margin[[j]] <- rb[[which.min(ib)]]
  }
  if (any(vapply(best_by_margin,is.null,logical(1)))) warning("At least one margin has no valid selected model.")

  # Compact selected-margin objects for downstream fit_multizimp().
  # Keep only fields needed to reuse the fitted marginal model and, when
  # necessary, recompute the interval-probability matrices I.
  selected_margins <- vector("list",k)
  for (j in seq_len(k)) {
    b <- best_by_margin[[j]]; if(is.null(b)) next
    selected_margins[[j]] <- list(
      margin=j,
      zero_inflation=isTRUE(b$zero_inflation),
      distribution=b$distribution,
      var.power=b$var.power,
      var_power=b$var.power,
      use_custom_var_power=if(!is.null(b$use_custom_var_power))
        isTRUE(b$use_custom_var_power) else TRUE,
      be=b$be,
      phi=b$phi,
      gam=b$gam,
      log_lik_value=b$log_lik_value
    )
    class(selected_margins[[j]]) <- c("mzimp_margin_fit",class(selected_margins[[j]]))
  }

  best_ll <- vapply(best_by_margin,function(x)if(is.null(x)) NA_real_ else x$log_lik_value,numeric(1))
  bench_ll <- vapply(benchmark_by_margin,function(x)if(is.null(x)) NA_real_ else x$log_lik_value,numeric(1))
  gain <- best_ll-bench_ll; ok <- is.finite(gain)
  benchmark_summary <- list(distribution=benchmark.distribution,per_margin_gain=gain,
    per_margin_best_loglik=best_ll,per_margin_benchmark_loglik=bench_ll,
    total_best_loglik=sum(best_ll[ok]),total_benchmark_loglik=sum(bench_ll[ok]),
    total_loglik_gain=sum(gain[ok]),mean_loglik_gain=if(any(ok))mean(gain[ok]) else NA_real_,
    n_compared=sum(ok),n_total=k)

  selected_summary <- data.frame(
    margin=seq_len(k),
    distribution=vapply(best_by_margin,function(x)if(is.null(x))NA_character_ else x$distribution,character(1)),
    var.power=vapply(best_by_margin,function(x)if(is.null(x))NA_real_ else as.numeric(x$var.power),numeric(1)),
    zero_inflation=vapply(best_by_margin,function(x)!is.null(x)&&isTRUE(x$zero_inflation),logical(1)),
    logLik=best_ll,AIC=vapply(best_by_margin,function(x)if(is.null(x))NA_real_ else x$AIC,numeric(1)),
    BIC=vapply(best_by_margin,function(x)if(is.null(x))NA_real_ else x$BIC,numeric(1)),
    benchmark.logLik=bench_ll,gain.vs.benchmark=gain,stringsAsFactors=FALSE)

  # Compact benchmark records before releasing full candidate fits.
  benchmark_by_margin_compact <- lapply(seq_len(k), function(j) {
    b <- benchmark_by_margin[[j]]
    if (is.null(b)) return(NULL)
    list(margin=j, distribution=b$distribution, var.power=b$var.power,
         zero_inflation=isTRUE(b$zero_inflation), logLik=b$log_lik_value,
         AIC=b$AIC, BIC=b$BIC)
  })

  computeI_t0 <- if (isTRUE(compute.I)) proc.time()[["elapsed"]] else NA_real_
  computeI_ram_start <- if (isTRUE(compute.I)) ram_snapshot_MB(cl) else NULL

  # The estimation workers have finished their job.  Stop them BEFORE compute-I
  # so their heaps, candidate fits, closures and temporary allocations disappear.
  if (run_parallel && !is.null(cl)) {
    parallel::stopCluster(cl)
    cl <- NULL
  }
  rm(fits, run1_fits, run2_fits, benchmark_fits, best_by_margin,
     benchmark_by_margin, run1_jobs, run2_jobs)
  invisible(gc())

  # Optional I cache.  A fresh PSOCK cluster is used.  Each task contains only
  # one response margin, its X/Z design, selected parameters, and I controls.
  I_y_list <- I_0_list <- pb_list <- y0_list <- NULL
  if (isTRUE(compute.I)) {
    if (any(vapply(selected_margins,is.null,logical(1))))
      stop("Cannot compute I: at least one selected margin failed.")
    msg("Computing I matrices for ",k," selected margins | d0 = ",I.d0)

    I_tasks <- lapply(seq_len(k), function(j) list(
      y=unname(Y[,j]),
      X=unname(as.matrix(get_Xj(j))),
      Z={zj <- get_Zj(j); if (is.null(zj)) NULL else unname(as.matrix(zj))},
      margin_fit=selected_margins[[j]],
      I.d0=I.d0,
      n_gh=n_gh,
      gh=gh
    ))

    if (run_parallel) {
      cl <- parallel::makeCluster(n_workers,outfile="")
      dev_pkg_path <- if (file.exists("DESCRIPTION") && dir.exists("R")) normalizePath(".",winslash="/",mustWork=TRUE) else NULL
      parallel::clusterCall(cl,function(dev_pkg_path){
        if (!is.null(dev_pkg_path)) {
          if (!requireNamespace("pkgload",quietly=TRUE)) stop("pkgload is required on development workers.")
          pkgload::load_all(dev_pkg_path,quiet=TRUE,recompile=FALSE,export_all=TRUE)
          return("development")
        }
        if (requireNamespace("mzimp",quietly=TRUE)) {
          suppressPackageStartupMessages(library(mzimp)); return("installed")
        }
        stop("mzimp is unavailable on workers.")
      },dev_pkg_path)
      if (!is.null(seed)) parallel::clusterSetRNGStream(cl,iseed=as.integer(seed))
      # Stand-alone helper: no enclosing select_multi_margin environment is captured.
      worker_fun <- .mzimp_compute_I_margin_worker
      Ires <- parallel::parLapplyLB(cl,I_tasks,worker_fun)
      parallel::stopCluster(cl); cl <- NULL
    } else {
      Ires <- lapply(I_tasks,.mzimp_compute_I_margin_worker)
    }
    rm(I_tasks); invisible(gc())
    I_y_list<-lapply(Ires,`[[`,"I_y")
    I_0_list<-lapply(Ires,`[[`,"I_0")
    pb_list<-lapply(Ires,`[[`,"pb")
    y0_list<-lapply(Ires,`[[`,"y0")
    rm(Ires); invisible(gc())
  }

  computeI_elapsed <- if (isTRUE(compute.I)) unname(proc.time()[["elapsed"]] - computeI_t0) else NA_real_
  computeI_ram_end <- if (isTRUE(compute.I)) ram_snapshot_MB(NULL) else NULL

  performance <- list(
    estimation=list(elapsed_sec=estimation_elapsed, ram_start_MB=estimation_ram_start$total,
                    ram_end_MB=estimation_ram_end$total, max_worker_end_MB=estimation_ram_end$max_worker),
    compute_I=if(isTRUE(compute.I)) list(elapsed_sec=computeI_elapsed,
                    ram_start_MB=computeI_ram_start$total, ram_end_MB=computeI_ram_end$total) else NULL
  )
  if (isTRUE(progress)) {
    cat(sprintf("\n[1] Marginal estimation | elapsed = %.2f s (%.2f min) | RAM start = %.1f MB | RAM end = %.1f MB | max worker end = %.1f MB\n",
                estimation_elapsed,estimation_elapsed/60,estimation_ram_start$total,
                estimation_ram_end$total,estimation_ram_end$max_worker))
    if (isTRUE(compute.I))
      cat(sprintf("[2] Compute I | elapsed = %.2f s (%.2f min) | RAM before cleanup = %.1f MB | RAM after = %.1f MB\n",
                  computeI_elapsed,computeI_elapsed/60,computeI_ram_start$total,computeI_ram_end$total))
  }

  summary_obj <- list(
    selected_summary=selected_summary,
    benchmark=benchmark_summary,
    benchmark.distribution=benchmark.distribution,
    selection.criterion=selection.criterion,
    n=n,k=k,response_names=colnames(Y),
    performance=performance
  )
  out <- list(Y=Y,
    X=X_store$data,Z=Z_store$data,
    X_is_list=isTRUE(X_store$is_list),Z_is_list=isTRUE(Z_store$is_list),
    X_columns=X_store$columns,Z_columns=Z_store$columns,
    selected_margins=selected_margins,
    benchmark.distribution=benchmark.distribution,
    benchmark_by_margin=benchmark_by_margin_compact,
    benchmark=benchmark_summary,selected_summary=selected_summary,
    I_precomputed=isTRUE(compute.I),I_d0=if(compute.I)I.d0 else NULL,
    I_y_list=I_y_list,I_0_list=I_0_list,pb_list=pb_list,y0_list=y0_list,
    mixing.distribution=mixing.distribution,var.power=var.power,candidate_grid=candidate_grid,
    check.zero=isTRUE(check.zero),selection.criterion=selection.criterion,
    n=n,k=k,response_names=colnames(Y),parallel=parallel,n_workers=n_workers,performance=performance,summary=summary_obj)
  class(out) <- "multi_margin_fit"
  out
}
