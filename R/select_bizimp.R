
# =============================================================================
# select_bizimp.R
# Bivariate ZIMP model selection
#
# Requires:
#   - fit_zimp(), select_zimp()
#   - helpers / dependence engine from the current fit_bizimp() implementation:
#       .bizimp_resolve_input()
#       .fit_bizimp_dependence_engine()
#       .bizimp_public_copula_name()
#
# Public API:
#   select_bizimp()
# =============================================================================


.bizimp_selection_zero_label <- function(option, theta.B = NA_real_) {
  if (is.na(option) || is.null(option)) {
    return(list(name = "Not active", parameter = NA_character_))
  }

  if (identical(option, "normal")) {
    return(list(
      name = "Gaussian copula",
      parameter = paste0("rho = ", format(signif(theta.B, 7), trim = TRUE))
    ))
  }

  if (identical(option, "comonotone")) {
    return(list(
      name = "Frechet upper",
      parameter = ""
    ))
  }

  if (identical(option, "countermonotone")) {
    return(list(
      name = "Frechet lower",
      parameter = ""
    ))
  }

  list(name = as.character(option), parameter = NA_character_)
}


.bizimp_selection_copula_parameter <- function(model) {
  fam <- as.integer(model$family)

  if (fam == 0L) return("none")

  if (fam == 2L) {
    return(paste0(
      "rho = ", format(signif(model$par, 7), trim = TRUE)
    ))
  }

  th <- model$theta.G
  if (is.null(th)) {
    th <- if (!is.null(model$par2) && is.finite(model$par2) &&
              abs(model$par2) > 0) {
      c(model$par, model$par2)
    } else {
      model$par
    }
  }

  th <- as.numeric(th)
  if (length(th) == 1L) {
    paste0("parameter = ", format(signif(th, 7), trim = TRUE))
  } else {
    paste0(
      "parameters = ",
      paste(format(signif(th, 7), trim = TRUE), collapse = ", ")
    )
  }
}


.bizimp_selection_margin_label <- function(fit) {
  vp <- fit$var.power
  vp_txt <- if (is.null(vp) || length(vp) == 0L) {
    "native"
  } else {
    as.character(vp[1L])
  }

  paste0(
    fit$mixing.distribution,
    " | var.power = ", vp_txt,
    " | logLik = ", format(signif(fit$log_lik_value, 8), trim = TRUE)
  )
}


.bizimp_selection_num3 <- function(x) {
  if (is.null(x) || !length(x) || !is.finite(x[1L])) return(NA_real_)
  round(as.numeric(x[1L]), 3L)
}


.bizimp_selection_display_var_power <- function(fit) {
  vp <- fit$var.power
  if (!is.null(vp) && length(vp) && is.finite(vp[1L])) {
    return(as.numeric(vp[1L]))
  }

  # Native/default variance powers used by the marginal parameterizations.
  defaults <- c(
    gamma = 2, lnorm = 2, pareto = 2, ig = 3,
    weibull = 2, loglogistic = 2, betaprime = 2
  )
  dist <- as.character(fit$mixing.distribution)[1L]
  if (!is.na(dist) && dist %in% names(defaults)) {
    return(unname(defaults[[dist]]))
  }

  NA_real_
}


.bizimp_selection_zero_short <- function(option, theta.B = NA_real_) {
  if (is.na(option) || is.null(option)) {
    return(list(name = "Not active", parameter = NA_real_))
  }
  if (identical(option, "normal")) {
    return(list(name = "Gaussian",
                parameter = .bizimp_selection_num3(theta.B)))
  }
  if (identical(option, "comonotone")) {
    return(list(name = "Frechet upper",
                parameter = NA_real_))
  }
  if (identical(option, "countermonotone")) {
    return(list(name = "Frechet lower",
                parameter = NA_real_))
  }
  list(name = as.character(option),
       parameter = .bizimp_selection_num3(theta.B))
}


.bizimp_selection_mixing_name <- function(model) {
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


.bizimp_selection_mixing_para3 <- function(model) {
  fam <- as.integer(model$family)
  if (fam == 0L) return(NA_character_)

  if (fam == 2L) {
    return(format(round(as.numeric(model$par), 3L), nsmall = 0, trim = TRUE))
  }

  th <- model$theta.G
  if (is.null(th)) {
    th <- if (!is.null(model$par2) && length(model$par2) &&
              is.finite(model$par2) && abs(model$par2) > 0) {
      c(model$par, model$par2)
    } else {
      model$par
    }
  }

  paste(
    format(round(as.numeric(th), 3L), nsmall = 0, trim = TRUE),
    collapse = ", "
  )
}



.bizimp_selection_margin_compact <- function(fit) {
  dist <- as.character(fit$mixing.distribution)[1L]
  vp <- fit$var.power

  # NULL means the native/default parameterization: show only distribution.
  if (is.null(vp) || !length(vp) || !is.finite(vp[1L])) {
    return(dist)
  }

  paste0(dist, " (v.p = ",
         format(as.numeric(vp[1L]), trim = TRUE), ")")
}

.bizimp_selection_zero_compact <- function(option, theta.B = NA_real_) {
  if (is.na(option) || is.null(option)) return("Not active")

  if (identical(option, "normal")) {
    return(paste0(
      "Gaussian (rho = ",
      format(round(as.numeric(theta.B), 3L), nsmall = 0, trim = TRUE),
      ")"
    ))
  }
  if (identical(option, "comonotone")) return("Frechet upper")
  if (identical(option, "countermonotone")) return("Frechet lower")

  as.character(option)
}

.bizimp_selection_mixing_compact <- function(model) {
  fam <- as.integer(model$family)
  nm <- .bizimp_selection_mixing_name(model)

  if (fam == 0L) return("Independence")

  if (fam == 1L) {
    return(paste0(
      nm, " (rho = ",
      format(round(as.numeric(model$par), 3L), nsmall = 0, trim = TRUE),
      ")"
    ))
  }

  if (fam == 2L) {
    df <- suppressWarnings(as.numeric(model$par2)[1L])
    return(paste0(
      "Student-t (df = ",
      format(round(df, 3L), nsmall = 0, trim = TRUE),
      ", rho = ",
      format(round(as.numeric(model$par), 3L), nsmall = 0, trim = TRUE),
      ")"
    ))
  }

  th <- model$theta.G
  if (is.null(th)) {
    th <- if (!is.null(model$par2) && length(model$par2) &&
              is.finite(model$par2) && abs(model$par2) > 0) {
      c(model$par, model$par2)
    } else {
      model$par
    }
  }
  th <- as.numeric(th)

  if (length(th) == 1L) {
    paste0(
      nm, " (theta = ",
      format(round(th, 3L), nsmall = 0, trim = TRUE),
      ")"
    )
  } else {
    paste0(
      nm, " (theta = ",
      paste(format(round(th, 3L), nsmall = 0, trim = TRUE),
            collapse = ", "),
      ")"
    )
  }
}


.bizimp_selection_print_margin_table <- function(sel, margin_no) {
  tab <- sel$summary_table

  # Keep the marginal selection table compact.
  drop_cols <- intersect(
    names(tab),
    c("negative_log_likelihood", "neg_log_lik", "neg.logLik", "NegLogLik")
  )
  if (length(drop_cols)) {
    tab <- tab[, setdiff(names(tab), drop_cols), drop = FALSE]
  }

  cat("\nMargin ", margin_no, " selection completed\n", sep = "")
  cat(strrep("-", 72), "\n", sep = "")
  print(tab, row.names = FALSE)
  invisible(NULL)
}


.bizimp_selection_print_pair <- function(pair_id, total_pairs, i, j, f1, f2) {
  cat(
    "\n[Copula pair ", pair_id, "/", total_pairs, "] ",
    "margin 1 rank = ", i,
    " | margin 2 rank = ", j, "\n",
    sep = ""
  )
  cat("  Margin 1: ", .bizimp_selection_margin_label(f1), "\n", sep = "")
  cat("  Margin 2: ", .bizimp_selection_margin_label(f2), "\n", sep = "")
  flush.console()
}






#' Select a Bivariate Mixed Poisson or Zero-Inflated Mixed Poisson Model
#'
#' @description
#' Performs model selection for bivariate Mixed Poisson (MP) and
#' Zero-Inflated Mixed Poisson (ZIMP) regression models.
#'
#' Marginal models are first selected independently using [select_zimp()].
#' The best `margin.top` models from each margin are retained, and dependence
#' structures are then evaluated for every retained marginal-model pair.
#' The final bivariate model is selected by the largest maximized
#' observed-data log-likelihood.
#'
#' @details
#' `select_bizimp()` extends [select_zimp()] from univariate marginal
#' selection to bivariate model selection. In contrast to [fit_bizimp()],
#' which fits one specified bivariate model, `select_bizimp()` searches over
#' marginal specifications and mixing-copula families.
#'
#' The selection procedure consists of three stages:
#'
#' \enumerate{
#'   \item Select candidate models independently for each margin using
#'     [select_zimp()].
#'   \item Retain the best `margin.top` models from each margin and perform
#'     a Gaussian dependence prescreen for every retained marginal pair.
#'   \item Fit the compatible user-requested mixing copulas, together with
#'     the admissible zero-dependence structures, and select the candidate
#'     with the largest maximized bivariate log-likelihood.
#' }
#'
#' @section Model specification:
#'
#' Two input interfaces are supported.
#'
#' **Formula interface**
#'
#' Supply `formula1`, `formula2`, and `data`. Each formula follows the
#' syntax used by [select_zimp()]. For example,
#'
#' ```
#' formula1 = y1 ~ 1 + x1 + x2 | z1 + z2
#' formula2 = y2 ~ 1 + x1 + x2
#' ```
#'
#' specifies a ZIMP first margin and an MP second margin.
#'
#' The expression before `|` specifies the count component and the
#' expression after `|` specifies the zero-inflation component.
#' Intercepts are included only when explicitly requested with `1`.
#'
#' Both formulas must refer to the same observations. If missing values
#' cause different observations to be retained for the two margins,
#' bivariate model selection is stopped.
#'
#' **Matrix interface**
#'
#' Alternatively, supply `Y`, `X`, and optionally `Z`.
#'
#' `Y` may be an \eqn{n \times 2} matrix/data frame or a list containing
#' two count vectors. `X` may be a design matrix shared by both margins or
#' a list containing two margin-specific design matrices.
#'
#' `Z = NULL` specifies two MP margins. A shared matrix/data frame `Z`
#' specifies two ZIMP margins. A list may be used for mixed configurations;
#' for example, `Z = list(Z1, NULL)` specifies ZIMP--MP and
#' `Z = list(NULL, Z2)` specifies MP--ZIMP.
#'
#' @section Marginal model selection:
#'
#' Each margin is selected independently using [select_zimp()] with the
#' candidate mixing distributions supplied in `mixing.distribution` and
#' the variance-power specifications supplied in `var.power`.
#'
#' Supported mixing distributions are:
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
#' If `var.power = NULL`, each distribution uses its native/default
#' variance-power specification. Otherwise, every requested mixing
#' distribution is crossed with the requested variance powers, following
#' the candidate construction used by [select_zimp()].
#'
#' Marginal candidates are ranked by decreasing maximized univariate
#' log-likelihood. At most `margin.top` models are retained from each margin.
#' If \eqn{m_1} and \eqn{m_2} models are retained from margins 1 and 2,
#' respectively, dependence selection is performed for
#'
#' \deqn{m_1 m_2}
#'
#' marginal-model pairs.
#'
#' Thus, `margin.top = 1` restricts the dependence search to the independently
#' best marginal model from each margin.
#'
#' @section Mixing-copula candidates:
#'
#' Dependence between the two latent mixing variables is specified using
#' bivariate `VineCopula` family codes. `mixing.copula` may contain one or
#' more candidate families.
#'
#' Supported family codes are:
#'
#' \describe{
#'   \item{`0L`}{Independence copula.}
#'   \item{`1L`}{Gaussian copula.}
#'   \item{`2L`}{Student-t copula.}
#'   \item{`3L`}{Clayton copula.}
#'   \item{`4L`}{Gumbel copula.}
#'   \item{`5L`}{Frank copula.}
#'   \item{`6L`}{Joe copula.}
#'   \item{`7L`}{BB1 copula.}
#'   \item{`8L`}{BB6 copula.}
#'   \item{`9L`}{BB7 copula.}
#'   \item{`10L`}{BB8 copula.}
#'   \item{`13L`, `14L`, `16L`--`20L`}{
#'     180-degree rotations of the corresponding Archimedean and BB
#'     families.
#'   }
#'   \item{`23L`, `24L`, `26L`--`30L`}{
#'     90-degree rotations of the corresponding Archimedean and BB
#'     families.
#'   }
#'   \item{`33L`, `34L`, `36L`--`40L`}{
#'     270-degree rotations of the corresponding Archimedean and BB
#'     families.
#'   }
#'   \item{`104L`, `114L`, `124L`, `134L`,
#'         `204L`, `214L`, `224L`, `234L`}{
#'     Supported Tawn families and their rotations.
#'   }
#' }
#'
#' The default is `mixing.copula = 1L`, corresponding to the Gaussian
#' copula.
#'
#' For example,
#'
#' ```
#' mixing.copula = 1L             # Gaussian only
#' mixing.copula = c(1L, 3L, 5L) # Gaussian, Clayton, and Frank
#' ```
#'
#' For the Student-t copula (`2L`), the correlation parameter is optimized
#' for each degrees-of-freedom value supplied in `t.df.grid`, and the
#' likelihood-maximizing specification is retained.
#'
#' @section Gaussian dependence prescreen:
#'
#' For every retained marginal pair, a Gaussian checkerboard copula is
#' fitted first to the latent mixing variables. This fit is used only as a
#' prescreen for the sign of the mixing dependence.
#'
#' If the absolute fitted Gaussian correlation exceeds
#' `gaussian.sign.tol`, user-requested copula families incompatible with
#' that dependence direction are removed before the full dependence search.
#' If the fitted correlation is not finite or does not exceed the tolerance,
#' no sign restriction is imposed.
#'
#' If sign filtering would remove all user-requested families, all requested
#' families are restored. Consequently, the prescreen never introduces a
#' mixing-copula family that was not requested by the user.
#'
#' The Gaussian family is fitted internally for this prescreen even when
#' `1L` is not included in `mixing.copula`, but it is eligible for final
#' selection only when explicitly requested.
#'
#' @section Zero dependence:
#'
#' When both retained margins are ZIMP models, dependence between their
#' Bernoulli activation indicators is also selected. For every compatible
#' mixing copula, the following three zero-dependence specifications are
#' evaluated:
#'
#' \describe{
#'   \item{Gaussian}{
#'     Gaussian copula dependence with an estimated correlation parameter.
#'   }
#'   \item{Frechet upper}{
#'     Comonotone dependence at the upper Frechet bound.
#'   }
#'   \item{Frechet lower}{
#'     Countermonotone dependence at the lower Frechet bound.
#'   }
#' }
#'
#' Thus, for a ZIMP--ZIMP marginal pair, the full dependence search crosses
#' each compatible mixing-copula family with all three zero-dependence
#' specifications.
#'
#' If either margin is an MP model, the zero-dependence layer is inactive.
#'
#' @section Final model selection:
#'
#' The complete candidate set consists of the fitted combinations of
#'
#' \deqn{
#' \text{retained margin 1}
#' \times
#' \text{retained margin 2}
#' \times
#' \text{compatible mixing copula}
#' \times
#' \text{admissible zero dependence}.
#' }
#'
#' Candidates are evaluated using the full bivariate observed-data
#' likelihood. `candidate_table` is sorted by decreasing maximized
#' log-likelihood, and its first row corresponds to the selected model.
#'
#' AIC and BIC are also calculated for each candidate. Their parameter
#' counts include both marginal models, the mixing-copula parameters, and
#' the Gaussian zero-copula correlation parameter when applicable.
#' AIC and BIC are reported for comparison but are not used for model
#' selection.
#'
#' @section Parallel computation:
#'
#' `parallel` is passed to [select_zimp()] and controls parallel fitting of
#' candidate marginal models within each margin.
#'
#' The two marginal selections themselves are performed sequentially.
#' Dependence searches across retained marginal pairs, mixing copulas, and
#' zero-dependence specifications are also performed sequentially by
#' `select_bizimp()`.
#'
#' @param formula1 Optional formula specifying the first marginal model.
#'   The syntax follows [select_zimp()]. Used together with `formula2` and
#'   `data`.
#'
#' @param formula2 Optional formula specifying the second marginal model.
#'   The syntax follows [select_zimp()]. Used together with `formula1` and
#'   `data`.
#'
#' @param data Optional data frame containing the variables referenced by
#'   `formula1` and `formula2`. Required when the formula interface is used.
#'
#' @param Y Response for the direct matrix interface. Either an
#'   \eqn{n \times 2} numeric matrix/data frame or a list containing two
#'   numeric vectors of finite, non-negative integer counts.
#'
#' @param X Count-component design for the direct matrix interface. May be
#'   one numeric matrix/data frame shared by both margins or a list
#'   containing two margin-specific design matrices. Include intercept
#'   columns explicitly when required.
#'
#' @param Z Optional zero-component design for the direct matrix interface.
#'   May be `NULL`, one numeric matrix/data frame shared by both margins,
#'   or a list containing two margin-specific designs. A `NULL` list
#'   element specifies an MP margin; a non-`NULL` element specifies a
#'   ZIMP margin.
#'
#' @param mixing.distribution Character vector specifying candidate marginal
#'   mixing distributions passed to [select_zimp()]. Supported values are
#'   `"gamma"`, `"lnorm"`, `"ig"`, `"weibull"`, `"loglogistic"`, and
#'   `"betaprime"`. By default, all supported distributions are considered.
#'
#' @param var.power Optional numeric vector specifying candidate
#'   variance-power parameterizations passed to [select_zimp()]. Values
#'   must be selected from `1`, `2`, and `3`. If `NULL`, each requested
#'   mixing distribution uses its native/default parameterization.
#'
#' @param margin.top Positive integer specifying the maximum number of
#'   independently ranked marginal models retained from each margin.
#'   The default is `1L`. If `m` models are retained from each margin, up
#'   to \eqn{m^2} marginal pairs are evaluated.
#'
#' @param mixing.copula Integer vector containing one or more candidate
#'   `VineCopula` family codes for dependence between the latent mixing
#'   variables. See the "Mixing-copula candidates" section for supported
#'   codes. The default is `1L` (Gaussian).
#'
#' @param d0 Positive integer specifying the checkerboard grid size used
#'   to discretize the mixing copulas. The default is `10L`.
#'
#' @param n_gh Positive integer specifying the number of Gauss--Hermite
#'   quadrature nodes used when required by the marginal calculations.
#'   The default is `20`.
#'
#' @param gh Optional precomputed Gauss--Hermite quadrature object passed
#'   to the marginal-selection and dependence calculations. If `NULL`,
#'   the required quadrature rule is constructed internally.
#'
#' @param ini.par `NULL` or a list of length two containing optional initial
#'   parameter vectors for the two marginal selections. The first element
#'   is passed to margin 1 and the second to margin 2.
#'
#' @param parallel `NULL` or a positive integer specifying the number of
#'   workers passed to [select_zimp()] for candidate-model fitting within
#'   each margin. `NULL` or `1` uses serial marginal selection.
#'
#' @param progress Logical compatibility argument retained for the public
#'   interface. `select_bizimp()` performs model selection silently.
#'
#' @param show.result Logical compatibility argument retained for the public
#'   interface. Results are stored in the returned object rather than
#'   printed during fitting.
#'
#' @param trace Numeric optimizer trace argument retained for compatibility.
#'   The current selection interface performs fitting silently.
#'
#' @param eval.max Positive integer giving the maximum number of
#'   objective-function evaluations used by the marginal and dependence
#'   optimizers. The default is `100`.
#'
#' @param rel.tol Positive numeric value giving the relative convergence
#'   tolerance passed to the marginal model-selection routines. The default
#'   is `1e-5`.
#'
#' @param x.tol Positive numeric value giving the parameter convergence
#'   tolerance passed to the marginal model-selection routines. The default
#'   is `1e-5`.
#'
#' @param iter.max Positive integer giving the maximum number of optimization
#'   iterations. The default is `100`.
#'
#' @param gaussian.sign.tol Non-negative numeric tolerance used to determine
#'   whether the Gaussian prescreen provides an informative dependence
#'   direction. The default is `1e-6`.
#'
#' @param t.df.grid Integer vector containing candidate degrees of freedom
#'   greater than two for Student-t mixing copulas. The default is
#'   `c(3, 5, 10, 30)`.
#'
#' @return
#' An object of class `"bizimp_selection"`. Important components include:
#'
#' \describe{
#'   \item{`margin_selection_1`, `margin_selection_2`}{
#'     Results from the independent [select_zimp()] selections for margins
#'     1 and 2.
#'   }
#'   \item{`margin_top`}{
#'     Requested maximum number of marginal models retained from each margin.
#'   }
#'   \item{`retained_margin_1`, `retained_margin_2`}{
#'     Lists of fitted marginal models retained for the dependence search.
#'   }
#'   \item{`pair_results`}{
#'     Detailed dependence-search results for every retained marginal pair,
#'     including the Gaussian prescreen, compatible mixing-copula families,
#'     zero-dependence candidates, and fitted dependence models.
#'   }
#'   \item{`candidate_table`}{
#'     Data frame containing all fitted bivariate candidates, with columns
#'     describing margin 1, margin 2, zero dependence, mixing dependence,
#'     maximized log-likelihood, AIC, and BIC. Rows are ordered by decreasing
#'     log-likelihood.
#'   }
#'   \item{`pair_summary`}{
#'     Data frame containing the best dependence model obtained for each
#'     retained marginal pair, ordered by decreasing bivariate log-likelihood.
#'   }
#'   \item{`best_pair`}{
#'     Detailed result for the retained marginal pair containing the globally
#'     selected bivariate model.
#'   }
#'   \item{`best_model`}{
#'     Selected dependence model. It contains the selected mixing-copula
#'     family and parameter estimates, zero-dependence parameter when
#'     applicable, maximized bivariate log-likelihood, convergence
#'     information, total parameter count, AIC, and BIC.
#'   }
#'   \item{`best_dependence`}{
#'     Complete dependence-engine result associated with `best_model`.
#'   }
#'   \item{`best_margin_1`, `best_margin_2`}{
#'     Fitted marginal models belonging to the globally selected bivariate
#'     model.
#'   }
#'   \item{`interface`}{
#'     Character string identifying whether the formula or direct matrix
#'     interface was used.
#'   }
#'   \item{`formula1`, `formula2`}{
#'     Marginal formulas when the formula interface was used.
#'   }
#'   \item{`mixing.copula.candidates`}{
#'     User-requested candidate mixing-copula family codes.
#'   }
#'   \item{`summary`}{
#'     Object of class `"bizimp_selection_summary"` containing the candidate
#'     table, selected candidate, pair summary, and matched call for formatted
#'     display.
#'   }
#'   \item{`call`}{
#'     Matched call to `select_bizimp()`.
#'   }
#' }
#'
#' @seealso
#' [fit_bizimp()] for fitting one specified bivariate model,
#' [select_zimp()] for univariate marginal model selection, and
#' [fit_zimp()] for fitting one specified univariate marginal model.
#'
#' @examples
#' ## Generate a small bivariate Mixed Poisson data set
#' set.seed(12345)
#' n <- 30
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' X <- cbind(1, x1, x2)
#'
#' mu1 <- exp(X %*% c(0.5, 0.3, 0.7))
#' mu2 <- exp(X %*% c(0.3, 0.5, 0.4))
#'
#' ## Gamma-mixed Poisson margins can be generated from negative-binomial
#' ## distributions for this short reproducible example.
#' y1 <- rnbinom(n, size = 1 / 0.2, mu = mu1)
#' y2 <- rnbinom(n, size = 1 / 0.3, mu = mu2)
#' Y <- cbind(y1, y2)
#'
#' ## Compare several mixing-copula families while retaining only the
#' ## independently best marginal model from each margin.
#' sel <- select_bizimp(
#'   Y = Y,
#'   X = X,
#'   Z = NULL,
#'   mixing.distribution = "gamma",
#'   var.power = 2,
#'   margin.top = 1L,
#'   mixing.copula = c(1L, 3L, 5L)
#' )
#'
#' ## Candidate models are ordered by decreasing log-likelihood
#' sel$candidate_table
#'
#' ## Selected marginal and dependence models
#' sel$best_margin_1
#' sel$best_margin_2
#' sel$best_model
#'
#' @export
select_bizimp <- function(
    formula1 = NULL, formula2 = NULL, data = NULL,
    Y = NULL, X = NULL, Z = NULL,
    mixing.distribution = c(
      "gamma", "lnorm", "ig",
      "weibull", "loglogistic", "betaprime"
    ),
    var.power = NULL,
    margin.top = 1L,
    mixing.copula = c(1L),
    d0 = 10L,
    n_gh = 20, gh = NULL,
    ini.par = NULL,
    parallel = NULL,
    progress = FALSE, show.result = FALSE, trace = 0,
    eval.max = 100, rel.tol = 1e-5,
    x.tol = 1e-5, iter.max = 100,
    gaussian.sign.tol = 1e-6,
    t.df.grid = c(3, 5, 10, 30)) {

  # ---------------------------------------------------------------------------
  # 0. Validate and normalize input
  # ---------------------------------------------------------------------------
  inp <- .bizimp_resolve_input(formula1, formula2, data, Y, X, Z)

  # select_bizimp() is intentionally silent. These arguments are retained
  # only for backward compatibility; reporting is stored in $summary.
  progress <- FALSE
  show.result <- FALSE
  trace <- 0

  if (length(margin.top) != 1L || !is.numeric(margin.top) ||
      !is.finite(margin.top) || margin.top < 1 ||
      abs(margin.top - round(margin.top)) > 1e-12) {
    stop("'margin.top' must be a positive integer.")
  }
  margin.top <- as.integer(round(margin.top))

  mixing.copula <- unique(as.integer(mixing.copula))
  if (!length(mixing.copula) || any(!is.finite(mixing.copula))) {
    stop("'mixing.copula' must contain at least one valid VineCopula family code.")
  }

  bad_family <- setdiff(mixing.copula, .vine_all_biv_families)
  if (length(bad_family)) {
    stop(
      "Unsupported VineCopula family code(s): ",
      paste(bad_family, collapse = ", ")
    )
  }

  if (is.null(ini.par)) ini.par <- list(NULL, NULL)
  if (!is.list(ini.par) || length(ini.par) != 2L) {
    stop("'ini.par' must be NULL or a list of length 2.")
  }

  if (!is.null(parallel)) {
    if (length(parallel) != 1L || !is.numeric(parallel) ||
        !is.finite(parallel) || parallel < 1 ||
        abs(parallel - round(parallel)) > 1e-12) {
      stop("'parallel' must be NULL or a positive integer.")
    }
    parallel <- as.integer(round(parallel))
  }

  if (length(show.result) != 1L || !is.logical(show.result) || is.na(show.result)) {
    stop("'show.result' must be TRUE or FALSE.")
  }

  # Number of public marginal candidates evaluated by select_zimp().
  dist_req <- unique(tolower(as.character(mixing.distribution)))
  if (is.null(var.power)) {
    n_margin_models <- length(dist_req)
  } else {
    vp_req <- unique(as.integer(var.power))
    n_margin_models <- length(dist_req) * length(vp_req)
  }

  run_txt <- if (is.null(parallel) || parallel == 1L) {
    "serial"
  } else {
    paste0("parallel | ", parallel, " cores")
  }

  if (isTRUE(progress)) {
    cat("\nBivariate model selection | ", run_txt, "\n", sep = "")
  }

  # ---------------------------------------------------------------------------
  # 1. Marginal selection -- same select_zimp() engine and output philosophy
  # ---------------------------------------------------------------------------
  select_one_margin <- function(j) {
    args <- list(
      mixing.distribution = mixing.distribution,
      var.power = var.power,
      n_gh = n_gh,
      gh = gh,
      ini.par = ini.par[[j]],
      parallel = parallel,
      progress = FALSE,
      show.result = FALSE,
      trace = trace,
      eval.max = eval.max,
      rel.tol = rel.tol,
      x.tol = x.tol,
      iter.max = iter.max,
      .return.all = TRUE
    )

    if (inp$interface == "formula") {
      do.call(
        select_zimp,
        c(list(formula = if (j == 1L) inp$formula1 else inp$formula2, data = inp$data), args)
      )
    } else {
      do.call(
        select_zimp,
        c(list(
          Y = inp$Y[[j]],
          X = inp$X[[j]],
          Z = inp$Z[[j]]
        ), args)
      )
    }
  }

  if (isTRUE(progress)) {
    cat("\nFitting margin 1: ", n_margin_models, " models | ", run_txt, "\n", sep = "")
    flush.console()
  }
  sel1 <- select_one_margin(1L)

  if (isTRUE(progress)) {
    cat("Fitting margin 2: ", n_margin_models, " models | ", run_txt, "\n", sep = "")
    flush.console()
  }
  sel2 <- select_one_margin(2L)

  # Both margins must refer to the same observations.
  if (length(sel1$best_model$y) != length(sel2$best_model$y)) {
    stop("The two margins must use the same number of observations.")
  }
  if (inp$interface == "formula" &&
      !identical(sel1$best_model$na_removed, sel2$best_model$na_removed)) {
    stop(
      "The two formulas removed different rows; ",
      "bivariate fitting requires identical observations."
    )
  }

  if (is.null(sel1$all_results) || !length(sel1$all_results) || is.null(sel2$all_results) || !length(sel2$all_results)) stop("Internal marginal candidate results are unavailable. Ensure select_zimp() and select_bizimp() are from the same package version.")

  m1 <- min(margin.top, length(sel1$all_results))
  m2 <- min(margin.top, length(sel2$all_results))

  top1 <- sel1$all_results[seq_len(m1)]
  top2 <- sel2$all_results[seq_len(m2)]

  if (isTRUE(progress)) {
    cat("\nSelecting dependence structure\n")
    cat("Marginal pairs to evaluate: ", m1, " x ", m2,
        " = ", m1 * m2, "\n", sep = "")
    cat("Candidate mixing copulas: ",
        paste(vapply(
          mixing.copula,
          function(f) .bizimp_public_copula_name(f),
          character(1)
        ), collapse = ", "),
        "\n", sep = "")
  }

  # ---------------------------------------------------------------------------
  # 2. Dependence selection for each retained marginal pair
  # ---------------------------------------------------------------------------
  pair_results <- vector("list", m1 * m2)
  pair_summary <- vector("list", m1 * m2)
  pair_id <- 0L

  for (i in seq_len(m1)) {
    for (j in seq_len(m2)) {
      pair_id <- pair_id + 1L
      f1 <- top1[[i]]
      f2 <- top2[[j]]

      if (isTRUE(progress)) {
        cat(
          "\n[Marginal pair ", pair_id, "/", m1 * m2, "] ",
          "M1 rank ", i, " x M2 rank ", j, "\n",
          sep = ""
        )
      }

      both_zero <- isTRUE(f1$zero_inflation) &&
                   isTRUE(f2$zero_inflation)

      # -----------------------------------------------------------------------
      # 2A. Gaussian prescreen for the SIGN of the mixing dependence only.
      #
      # IMPORTANT:
      #   The zero-copula choice is NOT fixed by this prescreen.
      #   For every retained mixing copula below we fit all admissible
      #   zero-copula specifications and keep them as separate candidates.
      # -----------------------------------------------------------------------
      zero_options <- if (both_zero) {
        c("normal", "comonotone", "countermonotone")
      } else {
        "normal"
      }

      pilot <- vector("list", length(zero_options))
      names(pilot) <- zero_options

      if (isTRUE(progress)) {
        cat("  Prescreen: Gaussian checkerboard mixing copula\n")
      }

      for (b in seq_along(zero_options)) {
        zopt <- zero_options[b]

        # Zero-copula variants are fitted internally but intentionally not
        # printed here; this keeps the dependence-search progress compact.
        if (isTRUE(progress) && b == 1L) {
          cat("  Fitting Gaussian checkerboard copula ...\n")
          flush.console()
        }

        pilot[[b]] <- .fit_bizimp_dependence_engine(
          fit_margin_1 = f1,
          fit_margin_2 = f2,
          n_gh = n_gh,
          gh = gh,
          exact.copula.family = 1L,
          B.depen.option = zopt,
          d0 = d0,
          gaussian.sign.tol = gaussian.sign.tol,
          t.df.grid = t.df.grid,
          progress = FALSE,
          trace = trace,
          eval.max = eval.max,
          iter.max = iter.max
        )
      }

      pilot_ll <- vapply(
        pilot,
        function(z) z$best_model$log_lik_value,
        numeric(1)
      )
      pilot_ll[!is.finite(pilot_ll)] <- -Inf

      # Use the best Gaussian pilot only to determine the sign used for
      # filtering incompatible mixing-copula families.
      pilot_best_id <- which.max(pilot_ll)
      pilot_best <- pilot[[pilot_best_id]]
      gaussian_rho <- pilot_best$best_model$par

      # Use the Gaussian pilot only as a sign prescreen. If the estimated
      # Gaussian dependence is not clearly different from zero, do not impose
      # any sign restriction. Otherwise retain only user-requested families
      # compatible with the estimated sign. If that filtering removes every
      # requested family, fall back to fitting all user-requested families.
      gaussian_rho <- as.numeric(gaussian_rho)[1L]
      gaussian_significant <- is.finite(gaussian_rho) && abs(gaussian_rho) > gaussian.sign.tol

        # Gaussian pilot is used only to determine the dependence sign.
      gaussian_rho <- as.numeric(gaussian_rho)[1L]

      # If Gaussian dependence is not finite or is too close to zero,
      # the sign is not informative: retain all user-requested families.
      if (!is.finite(gaussian_rho) || abs(gaussian_rho) <= gaussian.sign.tol) {
        gaussian_significant <- FALSE
        compatible_families <- mixing.copula
      } else {
        gaussian_significant <- TRUE

        # Filter only the families explicitly requested by the user.
        compatible_families <- .vine_filter_by_gaussian_sign(
          mixing.copula, gaussian_rho, gaussian.sign.tol
        )
        compatible_families <- as.integer(compatible_families)
        compatible_families <- compatible_families[
          !is.na(compatible_families) & compatible_families %in% mixing.copula
        ]
        compatible_families <- unique(compatible_families)

        # If sign filtering removes every requested family, fall back to all
        # user-requested families rather than introducing the Gaussian pilot.
        if (!length(compatible_families)) compatible_families <- mixing.copula
      }

      # Safety: the Gaussian pilot must never introduce an unrequested family.
      compatible_families <- intersect(compatible_families, mixing.copula)
      if (!length(compatible_families)) compatible_families <- mixing.copula

      # -----------------------------------------------------------------------
      # 2B. Full dependence search:
      #     requested/compatible mixing copula x zero copula.
      #
      # Gaussian family 1 is reused from the pilot only when family 1 was
      # explicitly requested by the user and survives the prescreen.
      # -----------------------------------------------------------------------
      candidate_dep <- list()
      candidate_zero <- character(0)

      if (1L %in% compatible_families) {
        for (b in seq_along(zero_options)) {
          key <- paste0("1:", zero_options[b])
          candidate_dep[[key]] <- pilot[[b]]
          candidate_zero[key] <- zero_options[b]
        }
      }

      other_families <- setdiff(compatible_families, 1L)
      if (length(other_families)) {
        for (fam in other_families) {
          for (zopt in zero_options) {
            # Fit all zero-copula variants internally, but print each
            # mixing checkerboard family only once.
            if (isTRUE(progress) && identical(zopt, zero_options[1L])) {
              cat(
                "  Fitting ",
                .bizimp_public_copula_name(fam),
                " checkerboard copula ...\n",
                sep = ""
              )
              flush.console()
            }

            key <- paste0(fam, ":", zopt)
            candidate_dep[[key]] <- .fit_bizimp_dependence_engine(
              fit_margin_1 = f1,
              fit_margin_2 = f2,
              n_gh = n_gh,
              gh = gh,
              exact.copula.family = fam,
              B.depen.option = zopt,
              d0 = d0,
              gaussian.sign.tol = gaussian.sign.tol,
              t.df.grid = t.df.grid,
              progress = FALSE,
              trace = trace,
              eval.max = eval.max,
              iter.max = iter.max
            )
            candidate_zero[key] <- zopt
          }
        }
      }

      if (!length(candidate_dep)) stop("No compatible dependence candidate could be fitted for a retained marginal pair.")

      candidate_ll <- vapply(
        candidate_dep,
        function(z) z$best_model$log_lik_value,
        numeric(1)
      )
      candidate_ll[!is.finite(candidate_ll) |
                   candidate_ll <= -1e12 + 1] <- -Inf
      best_candidate <- which.max(candidate_ll)
      dep <- candidate_dep[[best_candidate]]
      bm <- dep$best_model
      if (!(as.integer(bm$family) %in% mixing.copula)) {
        stop("Internal error: selected mixing copula is not among the user-requested candidate families.")
      }
      selected_zero <- if (both_zero) {
        candidate_zero[[names(candidate_dep)[best_candidate]]]
      } else {
        NA_character_
      }
      zero_info <- .bizimp_selection_zero_label(
        selected_zero,
        bm$theta.B
      )
      mixing_name <- .bizimp_public_copula_name(
        bm$family, fallback = bm$family_name
      )

      pair_results[[pair_id]] <- list(
        margin1_rank = i,
        margin2_rank = j,
        margin1 = f1,
        margin2 = f2,
        pilot = pilot,
        selected_zero = selected_zero,
        gaussian_rho = gaussian_rho,
        gaussian_significant = gaussian_significant,
        compatible_families = compatible_families,
        candidate_dependence = candidate_dep,
        candidate_zero = candidate_zero,
        dependence = dep
      )

      pair_summary[[pair_id]] <- data.frame(
        pair = pair_id,
        margin1_rank = i,
        margin2_rank = j,
        margin1_distribution = f1$mixing.distribution,
        margin1_var_power = .bizimp_selection_display_var_power(f1),
        margin2_distribution = f2$mixing.distribution,
        margin2_var_power = .bizimp_selection_display_var_power(f2),
        gaussian_prescreen_rho = gaussian_rho,
        zero_copula = zero_info$name,
        zero_parameter = zero_info$parameter,
        mixing_copula_family = bm$family,
        mixing_copula = mixing_name,
        mixing_parameter = .bizimp_selection_copula_parameter(bm),
        logLik = bm$log_lik_value,
        convergence = bm$convergence,
        stringsAsFactors = FALSE
      )

    }
  }

  # Full candidate table: retained marginal choices x zero copula x mixing copula.
  candidate_rows <- list()
  rr <- 0L
  for (pp in seq_along(pair_results)) {
    pr <- pair_results[[pp]]
    f1 <- pr$margin1
    f2 <- pr$margin2
    both_zero <- isTRUE(f1$zero_inflation) && isTRUE(f2$zero_inflation)

    dep_names <- names(pr$candidate_dependence)
    for (kk in seq_along(pr$candidate_dependence)) {
      dd <- pr$candidate_dependence[[kk]]
      mm <- dd$best_model
      actual_zero <- if (both_zero) {
        if (!is.null(pr$candidate_zero) &&
            dep_names[kk] %in% names(pr$candidate_zero)) {
          unname(pr$candidate_zero[[dep_names[kk]]])
        } else if (!is.null(dd$B.depen.option)) {
          dd$B.depen.option
        } else if (!is.null(mm$B.depen.option)) {
          mm$B.depen.option
        } else {
          NA_character_
        }
      } else {
        NA_character_
      }
      zi <- .bizimp_selection_zero_short(
        actual_zero,
        mm$theta.B
      )
      rr <- rr + 1L
      ll <- mm$log_lik_value
      if (!is.finite(ll) || ll <= -1e12 + 1) ll <- -Inf

      zero_option_for_ic <- if (both_zero) actual_zero else NA_character_
      k_total <- .bizimp_total_npar(
        f1, f2, mm, zero_option_for_ic
      )
      ic <- .bizimp_information_criteria(
        ll, length(f1$y), k_total
      )

      candidate_rows[[rr]] <- data.frame(
        margin.1 = .bizimp_selection_margin_compact(f1),
        margin.2 = .bizimp_selection_margin_compact(f2),
        zero.copula = .bizimp_selection_zero_compact(
          actual_zero, mm$theta.B
        ),
        mixing.copula = .bizimp_selection_mixing_compact(mm),
        logLik = ll,
        AIC = unname(ic[["AIC"]]),
        BIC = unname(ic[["BIC"]]),
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }
  if (!length(candidate_rows)) stop("No valid bivariate candidate models were fitted.")
  candidate_table <- do.call(rbind, candidate_rows)
  candidate_table <- candidate_table[
    order(candidate_table$logLik, decreasing = TRUE, na.last = TRUE),
    , drop = FALSE
  ]
  rownames(candidate_table) <- NULL

  pair_table <- do.call(rbind, pair_summary)
  ord <- order(pair_table$logLik, decreasing = TRUE, na.last = TRUE)
  pair_table <- pair_table[ord, , drop = FALSE]
  rownames(pair_table) <- NULL

  best_pair_id <- pair_table$pair[1L]
  best_pair <- pair_results[[best_pair_id]]
  best_dep <- best_pair$dependence
  best_model <- best_dep$best_model
  best_f1 <- best_pair$margin1
  best_f2 <- best_pair$margin2

  best_k <- .bizimp_total_npar(
    best_f1, best_f2, best_model,
    if (isTRUE(best_f1$zero_inflation) &&
        isTRUE(best_f2$zero_inflation)) best_pair$selected_zero else NA_character_
  )
  best_ic <- .bizimp_information_criteria(
    best_model$log_lik_value, length(best_f1$y), best_k
  )
  best_model$npar <- best_k
  best_model$AIC <- unname(best_ic[["AIC"]])
  best_model$BIC <- unname(best_ic[["BIC"]])

  both_zero_best <- isTRUE(best_f1$zero_inflation) &&
                    isTRUE(best_f2$zero_inflation)

  zero_info <- .bizimp_selection_zero_label(
    if (both_zero_best) best_pair$selected_zero else NA_character_,
    best_model$theta.B
  )
  mixing_name <- .bizimp_public_copula_name(
    best_model$family,
    fallback = best_model$family_name
  )

  # ---------------------------------------------------------------------------
  # 3. Public result -- silent; report is stored in $summary
  # ---------------------------------------------------------------------------
  selection_summary <- list(
    summary_table = candidate_table,
    best_model = candidate_table[1L, , drop = FALSE],
    pair_summary = pair_table,
    call = match.call()
  )
  class(selection_summary) <- "bizimp_selection_summary"

  out <- list(
    margin_selection_1 = sel1,
    margin_selection_2 = sel2,
    margin_top = margin.top,
    retained_margin_1 = top1,
    retained_margin_2 = top2,
    pair_results = pair_results,
    candidate_table = candidate_table,
    summary_table = candidate_table,
    pair_summary = pair_table,
    best_pair = best_pair,
    best_model = best_model,
    best_dependence = best_dep,
    best_margin_1 = best_f1,
    best_margin_2 = best_f2,
    interface = inp$interface,
    formula1 = inp$formula1,
    formula2 = inp$formula2,
    formula = if (inp$interface == "formula") list(inp$formula1, inp$formula2) else NULL,
    mixing.copula.candidates = mixing.copula,
    summary = selection_summary,
    call = match.call()
  )
  class(out) <- "bizimp_selection"
  out
}

#' Print a select_bizimp summary
#'
#' @param x An object of class `bizimp_selection_summary`.
#' @param ... Additional arguments (currently unused).
#' @return The summary object, invisibly.
#' @export
print.bizimp_selection_summary <- function(x, ...) {
  cat("\nBivariate model selection completed\n")
  cat(strrep("=", 72), "\n", sep = "")
  cat("\nCandidate models\n")
  cat(strrep("-", 72), "\n", sep = "")
  print(x$summary_table, row.names = FALSE)
  cat("\nBest model\n")
  cat(strrep("-", 72), "\n", sep = "")
  print(x$best_model, row.names = FALSE)
  cat(strrep("=", 72), "\n", sep = "")
  invisible(x)
}
