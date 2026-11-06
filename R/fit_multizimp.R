# ============================================================
# make_factor_loadings
#
# Convert an identifiable raw parameter vector into a
# normalized lower-triangular factor-loading matrix.
#
# A[j,r] = 0           for r > j
# A[j,j] = exp(raw)    positive diagonal for identification
# A[j,r] = raw         for r < j
#
# Psi[j,] = A[j,] / sqrt(1 + ||A[j,]||^2)
#
# Therefore every row satisfies:
#   sum(Psi[j,]^2) < 1
# ============================================================

make_factor_loadings <- function(
    par,
    k,
    q,
    return_all = FALSE
) {

  if (!is.numeric(k) ||
      length(k) != 1L ||
      !is.finite(k) ||
      k < 1L ||
      abs(k - round(k)) > 1e-12) {
    stop("`k` must be a positive integer.")
  }

  if (!is.numeric(q) ||
      length(q) != 1L ||
      !is.finite(q) ||
      q < 1L ||
      abs(q - round(q)) > 1e-12) {
    stop("`q` must be a positive integer.")
  }

  k <- as.integer(round(k))
  q <- as.integer(round(q))

  if (q > k) {
    stop("`q` must not exceed `k`.")
  }

  npar <- k * q - q * (q - 1L) / 2L

  if (!is.numeric(par) ||
      length(par) != npar ||
      any(!is.finite(par))) {
    stop(
      "`par` must be a finite numeric vector of length ",
      npar, "."
    )
  }

  A <- matrix(
    0,
    nrow = k,
    ncol = q
  )

  map <- matrix(
    NA_integer_,
    nrow = npar,
    ncol = 2L,
    dimnames = list(
      NULL,
      c("margin", "factor")
    )
  )

  id <- 1L

  for (j in seq_len(k)) {

    for (r in seq_len(q)) {

      if (r <= j) {

        if (j == r) {

          # Positive diagonal for identification.
          A[j, r] <- exp(par[id])

        } else {

          A[j, r] <- par[id]
        }

        map[id, ] <- c(j, r)

        id <- id + 1L
      }
    }
  }

  denom <- sqrt(
    1 + rowSums(A^2)
  )

  Psi <- A / denom

  if (any(!is.finite(Psi))) {
    stop("Non-finite factor loadings were generated.")
  }

  if (any(rowSums(Psi^2) >= 1)) {
    stop("Invalid factor loadings: row norm must be strictly below 1.")
  }

  if (isTRUE(return_all)) {

    return(
      list(
        Psi = Psi,
        A = A,
        map = map,
        npar = npar
      )
    )
  }

  Psi
}


# ============================================================
# make_factor_start
#
# Construct an identifiable starting raw parameter vector.
#
# Diagonal loadings of the first q rows are initialized to
# `loading0`; all off-diagonal loadings start at zero.
# ============================================================

make_factor_start <- function(
    k,
    q,
    loading0 = 0.05
) {

  if (!is.numeric(k) ||
      length(k) != 1L ||
      !is.finite(k) ||
      k < 1L ||
      abs(k - round(k)) > 1e-12) {
    stop("`k` must be a positive integer.")
  }

  if (!is.numeric(q) ||
      length(q) != 1L ||
      !is.finite(q) ||
      q < 1L ||
      abs(q - round(q)) > 1e-12) {
    stop("`q` must be a positive integer.")
  }

  k <- as.integer(round(k))
  q <- as.integer(round(q))

  if (q > k) {
    stop("`q` must not exceed `k`.")
  }

  if (!is.numeric(loading0) ||
      length(loading0) != 1L ||
      !is.finite(loading0) ||
      loading0 <= 0 ||
      loading0 >= 1) {
    stop("`loading0` must lie strictly between 0 and 1.")
  }

  npar <- k * q - q * (q - 1L) / 2L

  par0 <- numeric(npar)

  # If Psi_jj = loading0, then
  #
  # loading0 = A_jj / sqrt(1 + A_jj^2)
  #
  # hence:
  A0_diag <-
    loading0 /
    sqrt(1 - loading0^2)

  id <- 1L

  for (j in seq_len(k)) {

    for (r in seq_len(q)) {

      if (r <= j) {

        if (j == r) {

          # Since A_jj = exp(par_jj).
          par0[id] <- log(A0_diag)

        } else {

          par0[id] <- 0
        }

        id <- id + 1L
      }
    }
  }

  par0
}



# ============================================================
# make_factor_start_nonzero
#
# Analytic Hermite-polynomial gradients use psi_j = r_j u_j and therefore
# require r_j = ||psi_j|| > 0 at the starting point.
#
# The generic triangular start leaves rows j > q exactly zero.
# This helper keeps the identifying diagonal start but gives every
# free row a tiny deterministic non-zero loading pattern.
# ============================================================

make_factor_start_nonzero <- function(
    k,
    q,
    loading0 = 0.05,
    offdiag0 = 0.02
) {

  par0 <- make_factor_start(
    k = k,
    q = q,
    loading0 = loading0
  )

  id <- 1L

  for (j in seq_len(k)) {
    for (r in seq_len(q)) {

      if (r <= j) {

        ## For rows beyond the identifying q x q block, every free
        ## element is off-diagonal under the triangular convention.
        ## Give them a small deterministic alternating pattern so
        ## ||psi_j|| is strictly positive from the first evaluation.
        if (j > q) {
          par0[id] <- offdiag0 * ifelse((j + r) %% 2L == 0L, 1, -1)
        }

        id <- id + 1L
      }
    }
  }

  par0
}



# ============================================================
# upgrade_factor_start
#
# Warm-start a larger Gaussian factor model from a fitted smaller
# loading matrix.  Existing fitted columns are retained exactly at
# the starting point; newly added columns receive small deterministic
# non-zero loadings.  The complete enlarged Psi is then transformed
# back to the raw lower-triangular parameterization used by
# make_factor_loadings().
# ============================================================

upgrade_factor_start <- function(
    Psi_old,
    q_new,
    loading0 = 0.05,
    offdiag0 = 0.02,
    max_row_norm = 0.98
) {

  Psi_old <- as.matrix(Psi_old)
  k <- nrow(Psi_old)
  q_old <- ncol(Psi_old)
  q_new <- as.integer(round(q_new))

  if (q_new <= q_old) {
    stop("`q_new` must be larger than the number of columns in `Psi_old`.")
  }
  if (q_new > k) stop("`q_new` must not exceed k.")
  if (any(!is.finite(Psi_old))) stop("`Psi_old` contains non-finite values.")
  if (any(rowSums(Psi_old^2) >= 1)) stop("Rows of `Psi_old` must have norm < 1.")

  if (!is.numeric(max_row_norm) || length(max_row_norm) != 1L ||
      !is.finite(max_row_norm) || max_row_norm <= 0 || max_row_norm >= 1) {
    stop("`max_row_norm` must lie strictly between 0 and 1.")
  }

  Psi_new <- matrix(0, nrow = k, ncol = q_new)
  Psi_new[, seq_len(q_old)] <- Psi_old

  # Add small deterministic non-zero loadings for the new factors.
  # Positive diagonal entries preserve the triangular identification.
  for (r in (q_old + 1L):q_new) {
    for (j in r:k) {
      if (j == r) {
        Psi_new[j, r] <- loading0
      } else {
        Psi_new[j, r] <- offdiag0 * ifelse((j + r) %% 2L == 0L, 1, -1)
      }
    }
  }

  # Keep every enlarged row safely inside the unit sphere.
  #
  # Important:
  # For rows q_old < j <= q_new, the newly added diagonal Psi_new[j,j]
  # MUST stay strictly positive for triangular identification.
  #
  # The previous implementation protected the old coordinates at all costs
  # and, if an old fitted row already exceeded `max_row_norm`, set ALL new
  # coordinates to zero.  That can zero out a newly required diagonal and
  # causes:
  #
  #   "Invalid non-positive diagonal while upgrading factor start."
  #
  # New rule:
  #   1) If the old part still leaves room, shrink only the new part.
  #   2) If the old part already uses the whole norm budget, retain the
  #      newly initialized coordinates and shrink the OLD part just enough
  #      to make room for them.
  #
  # Hence the new identifying diagonal is never deliberately zeroed.
  old_idx <- seq_len(q_old)
  new_idx <- (q_old + 1L):q_new
  max2 <- max_row_norm^2

  for (j in seq_len(k)) {

    old2 <- sum(Psi_new[j, old_idx, drop = TRUE]^2)
    new2 <- sum(Psi_new[j, new_idx, drop = TRUE]^2)

    if (new2 <= 0 || old2 + new2 < max2) next

    if (old2 < max2) {

      # There is still some room after retaining the fitted old loadings.
      # Scale only the newly added coordinates.
      available2 <- max(max2 - old2, 0)

      if (available2 > 0) {
        Psi_new[j, new_idx] <-
          Psi_new[j, new_idx] * sqrt(available2 / new2)
      } else {
        # This branch should be unreachable when old2 < max2, but avoid
        # ever zeroing a required new diagonal because of round-off.
        target_new2 <- min(new2, max2 * 1e-6)
        max_old2 <- max(max2 - target_new2, 0)
        if (old2 > 0 && max_old2 < old2) {
          Psi_new[j, old_idx] <-
            Psi_new[j, old_idx] * sqrt(max_old2 / old2)
        }
        if (new2 > 0 && target_new2 < new2) {
          Psi_new[j, new_idx] <-
            Psi_new[j, new_idx] * sqrt(target_new2 / new2)
        }
      }

    } else {

      # The fitted old coordinates already consume the norm budget.
      # Preserve the newly initialized block (especially its positive
      # identifying diagonal) and shrink the OLD block just enough.
      #
      # `new2` is small by construction, but guard against an extreme
      # user-supplied loading0/offdiag0 nevertheless.
      if (new2 >= max2) {
        target_new2 <- 0.25 * max2
        Psi_new[j, new_idx] <-
          Psi_new[j, new_idx] * sqrt(target_new2 / new2)
        new2 <- target_new2
      }

      target_old2 <- max(max2 - new2, 0)

      if (old2 > 0 && target_old2 < old2) {
        Psi_new[j, old_idx] <-
          Psi_new[j, old_idx] * sqrt(target_old2 / old2)
      }
    }
  }

  # Enforce lower-triangular zeros exactly and verify that all newly
  # introduced identifying diagonals remain strictly positive.
  for (j in seq_len(k)) {
    if (j < q_new) {
      Psi_new[j, (j + 1L):q_new] <- 0
    }
  }

  new_diag_rows <- (q_old + 1L):q_new
  new_diag <- Psi_new[cbind(new_diag_rows, new_diag_rows)]

  if (any(!is.finite(new_diag)) || any(new_diag <= 0)) {
    stop(
      "Unable to preserve a positive identifying diagonal while upgrading ",
      "factor start from q=", q_old, " to q=", q_new, "."
    )
  }

  # Numerical safety: if a row lands infinitesimally above the requested
  # cap because of floating-point rounding, scale the complete row by a
  # negligible amount.  Positive diagonals stay positive.
  row_norm2 <- rowSums(Psi_new^2)
  ii_cap <- which(row_norm2 >= max2)
  if (length(ii_cap) > 0L) {
    target_norm <- max_row_norm * (1 - 1e-10)
    for (j in ii_cap) {
      rn <- sqrt(row_norm2[j])
      if (is.finite(rn) && rn > 0) {
        Psi_new[j, ] <- Psi_new[j, ] * target_norm / rn
      }
    }
  }

  # Convert the desired enlarged Psi exactly back to A, row by row:
  #   Psi_j = A_j / sqrt(1 + ||A_j||^2)
  # => A_j = Psi_j / sqrt(1 - ||Psi_j||^2).
  row_norm2 <- rowSums(Psi_new^2)

  if (any(!is.finite(row_norm2)) || any(row_norm2 >= 1)) {
    stop("Unable to construct a valid enlarged loading start.")
  }

  A_new <- Psi_new / sqrt(1 - row_norm2)

  npar <- k * q_new - q_new * (q_new - 1L) / 2L
  par_new <- numeric(npar)
  id <- 1L

  for (j in seq_len(k)) {
    for (r in seq_len(q_new)) {
      if (r <= j) {
        if (j == r) {
          if (!is.finite(A_new[j, r]) || A_new[j, r] <= 0) {
            stop("Invalid non-positive diagonal while upgrading factor start.")
          }
          par_new[id] <- log(A_new[j, r])
        } else {
          par_new[id] <- A_new[j, r]
        }
        id <- id + 1L
      }
    }
  }

  attr(par_new, "Psi_start") <- Psi_new
  attr(par_new, "upgrade_max_row_norm") <- max(sqrt(row_norm2))
  attr(par_new, "upgrade_old_q") <- q_old
  attr(par_new, "upgrade_new_q") <- q_new

  par_new
}


# ============================================================
# factor_start_from_Psi
#
# Convert an already fitted Psi at the SAME factor dimension back to
# the raw lower-triangular parameterization.  Unlike
# upgrade_factor_start(), this does not add a new factor column.
# ============================================================
factor_start_from_Psi <- function(Psi, max_row_norm = 0.999999) {
  Psi <- as.matrix(Psi)
  k <- nrow(Psi)
  q <- ncol(Psi)
  if (q < 1L || q > k) stop("Invalid factor dimension in `Psi`.")
  if (any(!is.finite(Psi))) stop("`Psi` contains non-finite values.")

  # Enforce the identified triangular zeros numerically.
  for (j in seq_len(k)) {
    if (j < q) Psi[j, (j + 1L):q] <- 0
  }

  rn2 <- rowSums(Psi^2)
  if (any(rn2 >= 1)) stop("Rows of `Psi` must have norm < 1.")
  if (any(rn2 > max_row_norm^2)) {
    ii <- which(rn2 > max_row_norm^2)
    for (j in ii) Psi[j, ] <- Psi[j, ] * max_row_norm / sqrt(rn2[j])
    rn2 <- rowSums(Psi^2)
  }

  A <- Psi / sqrt(1 - rn2)
  npar <- k * q - q * (q - 1L) / 2L
  par <- numeric(npar)
  id <- 1L
  for (j in seq_len(k)) {
    for (r in seq_len(q)) {
      if (r <= j) {
        if (j == r) {
          if (!is.finite(A[j, r]) || A[j, r] <= 0)
            stop("Invalid non-positive diagonal while encoding factor start.")
          par[id] <- log(A[j, r])
        } else {
          par[id] <- A[j, r]
        }
        id <- id + 1L
      }
    }
  }
  attr(par, "Psi_start") <- Psi
  par
}



# ============================================================
# .make_biexp_factor_backend
#
# Internal computational backend for standardized bilateral
# exponential (BiExp) factor copulas.
#
# Parameterization:
#
#   Z = a E1 - b E2 - (a-b),
#   a = lambda,
#   b = sqrt(1-lambda^2),
#
# with E(Z)=0 and Var(Z)=1.
#
# This factory reuses the exact finite-sum distribution and
# analytic loading derivatives from the validated BiExp code.
# ============================================================

.make_biexp_factor_backend <- function(
    lambda,
    q,
    d0,
    gl_order = 15,
    root_tol = 1e-10
) {

  if (!is.numeric(lambda) ||
      length(lambda) != 1L ||
      !is.finite(lambda) ||
      lambda <= 0 ||
      lambda >= 1) {
    stop("`lambda` must lie strictly in (0,1).")
  }

  if (!is.numeric(q) ||
      length(q) != 1L ||
      !is.finite(q) ||
      q < 1L ||
      abs(q - round(q)) > 1e-12) {
    stop("`q` must be a positive integer.")
  }

  q <- as.integer(round(q))

  # Keep the tested code's notation internally.
  J <- q

  a0 <- lambda
  b0 <- sqrt(1 - lambda^2)
  mu0 <- a0 - b0

  z0 <- -mu0
  p0 <- b0 / (a0 + b0)

    preserve_dim <- function(out, x) {
      if (!is.null(dim(x))) {
        dim(out) <- dim(x)
      }
      out
    }
  
  
    ## ------------------------------------------------------------
    ## CDF
    ## ------------------------------------------------------------
  
    F_sd_biexp <- function(z) {
    
      zz <- as.numeric(z)
      out <- numeric(length(zz))
    
      neg <- zz < z0
      pos <- !neg
    
      if (any(neg)) {
        out[neg] <-
          b0 / (a0 + b0) *
          exp((zz[neg] + mu0) / b0)
      }
    
      if (any(pos)) {
        out[pos] <-
          1 -
          a0 / (a0 + b0) *
          exp(-(zz[pos] + mu0) / a0)
      }
    
      out <- pmin(pmax(out, 0), 1)
    
      preserve_dim(out, z)
    }
  
  
    ## ------------------------------------------------------------
    ## Density
    ## ------------------------------------------------------------
  
    f_sd_biexp <- function(z) {
    
      zz <- as.numeric(z)
      out <- numeric(length(zz))
    
      neg <- zz < z0
      pos <- !neg
    
      if (any(neg)) {
        out[neg] <-
          1 / (a0 + b0) *
          exp((zz[neg] + mu0) / b0)
      }
    
      if (any(pos)) {
        out[pos] <-
          1 / (a0 + b0) *
          exp(-(zz[pos] + mu0) / a0)
      }
    
      preserve_dim(out, z)
    }
  
  
    ## ------------------------------------------------------------
    ## Quantile
    ## ------------------------------------------------------------
  
    Q_sd_biexp <- function(u) {
    
      uu <- as.numeric(u)
    
      ## Gauss-Legendre nodes never equal 0 or 1,
      ## but clamp for safety.
      uu <- pmin(pmax(uu, 1e-15), 1 - 1e-15)
    
      out <- numeric(length(uu))
    
      neg <- uu < p0
      pos <- !neg
    
      if (any(neg)) {
        out[neg] <-
          b0 *
          log(
            (a0 + b0) * uu[neg] / b0
          ) -
          mu0
      }
    
      if (any(pos)) {
        out[pos] <-
          -a0 *
          log(
            (a0 + b0) *
              (1 - uu[pos]) / a0
          ) -
          mu0
      }
    
      preserve_dim(out, u)
    }
  
  

    make_gl_grid <- function(J, gl_order) {
    
      gl <- statmod::gauss.quad(
        gl_order,
        kind = "legendre"
      )
    
      ## [-1,1] -> [0,1]
      u1 <- (gl$nodes + 1) / 2
      w1 <- gl$weights / 2
    
      idx <- expand.grid(
        rep(list(seq_len(gl_order)), J)
      )
    
      U <- as.matrix(
        sapply(
          seq_len(J),
          function(r) u1[idx[[r]]]
        )
      )
    
      W <- apply(
        as.matrix(
          sapply(
            seq_len(J),
            function(r) w1[idx[[r]]]
          )
        ),
        1,
        prod
      )
    
      Z <- Q_sd_biexp(U)
    
      list(
        u = U,
        z = Z,
        w = W
      )
    }
  
  
    gl_grid <- make_gl_grid(
      J = J,
      gl_order = gl_order
    )
  
    Znode <- gl_grid$z
    Wnode <- gl_grid$w
  
    n_node <- length(Wnode)
  
  

    get_pf_coef_deriv <- function(
      scale,
      dscale,
      sep_tol = 1e-9
    ) {
    
      L <- length(scale)
    
      if (L == 0L) {
        stop("No active exponential scales.")
      }
    
      if (any(!is.finite(scale)) ||
          any(scale <= 0)) {
        stop("All transformed scales must be positive.")
      }
    
    
      ## ----------------------------------------------------------
      ## One scale only
      ## ----------------------------------------------------------
    
      if (L == 1L) {
      
        return(
          list(
            coef = 1,
            dcoef = matrix(
              0,
              nrow = 1,
              ncol = J
            )
          )
        )
      }
    
    
      ## ----------------------------------------------------------
      ## Exact partial-fraction formula assumes distinct scales
      ## ----------------------------------------------------------
    
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
            "The simple partial-fraction representation becomes ",
            "numerically unstable."
          )
        )
      }
    
    
      ## ----------------------------------------------------------
      ## Coefficients
      ## ----------------------------------------------------------
    
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
    
    
      ## ----------------------------------------------------------
      ## G[i,r] =
      ##
      ## d log(A_i) / d scale_r
      ## ----------------------------------------------------------
    
      G <- matrix(
        0,
        nrow = L,
        ncol = L
      )
    
      for (i in seq_len(L)) {
      
        G[i, i] <-
          sum(
            1 / scale[i] -
              1 /
              (
                scale[i] -
                  scale[-i]
              )
          )
      
        for (r in seq_len(L)) {
        
          if (r != i) {
          
            G[i, r] <-
              1 /
              (
                scale[i] -
                  scale[r]
              )
          }
        }
      }
    
    
      ## ----------------------------------------------------------
      ## Chain rule:
      ##
      ## dA_i/dpsi
      ## =
      ## A_i *
      ## sum_r
      ## [d log(A_i)/d scale_r] *
      ## [d scale_r/dpsi]
      ## ----------------------------------------------------------
    
      dcoef <-
        coef *
        (
          G %*% dscale
        )
    
    
      list(
        coef = coef,
        dcoef = dcoef
      )
    }
  
  
    ## ============================================================
    ## Exact SumSdBiExp object
    ##
    ## Returns exact:
    ##
    ## F_S(x)
    ## f_S(x)
    ## dF_S(x)/dpsi
    ##
    ## for
    ##
    ## S =
    ## sum_l v_l Z_l
    ##
    ## v =
    ## (psi_1,...,psi_J,sigma)
    ##
    ## ============================================================
  
    make_exact_sum_sd_biexp <- function(
      psi_j,
      zero_tol = 1e-12,
      sep_tol = 1e-9
    ) {
    
      psi_j <- as.numeric(psi_j)
    
      norm2 <- sum(psi_j^2)
    
      if (!is.finite(norm2) ||
          norm2 >= 1) {
        stop("||psi_j|| must be smaller than 1.")
      }
    
    
      ## ----------------------------------------------------------
      ## Idiosyncratic loading
      ## ----------------------------------------------------------
    
      sigma_j <-
        sqrt(
          1 - norm2
        )
    
    
      ## ----------------------------------------------------------
      ## v =
      ## (psi_1,...,psi_J,sigma)
      ## ----------------------------------------------------------
    
      v_full <-
        c(
          psi_j,
          sigma_j
        )
    
    
      ## ----------------------------------------------------------
      ## dv_l / dpsi_r
      ##
      ## for factor coefficients:
      ##
      ## d psi_l / d psi_r = I(l=r)
      ##
      ## and
      ##
      ## d sigma / d psi_r
      ## =
      ## -psi_r/sigma
      ## ----------------------------------------------------------
    
      dv_full <- rbind(
        diag(J),
        -psi_j / sigma_j
      )
    
    
      ## ----------------------------------------------------------
      ## Terms with coefficient exactly zero contribute nothing.
      ##
      ## Structural zeros in the loading matrix are therefore
      ## removed from the SumSdBiExp representation.
      ## ----------------------------------------------------------
    
      active <-
        abs(v_full) >
        zero_tol
    
      v <-
        v_full[active]
    
      dv <-
        dv_full[
          active,
          ,
          drop = FALSE
        ]
    
    
      if (length(v) == 0L) {
        stop("No active coefficients in SumSdBiExp.")
      }
    
    
      ## ==========================================================
      ## Transform each coefficient v_l:
      ##
      ## if v_l > 0:
      ##
      ## alpha_l = a0 * v_l
      ## beta_l  = b0 * v_l
      ##
      ## if v_l < 0:
      ##
      ## alpha_l = -b0 * v_l
      ## beta_l  = -a0 * v_l
      ##
      ## Thus
      ##
      ## Y =
      ## sum alpha_l E_l^+
      ## -
      ## sum beta_l E_l^-
      ##
      ## and
      ##
      ## S =
      ## Y -
      ## mu0 * sum(v_l)
      ## ==========================================================
    
      positive <-
        v > 0
    
    
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
    
    
      ## ----------------------------------------------------------
      ## Analytic derivatives of transformed scales
      ## ----------------------------------------------------------
    
      alpha_multiplier <-
        ifelse(
          positive,
          a0,
          -b0
        )
    
      beta_multiplier <-
        ifelse(
          positive,
          b0,
          -a0
        )
    
    
      dalpha <-
        dv *
        alpha_multiplier
    
      dbeta <-
        dv *
        beta_multiplier
    
    
      ## ==========================================================
      ## Partial-fraction coefficients
      ## ==========================================================
    
      aa <-
        get_pf_coef_deriv(
          scale = alpha,
          dscale = dalpha,
          sep_tol = sep_tol
        )
    
      bb <-
        get_pf_coef_deriv(
          scale = beta,
          dscale = dbeta,
          sep_tol = sep_tol
        )
    
    
      Acoef <- aa$coef
      Bcoef <- bb$coef
    
      dAcoef <- aa$dcoef
      dBcoef <- bb$dcoef
    
    
      Lp <- length(alpha)
      Lm <- length(beta)
    
    
      ## ==========================================================
      ## Standardization shift
      ##
      ## S =
      ## Y - mu0 * sum(v)
      ##
      ## therefore
      ##
      ## y =
      ## x + mu0 * sum(v)
      ## ==========================================================
    
      shift <-
        mu0 *
        sum(v)
    
    
      ## ----------------------------------------------------------
      ## derivative of shift w.r.t psi
      ##
      ## Active terms only.
      ##
      ## For structural zero loadings the first derivative of the
      ## corresponding mean-zero factor contribution at zero is zero.
      ## ----------------------------------------------------------
    
      dshift <-
        mu0 *
        colSums(dv)
    
    
      ## ==========================================================
      ## Exact evaluator
      ##
      ## Returns:
      ##
      ## F = F_S(x)
      ## f = f_S(x)
      ## dF = d F_S(x) / d psi
      ##
      ## ==========================================================
    
      eval_all <- function(x) {
      
        if (length(x) != 1L ||
            !is.finite(x)) {
          stop("x must be a finite scalar.")
        }
      
      
        y <-
          x +
          shift
      
        dy <-
          dshift
      
      
        FS <- 0
        fS <- 0
      
        dFS <-
          numeric(J)
      
      
        ## ========================================================
        ## Left branch: y < 0
        ## ========================================================
      
        if (y < 0) {
        
          for (ii in seq_len(Lp)) {
          
            ai <- alpha[ii]
          
            dai <-
              dalpha[
                ii,
                ,
                drop = TRUE
              ]
          
          
            for (jj in seq_len(Lm)) {
            
              bj <- beta[jj]
            
              dbj <-
                dbeta[
                  jj,
                  ,
                  drop = TRUE
                ]
            
            
              den <-
                ai + bj
            
            
              Ai <- Acoef[ii]
              Bj <- Bcoef[jj]
            
              dAi <-
                dAcoef[
                  ii,
                  ,
                  drop = TRUE
                ]
            
              dBj <-
                dBcoef[
                  jj,
                  ,
                  drop = TRUE
                ]
            
            
              ## --------------------------------------------------
              ## CDF contribution:
              ##
              ## A_i B_j *
              ## beta_j/(alpha_i+beta_j) *
              ## exp(y/beta_j)
              ## --------------------------------------------------
            
              g_cdf <-
                bj /
                den *
                exp(
                  y / bj
                )
            
              term_cdf <-
                Ai *
                Bj *
                g_cdf
            
              FS <-
                FS +
                term_cdf
            
            
              ## --------------------------------------------------
              ## derivative log(g_cdf)
              ##
              ## =
              ## db/b
              ## -
              ## (da+db)/(a+b)
              ## +
              ## dy/b
              ## -
              ## y db/b^2
              ## --------------------------------------------------
            
              dlog_g_cdf <-
                dbj / bj -
                (
                  dai + dbj
                ) / den +
                dy / bj -
                y *
                dbj /
                bj^2
            
            
              dFS <-
                dFS +
                (
                  dAi * Bj +
                    Ai * dBj
                ) *
                g_cdf +
                Ai *
                Bj *
                g_cdf *
                dlog_g_cdf
            
            
              ## --------------------------------------------------
              ## Density contribution:
              ##
              ## A_i B_j /(a_i+b_j) exp(y/b_j)
              ## --------------------------------------------------
            
              g_den <-
                exp(
                  y / bj
                ) /
                den
            
              fS <-
                fS +
                Ai *
                Bj *
                g_den
            }
          }
        
        
          ## ========================================================
          ## Right branch: y >= 0
          ## ========================================================
        
        } else {
        
          tail_prob <- 0
          dtail <- numeric(J)
        
        
          for (ii in seq_len(Lp)) {
          
            ai <- alpha[ii]
          
            dai <-
              dalpha[
                ii,
                ,
                drop = TRUE
              ]
          
          
            for (jj in seq_len(Lm)) {
            
              bj <- beta[jj]
            
              dbj <-
                dbeta[
                  jj,
                  ,
                  drop = TRUE
                ]
            
            
              den <-
                ai + bj
            
            
              Ai <- Acoef[ii]
              Bj <- Bcoef[jj]
            
              dAi <-
                dAcoef[
                  ii,
                  ,
                  drop = TRUE
                ]
            
              dBj <-
                dBcoef[
                  jj,
                  ,
                  drop = TRUE
                ]
            
            
              ## --------------------------------------------------
              ## Survival contribution:
              ##
              ## A_i B_j *
              ## alpha_i/(alpha_i+beta_j) *
              ## exp(-y/alpha_i)
              ## --------------------------------------------------
            
              g_tail <-
                ai /
                den *
                exp(
                  -y / ai
                )
            
              term_tail <-
                Ai *
                Bj *
                g_tail
            
              tail_prob <-
                tail_prob +
                term_tail
            
            
              ## --------------------------------------------------
              ## derivative log(g_tail)
              ##
              ## =
              ## da/a
              ## -
              ## (da+db)/(a+b)
              ## -
              ## dy/a
              ## +
              ## y da/a^2
              ## --------------------------------------------------
            
              dlog_g_tail <-
                dai / ai -
                (
                  dai + dbj
                ) / den -
                dy / ai +
                y *
                dai /
                ai^2
            
            
              dtail <-
                dtail +
                (
                  dAi * Bj +
                    Ai * dBj
                ) *
                g_tail +
                Ai *
                Bj *
                g_tail *
                dlog_g_tail
            
            
              ## --------------------------------------------------
              ## Exact density contribution
              ## --------------------------------------------------
            
              g_den <-
                exp(
                  -y / ai
                ) /
                den
            
              fS <-
                fS +
                Ai *
                Bj *
                g_den
            }
          }
        
        
          FS <-
            1 -
            tail_prob
        
          dFS <-
            -dtail
        }
      
      
        ## ========================================================
        ## Numerical checks only.
        ##
        ## These do NOT approximate F or f.
        ## They only protect against floating-point cancellation.
        ## ========================================================
      
        if (!is.finite(FS) ||
            !is.finite(fS) ||
            any(!is.finite(dFS))) {
          stop("Non-finite exact SumSdBiExp evaluation.")
        }
      
      
        if (FS < -1e-7 ||
            FS > 1 + 1e-7) {
          stop(
            "Exact SumSdBiExp CDF outside [0,1]; possible numerical instability."
          )
        }
      
      
        if (fS <= 0) {
          stop(
            "Non-positive SumSdBiExp density; possible numerical instability."
          )
        }
      
      
        ## Only round-off correction
        FS <-
          min(
            max(FS, 0),
            1
          )
      
      
        list(
          F = FS,
          f = fS,
          dF = dFS
        )
      }
    
    
      ## ----------------------------------------------------------
      ## Convenience wrappers
      ## ----------------------------------------------------------
    
      cdf <- function(x) {
      
        vapply(
          x,
          function(xx) {
            eval_all(xx)$F
          },
          numeric(1)
        )
      }
    
    
      density <- function(x) {
      
        vapply(
          x,
          function(xx) {
            eval_all(xx)$f
          },
          numeric(1)
        )
      }
    
    
      list(
        psi = psi_j,
        sigma = sigma_j,
        v = v,
        alpha = alpha,
        beta = beta,
        Acoef = Acoef,
        Bcoef = Bcoef,
        dAcoef = dAcoef,
        dBcoef = dBcoef,
        shift = shift,
        dshift = dshift,
        eval = eval_all,
        cdf = cdf,
        density = density
      )
    }
  
  
    ## ============================================================
    ## 6. Quantiles Q_S(m/d) and exact analytic dQ/dpsi
    ##
    ## IMPORTANT:
    ##
    ## Q_S itself generally has no closed form.
    ##
    ## It is obtained numerically by solving
    ##
    ## F_S(q;psi) = p,
    ##
    ## but F_S is evaluated by the EXACT analytical finite-sum
    ## expression above.
    ##
    ## No quadrature is used for F_S or f_S.
    ##
    ## Analytic quantile derivative:
    ##
    ## dQ/dpsi
    ## =
    ## -
    ## [dF_S(q;psi)/dpsi] /
    ## f_S(q;psi)
    ## ============================================================
  
    compute_Q_and_dQ <- function(psi_j) {
    
      norm2 <- sum(psi_j^2)
    
      if (!is.finite(norm2) ||
          norm2 >= 1) {
        stop("||psi_j|| must be smaller than 1.")
      }
    
    
      sigma_j <-
        sqrt(
          1 - norm2
        )
    
    
      ## ----------------------------------------------------------
      ## Exact analytic distribution object
      ## ----------------------------------------------------------
    
      Sobj <-
        make_exact_sum_sd_biexp(
          psi_j = psi_j
        )
    
    
      prob_inner <-
        seq_len(d0 - 1) /
        d0
    
    
      Q_inner <-
        numeric(
          d0 - 1
        )
    
    
      dQ_inner <-
        matrix(
          0,
          nrow = d0 - 1,
          ncol = J
        )
    
    
      ## ----------------------------------------------------------
      ## Quantiles are increasing.
      ##
      ## We retain a robust adaptive bracket.
      ## ----------------------------------------------------------
    
      lower_start <- -8
      upper_start <-  8
    
    
      for (mm in seq_len(d0 - 1)) {
      
        p <-
          prob_inner[mm]
      
      
        root_fun <- function(x) {
        
          Sobj$cdf(x) -
            p
        }
      
      
        ## --------------------------------------------------------
        ## Use the previous quantile as a natural lower bracket
        ## when possible.
        ## --------------------------------------------------------
      
        if (mm == 1L) {
        
          lower <-
            lower_start
        
        } else {
        
          lower <-
            Q_inner[mm - 1]
        }
      
      
        upper <-
          upper_start
      
      
        ## --------------------------------------------------------
        ## Expand lower tail if necessary
        ## --------------------------------------------------------
      
        while (
          root_fun(lower) > 0
        ) {
        
          if (lower >= 0) {
            lower <- -1
          } else {
            lower <- 2 * lower
          }
        }
      
      
        ## --------------------------------------------------------
        ## Expand upper tail if necessary
        ## --------------------------------------------------------
      
        while (
          root_fun(upper) < 0
        ) {
        
          if (upper <= 0) {
            upper <- 1
          } else {
            upper <- 2 * upper
          }
        }
      
      
        ## --------------------------------------------------------
        ## Numerical inversion of an EXACT analytical CDF
        ## --------------------------------------------------------
      
        q <-
          uniroot(
            root_fun,
            lower = lower,
            upper = upper,
            tol = root_tol
          )$root
      
      
        Q_inner[mm] <-
          q
      
      
        ## --------------------------------------------------------
        ## Exact F_S(q), f_S(q), dF_S(q)/dpsi
        ## --------------------------------------------------------
      
        tmp <-
          Sobj$eval(q)
      
      
        if (!is.finite(tmp$f) ||
            tmp$f <= 0) {
          stop(
            "Invalid exact density at S quantile."
          )
        }
      
      
        ## --------------------------------------------------------
        ## Implicit differentiation:
        ##
        ## F_S(Q(p);psi) = p
        ##
        ## =>
        ##
        ## dQ/dpsi
        ## =
        ## -(dF/dpsi)/f
        ## --------------------------------------------------------
      
        dQ_inner[mm, ] <-
          -tmp$dF /
          tmp$f
      }
    
    
      Q_full <-
        c(
          -Inf,
          Q_inner,
          Inf
        )
    
    
      ## ----------------------------------------------------------
      ## Endpoint probabilities are constants 0 and 1,
      ## hence their derivatives in the conditional-cell
      ## calculations are zero.
      ## ----------------------------------------------------------
    
      dQ_full <-
        rbind(
          rep(0, J),
          dQ_inner,
          rep(0, J)
        )
    
    
      ## ----------------------------------------------------------
      ## eta_j is required ONLY for the OUTER likelihood integral.
      ##
      ## This is the only place where numerical quadrature nodes
      ## enter after F_S and f_S have been made exact.
      ## ----------------------------------------------------------
    
      eta_j <-
        as.vector(
          Znode %*%
            psi_j
        )
    
    
      list(
        Q = Q_full,
        dQ = dQ_full,
        sigma = sigma_j,
        eta = eta_j,
        exact_S_distribution = Sobj
      )
    }
  


  # ----------------------------------------------------------
  # Generic quantile + analytic derivative for arbitrary
  # probabilities. This is needed for the zero layer because
  # p_B can vary by observation.
  # ----------------------------------------------------------

  compute_prob_Q_and_dQ <- function(
      loading,
      prob
  ) {

    loading <- as.numeric(loading)
    prob <- as.numeric(prob)

    if (length(loading) != J) {
      stop("`loading` must have length q.")
    }

    if (any(!is.finite(prob)) ||
        any(prob <= 0) ||
        any(prob >= 1)) {
      stop("All probabilities must lie strictly in (0,1).")
    }

    norm2 <- sum(loading^2)

    if (!is.finite(norm2) ||
        norm2 >= 1) {
      stop("||loading|| must be smaller than 1.")
    }

    sigma_j <- sqrt(1 - norm2)

    Sobj <- make_exact_sum_sd_biexp(
      psi_j = loading
    )

    # Solve only once for repeated probability values.
    up <- unique(prob)

    Q_u <- numeric(length(up))
    dQ_u <- matrix(
      0,
      nrow = length(up),
      ncol = J
    )

    ord <- order(up)
    up_sorted <- up[ord]

    Q_sorted <- numeric(length(up_sorted))
    dQ_sorted <- matrix(
      0,
      nrow = length(up_sorted),
      ncol = J
    )

    lower_start <- -8
    upper_start <- 8

    for (ii in seq_along(up_sorted)) {

      p <- up_sorted[ii]

      root_fun <- function(x) {
        Sobj$cdf(x) - p
      }

      lower <- if (ii == 1L) {
        lower_start
      } else {
        Q_sorted[ii - 1L]
      }

      upper <- upper_start

      while (root_fun(lower) > 0) {
        lower <- lower * 2
        if (!is.finite(lower) || lower < -1e8) {
          stop("Unable to bracket BiExp quantile from below.")
        }
      }

      while (root_fun(upper) < 0) {
        upper <- upper * 2
        if (!is.finite(upper) || upper > 1e8) {
          stop("Unable to bracket BiExp quantile from above.")
        }
      }

      qv <- uniroot(
        root_fun,
        interval = c(lower, upper),
        tol = root_tol
      )$root

      ev <- Sobj$eval(qv)

      if (!is.finite(ev$f) || ev$f <= 0) {
        stop("Invalid exact density at BiExp quantile.")
      }

      Q_sorted[ii] <- qv
      dQ_sorted[ii, ] <- -ev$dF / ev$f
    }

    # Undo sorted order, then map back to original probabilities.
    Q_u[ord] <- Q_sorted
    dQ_u[ord, ] <- dQ_sorted

    idx <- match(prob, up)

    list(
      Q = Q_u[idx],
      dQ = dQ_u[idx, , drop = FALSE],
      sigma = sigma_j,
      exact_S_distribution = Sobj
    )
  }


  # Construct the outer integration grid once.
  gl_grid <- make_gl_grid(
    J = J,
    gl_order = gl_order
  )

  list(
    lambda = lambda,
    skewness = 2 * (
      lambda^3 -
        (1 - lambda^2)^(3 / 2)
    ),
    F = F_sd_biexp,
    f = f_sd_biexp,
    Q = Q_sd_biexp,
    grid = gl_grid,
    Znode = gl_grid$z,
    Wnode = gl_grid$w,
    n_node = length(gl_grid$w),
    make_exact_sum = make_exact_sum_sd_biexp,
    compute_Q_and_dQ = compute_Q_and_dQ,
    compute_prob_Q_and_dQ = compute_prob_Q_and_dQ
  )
}


# ============================================================
# .biexp_lambda_from_skewness
#
# Inverts
#
#   gamma_1(lambda)
#     = 2 [lambda^3 - (1-lambda^2)^(3/2)]
#
# for lambda in (0,1).
# ============================================================

.biexp_lambda_from_skewness <- function(
    skewness,
    tol = 1e-12
) {

  if (!is.numeric(skewness) ||
      length(skewness) != 1L ||
      !is.finite(skewness) ||
      skewness <= -2 ||
      skewness >= 2) {
    stop(
      "`skewness` must be a finite scalar strictly between -2 and 2."
    )
  }

  if (abs(skewness) < tol) {
    return(1 / sqrt(2))
  }

  skew_fun <- function(lambda) {
    2 * (
      lambda^3 -
        (1 - lambda^2)^(3 / 2)
    )
  }

  uniroot(
    function(lambda) skew_fun(lambda) - skewness,
    interval = c(1e-10, 1 - 1e-10),
    tol = tol
  )$root
}



prepare_multi_zipm_I <- function(
    Y = NULL,
    X = NULL,
    Z = NULL,
    margin = NULL,
    mixing.distribution = "gamma",
    var.power = NULL,
    margin.ini.par = NULL,
    margin.n_workers = NULL,
    margin.parallel = TRUE,
    margin.block_size = NULL,
    d0 = 10,
    progress = TRUE
) {

  .progress_msg <- function(...) {
    if (isTRUE(progress)) {
      cat(..., "\n", sep = "")
      flush.console()
    }
  }

  allowed_dist <- c(
    "gamma", "lnorm", "ig",
    "weibull", "loglogistic", "betaprime"
  )
  default_vp <- c(
    gamma = 2L, lnorm = 2L,  ig = 3L,
    weibull = 2L, loglogistic = 2L, betaprime = 2L
  )

  # ------------------------------------------------------------------
  # Route 1: reuse the best marginal models from select_multi_margin().
  # ------------------------------------------------------------------
  if (!is.null(margin)) {
    if (!inherits(margin, "multi_margin_fit")) {
      stop("`margin` must be NULL or an object returned by select_multi_margin().")
    }
    # Compact select_multi_margin() objects deliberately omit X/Z to
    # reduce RAM.  They are valid inputs when the marginal I-cache is present.
    # X/Z are required only if prepare_multi_zipm_I() must recompute I.
    required <- c("Y", "n", "k", "selected_margins")
    miss <- required[!vapply(required, function(nm) !is.null(margin[[nm]]), logical(1))]
    if (length(miss) > 0L) {
      stop("`margin` is missing: ", paste(miss, collapse = ", "), ".")
    }
    if (!is.list(margin$selected_margins) ||
        length(margin$selected_margins) != margin$k ||
        any(vapply(margin$selected_margins, is.null, logical(1)))) {
      stop("`margin$selected_margins` must contain k successful best marginal fits.")
    }

    margin_fit <- margin
    .progress_msg("\n============================================================")
    .progress_msg("[fit_multizimp] Reusing best models from select_multi_margin().")

  } else {

    # ---------------------------------------------------------------
    # Route 2: direct user input. Fit exactly ONE model per margin.
    # There is deliberately no model-selection step in fit_multizimp().
    # ---------------------------------------------------------------
    if (is.null(Y)) stop("Supply either `margin` or `Y`.")
    if (is.null(X)) stop("When `margin = NULL`, `X` must be supplied.")

    if (!(is.matrix(Y) || is.data.frame(Y))) {
      stop("`Y` must be an n x k matrix/data.frame of count responses.")
    }
    yn <- colnames(Y)
    if (is.data.frame(Y) &&
        any(!vapply(Y, function(z) is.numeric(z) || is.integer(z), logical(1)))) {
      stop("Every column of `Y` must be numeric/integer.")
    }
    Y <- as.matrix(Y)
    if (!is.numeric(Y) || any(!is.finite(Y)) || any(Y < 0) ||
        any(abs(Y - round(Y)) > 1e-8)) {
      stop("Every element of `Y` must be a finite non-negative integer count.")
    }
    storage.mode(Y) <- "integer"
    n <- nrow(Y)
    k <- ncol(Y)
    if (n < 1L || k < 2L) stop("`Y` must have n >= 1 and k >= 2.")
    if (is.null(yn)) yn <- paste0("Y", seq_len(k))
    colnames(Y) <- yn

    # A vector is interpreted as one common covariate and is shared by
    # every margin. Data-frame factors/characters are one-hot encoded.
    .design_one <- function(A, name) {
      if (is.vector(A) && !is.list(A)) {
        if (length(A) != n) stop(name, " vector must have length nrow(Y).")
        if (!is.numeric(A) && !is.integer(A) && !is.logical(A)) {
          stop(name, " vector must be numeric/integer/logical.")
        }
        A <- matrix(as.numeric(A), ncol = 1L,
                    dimnames = list(NULL, name))
      }

      if (is.data.frame(A)) {
        if (nrow(A) != n) stop(name, " must have nrow(Y) rows.")
        if (!ncol(A)) stop(name, " must contain at least one variable.")
        A[] <- lapply(A, function(z) if (is.character(z)) factor(z) else z)
        ok <- vapply(
          A,
          function(z) is.numeric(z) || is.integer(z) ||
            is.logical(z) || is.factor(z),
          logical(1)
        )
        if (any(!ok)) {
          stop(name, " contains unsupported column type(s): ",
               paste(names(A)[!ok], collapse = ", "), ".")
        }
        if (anyNA(A)) stop(name, " contains NA.")

        has_int <- any(vapply(A, function(z) {
          (is.numeric(z) || is.integer(z) || is.logical(z)) &&
            all(is.finite(as.numeric(z))) &&
            all(abs(as.numeric(z) - 1) < 1e-12)
        }, logical(1)))

        pieces <- list()
        for (nm in names(A)) {
          z <- A[[nm]]
          if (is.factor(z)) {
            if (has_int) {
              mm <- model.matrix(~ z)[, -1L, drop = FALSE]
              if (ncol(mm)) colnames(mm) <- paste0(nm, levels(z)[-1L])
            } else {
              mm <- model.matrix(~ z - 1)
              colnames(mm) <- paste0(nm, levels(z))
            }
            if (ncol(mm)) pieces[[length(pieces) + 1L]] <- mm
          } else {
            vv <- as.numeric(z)
            if (any(!is.finite(vv))) stop(name, "$", nm, " contains NA/NaN/Inf.")
            pieces[[length(pieces) + 1L]] <-
              matrix(vv, ncol = 1L, dimnames = list(NULL, nm))
          }
        }
        A <- do.call(cbind, pieces)
      }

      if (!is.matrix(A)) {
        stop(name, " must be a vector, matrix, data.frame, or (outside this helper) a list.")
      }
      if (nrow(A) != n) stop(name, " must have nrow(Y) rows.")
      if (!is.numeric(A) || any(!is.finite(A)) || ncol(A) < 1L) {
        stop(name, " must produce a finite numeric design matrix.")
      }
      storage.mode(A) <- "double"
      if (is.null(colnames(A))) colnames(A) <- paste0(name, seq_len(ncol(A)))
      A
    }

    .design_list <- function(A, name, allow_null = FALSE) {
      if (is.null(A)) {
        if (allow_null) return(rep(list(NULL), k))
        stop(name, " must not be NULL.")
      }

      # A single vector/matrix/data.frame is shared by all k margins.
      if ((is.vector(A) && !is.list(A)) || is.matrix(A) || is.data.frame(A)) {
        M <- .design_one(A, name)
        return(rep(list(M), k))
      }

      # A list supplies one design per margin.
      if (is.list(A)) {
        if (length(A) != k) stop(name, " list must have length k = ", k, ".")
        out <- vector("list", k)
        for (j in seq_len(k)) {
          if (is.null(A[[j]])) {
            if (allow_null) {
              out[[j]] <- NULL
              next
            }
            stop(name, "[[", j, "]] is NULL.")
          }
          out[[j]] <- .design_one(A[[j]], paste0(name, "[[", j, "]]"))
        }
        return(out)
      }

      stop(name, " must be a vector/matrix/data.frame or a list of length k.")
    }

    X_list <- .design_list(X, "X", allow_null = FALSE)
    Z_list <- .design_list(Z, "Z", allow_null = TRUE)

    # mixing.distribution: scalar -> recycle; otherwise exactly k values.
    if (!is.character(mixing.distribution) || !length(mixing.distribution)) {
      stop("`mixing.distribution` must be a character value or vector.")
    }
    if (length(mixing.distribution) == 1L) {
      mixing.distribution <- rep(mixing.distribution, k)
    } else if (length(mixing.distribution) != k) {
      stop("`mixing.distribution` must have length 1 or k = ", k, ".")
    }
    bad_dist <- !(mixing.distribution %in% allowed_dist)
    if (any(bad_dist)) {
      stop("Unsupported `mixing.distribution`: ",
           paste(unique(mixing.distribution[bad_dist]), collapse = ", "), ".")
    }

    # var.power: NULL -> native/default separately for each distribution;
    # scalar -> recycle; otherwise exactly k values, all in {1,2,3}.
    if (is.null(var.power)) {
      vp_fit <- rep(NA_integer_, k)
      vp_report <- unname(default_vp[mixing.distribution])
      use_custom_vp <- rep(FALSE, k)
    } else {
      if (!is.numeric(var.power) || any(!is.finite(var.power)) ||
          any(abs(var.power - round(var.power)) > 1e-12) ||
          any(!(as.integer(round(var.power)) %in% 1:3))) {
        stop("`var.power` must be NULL or contain only 1, 2, or 3.")
      }
      var.power <- as.integer(round(var.power))
      if (length(var.power) == 1L) {
        var.power <- rep(var.power, k)
      } else if (length(var.power) != k) {
        stop("`var.power` must have length 1 or k = ", k, ".")
      }
      vp_fit <- var.power
      vp_report <- var.power
      use_custom_vp <- rep(TRUE, k)

    }

    # User starts: NULL, list length k, or a single start only when k=1.
    if (!is.null(margin.ini.par) &&
        (!is.list(margin.ini.par) || length(margin.ini.par) != k)) {
      stop("`margin.ini.par` must be NULL or a list of length k.")
    }

    # Public parallel convention: margin.parallel=FALSE -> serial.
    # If TRUE, margin.n_workers controls PSOCK workers; NULL defaults to
    # min(k, detected cores - 1).
    if (!is.logical(margin.parallel) || length(margin.parallel) != 1L ||
        is.na(margin.parallel)) {
      stop("`margin.parallel` must be TRUE or FALSE.")
    }

    if (!isTRUE(margin.parallel)) {
      n_workers <- 1L
    } else if (is.null(margin.n_workers)) {
      dc <- parallel::detectCores(logical = TRUE)
      if (!is.finite(dc)) dc <- 2L
      n_workers <- max(1L, min(k, as.integer(dc) - 1L))
    } else {
      if (length(margin.n_workers) != 1L || !is.numeric(margin.n_workers) ||
          !is.finite(margin.n_workers) || margin.n_workers < 1 ||
          margin.n_workers != as.integer(margin.n_workers)) {
        stop("`margin.n_workers` must be NULL or a positive integer.")
      }
      n_workers <- min(k, as.integer(margin.n_workers))
    }

    jobs <- lapply(seq_len(k), function(j) list(
      j = j,
      dist = mixing.distribution[j],
      vp_fit = if (is.na(vp_fit[j])) NULL else vp_fit[j],
      vp_report = vp_report[j],
      use_custom_vp = use_custom_vp[j],
      start = if (is.null(margin.ini.par)) NULL else margin.ini.par[[j]]
    ))

    .fit_margin_job <- function(job, Y, X_list, Z_list) {
      j <- job$j
      t0 <- proc.time()[["elapsed"]]
      ans <- fit_zimp(
        Y = Y[, j],
        X = X_list[[j]],
        Z = Z_list[[j]],
        offset.count = NULL,
        mixing.distribution = job$dist,
        var.power = job$vp_fit,
        ini.par = job$start,
        progress = FALSE
      )

      # Normalize the margin object for the existing copula backend.
      ans$margin <- j
      ans$y <- as.integer(Y[, j])
      ans$X <- X_list[[j]]
      ans$Z <- Z_list[[j]]
      ans$response_name <- colnames(Y)[j]
      ans$beta_names <- colnames(X_list[[j]])
      ans$gamma_names <- if (is.null(Z_list[[j]])) NULL else colnames(Z_list[[j]])
      ans$distribution <- job$dist
      ans$dist <- job$dist
      ans$mixing.distribution <- job$dist
      ans$var.power <- job$vp_report
      ans$var_power <- job$vp_report
      ans$use_custom_var_power <- isTRUE(job$use_custom_vp)
      ans$zero_inflation <- !is.null(Z_list[[j]])
      if (is.null(ans$be) && !is.null(ans$beta)) ans$be <- ans$beta
      if (is.null(ans$beta) && !is.null(ans$be)) ans$beta <- ans$be
      if (is.null(ans$gam) && !is.null(ans$al)) ans$gam <- ans$al
      if (is.null(ans$al) && !is.null(ans$gam)) ans$al <- ans$gam
      ans$alpha <- ans$gam
      ans$elapsed_sec <- unname(proc.time()[["elapsed"]] - t0)
      class(ans) <- unique(c("mzimp_margin_fit", class(ans)))
      ans
    }

    .progress_msg("\n============================================================")
    .progress_msg(
      "[fit_multizimp] Fitting exactly one marginal model per response",
      " | k=", k,
      " | ", if (n_workers > 1L) paste0("parallel | ", n_workers, " workers") else "serial"
    )

    if (n_workers <= 1L) {
      M <- lapply(jobs, .fit_margin_job, Y = Y, X_list = X_list, Z_list = Z_list)
    } else {
      cl <- parallel::makeCluster(n_workers, outfile = "")
      on.exit(parallel::stopCluster(cl), add = TRUE)

      dev_pkg_path <-
        if (file.exists("DESCRIPTION") && dir.exists("R"))
          normalizePath(".", winslash = "/", mustWork = TRUE)
        else NULL

      parallel::clusterCall(cl, function(dev_pkg_path) {
        if (!is.null(dev_pkg_path)) {
          if (!requireNamespace("pkgload", quietly = TRUE))
            stop("pkgload is required on development workers.")
          pkgload::load_all(
            dev_pkg_path, quiet = TRUE, recompile = FALSE, export_all = TRUE
          )
          return(invisible("development"))
        }
        if (requireNamespace("mzimp", quietly = TRUE)) {
          suppressPackageStartupMessages(library(mzimp))
          return(invisible("installed"))
        }
        stop("mzimp is unavailable on workers.")
      }, dev_pkg_path)

      M <- parallel::parLapplyLB(
        cl, jobs, .fit_margin_job,
        Y = Y, X_list = X_list, Z_list = Z_list
      )
    }

    bad <- which(vapply(M, function(x) {
      ll <- x$log_lik_value
      is.null(ll) || length(ll) != 1L || !is.finite(ll)
    }, logical(1)))
    if (length(bad)) {
      stop("Direct marginal fitting failed for margin(s): ",
           paste(bad, collapse = ", "), ".")
    }

    selected_summary <- data.frame(
      margin = seq_len(k),
      distribution = mixing.distribution,
      var.power = vp_report,
      zero_inflation = vapply(Z_list, Negate(is.null), logical(1)),
      logLik = vapply(M, function(x) x$log_lik_value, numeric(1)),
      stringsAsFactors = FALSE
    )

    margin_fit <- list(
      Y = Y,
      X = X_list,
      Z = Z_list,
      n = n,
      k = k,
      selected_margins = M,
      selected_summary = selected_summary,
      mixing.distribution = mixing.distribution,
      var.power = if (is.null(var.power)) NULL else var.power,
      I_precomputed = FALSE,
      I_d0 = NULL,
      I_y_list = NULL,
      I_0_list = NULL,
      pb_list = NULL,
      y0_list = NULL,
      response_names = colnames(Y),
      direct_fit = TRUE
    )
    class(margin_fit) <- "multi_margin_fit"
  }

  Y <- as.matrix(margin_fit$Y)
  M <- margin_fit$selected_margins
  n <- nrow(Y)
  k <- ncol(Y)

  # ---------------------------------------------------------------
  # Normalize compact V6 X/Z storage to the legacy per-margin view
  # expected by the copula backend.  V6 stores a common matrix only
  # once and optionally stores per-margin column indices; true list
  # inputs remain lists.  rep(list(M), k) only creates a lightweight
  # list of references to the same common matrix (no eager k-fold copy).
  # ---------------------------------------------------------------
  .expand_margin_design <- function(data, is_list = FALSE, columns = NULL,
                                    k, name, allow_null = FALSE) {
    if (is.null(data)) {
      if (allow_null) return(rep(list(NULL), k))
      stop("margin$", name, " is missing.")
    }

    if (isTRUE(is_list)) {
      if (!is.list(data) || length(data) != k)
        stop("margin$", name, " must be a list of length k when `",
             name, "_is_list = TRUE`.")
      return(data)
    }

    # Backward compatibility: old select objects already stored a list.
    if (is.list(data)) {
      if (length(data) != k)
        stop("margin$", name, " must be a list of length k.")
      return(data)
    }

    M0 <- as.matrix(data)
    if (nrow(M0) != n)
      stop("margin$", name, " must have nrow(Y) rows.")

    if (is.null(columns)) {
      return(rep(list(M0), k))
    }

    if (!is.list(columns) || length(columns) != k)
      stop("margin$", name, "_columns must be NULL or a list of length k.")

    lapply(seq_len(k), function(j) {
      jj <- columns[[j]]
      if (is.null(jj)) return(M0)
      M0[, jj, drop = FALSE]
    })
  }

  X_list <- .expand_margin_design(
    data = margin_fit$X,
    is_list = isTRUE(margin_fit$X_is_list),
    columns = margin_fit$X_columns,
    k = k, name = "X", allow_null = FALSE
  )

  Z_list <- .expand_margin_design(
    data = margin_fit$Z,
    is_list = isTRUE(margin_fit$Z_is_list),
    columns = margin_fit$Z_columns,
    k = k, name = "Z", allow_null = TRUE
  )

  if (k < 2L) stop("prepare_multi_zipm_I() requires k >= 2.")
  d0 <- as.integer(d0)
  if (d0 < 2L) stop("d0 must be >= 2.")

  get_I_vec <- function(mj, z, d, Xj) {
    dist <- as.character(mj$distribution)[1L]
    Xj <- as.matrix(Xj)
    be <- mj$be
    phi <- mj$phi
    arg_vp <- if (isTRUE(mj$use_custom_var_power)) list(var_power = mj$var_power) else list()
    base_args <- list(z = z, X = Xj, be = be, phi = phi, d = d)

    if (dist == "gamma") {
      out <- do.call(I_gamma, c(base_args, arg_vp))
    } else if (dist == "lnorm") {
      out <- do.call(I_lnorm, c(base_args, list(n_gh = 20, gh = NULL), arg_vp))
    } else if (dist == "ig") {
      out <- do.call(I_ig, c(base_args, arg_vp))
    } else if (dist == "weibull") {
      out <- do.call(I_weibull, c(base_args, arg_vp))
    } else if (dist == "loglogistic") {
      if (!exists("I_loglogistic", mode = "function", inherits = TRUE))
        stop("Selected loglogistic margin requires I_loglogistic().")
      out <- do.call(get("I_loglogistic", mode = "function", inherits = TRUE), c(base_args, arg_vp))
    } else if (dist == "betaprime") {
      if (!exists("I_betaprime", mode = "function", inherits = TRUE))
        stop("Selected betaprime margin requires I_betaprime().")
      out <- do.call(get("I_betaprime", mode = "function", inherits = TRUE), c(base_args, arg_vp))
    } else {
      stop("Unsupported selected marginal distribution: ", dist)
    }
    as.matrix(out)
  }

  marginal.dist <- vapply(M, function(x) as.character(x$distribution)[1L], character(1))
  zero_inflation <- vapply(M, function(x) isTRUE(x$zero_inflation), logical(1))
  zi_index <- which(zero_inflation)
  n_zi <- length(zi_index)

  I_y_list <- vector("list", k)
  I_0_list <- vector("list", k)
  pb_list <- vector("list", k)
  y0_list <- vector("list", k)
  zi_position <- integer(k)
  if (n_zi > 0L) zi_position[zi_index] <- seq_len(n_zi)

  cache_ok <-
    !is.null(margin) &&
    isTRUE(margin_fit$I_precomputed) &&
    !is.null(margin_fit$I_d0) &&
    identical(as.integer(margin_fit$I_d0), as.integer(d0)) &&
    is.list(margin_fit$I_y_list) && length(margin_fit$I_y_list) == k &&
    is.list(margin_fit$I_0_list) && length(margin_fit$I_0_list) == k &&
    is.list(margin_fit$pb_list) && length(margin_fit$pb_list) == k &&
    is.list(margin_fit$y0_list) && length(margin_fit$y0_list) == k

  if (isTRUE(cache_ok)) {
    dim_ok <- TRUE
    for (j in seq_len(k)) {
      Iy <- margin_fit$I_y_list[[j]]
      I0 <- margin_fit$I_0_list[[j]]
      if (!is.matrix(Iy) || !is.matrix(I0) ||
          nrow(Iy) != n || ncol(Iy) != d0 ||
          nrow(I0) != n || ncol(I0) != d0 ||
          length(margin_fit$pb_list[[j]]) != n ||
          (!is.null(margin_fit$y0_list[[j]]) &&
           length(margin_fit$y0_list[[j]]) != n)) {
        dim_ok <- FALSE
        break
      }
    }
    cache_ok <- dim_ok
  }

  if (isTRUE(cache_ok)) {
    I_y_list <- margin_fit$I_y_list
    I_0_list <- margin_fit$I_0_list
    pb_list <- margin_fit$pb_list
    y0_list <- margin_fit$y0_list

    # The fused GH L=1 C++ kernel expects a logical y0 vector for EVERY
    # margin, including PM margins.  Compact marginal objects may store
    # y0_list[[j]] = NULL for PM margins because it is unnecessary for
    # marginal selection.  Reconstruct this cheap indicator from Y here;
    # this does not require X/Z and costs only n*k logical values.
    for (j in seq_len(k)) {
      if (is.null(y0_list[[j]])) {
        y0_list[[j]] <- as.logical(Y[, j] == 0L)
      } else {
        y0_list[[j]] <- as.logical(y0_list[[j]])
      }
      if (is.null(pb_list[[j]])) {
        if (zero_inflation[j]) {
          stop("Precomputed cache is missing `pb_list[[", j,
               "]]` for a ZIMP margin.")
        }
        pb_list[[j]] <- rep(1, n)
      }
    }

    .progress_msg("Marginal interval-probability matrices (I) are already available in `margin`; reusing them.")
  } else {
    # A compact margin object without X/Z can only be used through its
    # precomputed I-cache.  Do not silently attempt to reconstruct designs.
    if (!is.null(margin) &&
        any(vapply(X_list, is.null, logical(1)))) {    
      stop(
        "`margin` does not contain X/Z and its precomputed I-cache cannot be reused. ",
        "Use select_multi_margin(..., compute.I = TRUE) with the same d0 as ",
        "fit_multizimp(), or supply a non-compact margin object."
      )
    }

    .progress_msg("Computing marginal interval-probability matrices (I)...")
    for (j in seq_len(k)) {
      mj <- M[[j]]
      Xj <- X_list[[j]]
      Zj <- Z_list[[j]]

      I_y_list[[j]] <- as.matrix(
        get_I_vec(mj, Y[, j], d0, Xj)
      )

      I_0_list[[j]] <- as.matrix(
        get_I_vec(mj, rep(0L, n), d0, Xj)
      )

     if (zero_inflation[j]) {
        gam <- if (!is.null(mj$gam)) mj$gam else mj$al
        if (is.null(Zj) || is.null(gam))
          stop(
            "ZIMP margin ", j,
            " is missing Z or zero-inflation coefficients."
          )
        pb_list[[j]] <- plogis(as.vector(as.matrix(Zj) %*% gam))
        y0_list[[j]] <- Y[, j] == 0
      } else {
        pb_list[[j]] <- rep(1, n)
        # Keep the downstream invariant: y0_list contains one logical
        # vector of length n for every margin, PM as well as ZIMP.
        y0_list[[j]] <- as.logical(Y[, j] == 0L)
      }

      # Per-margin computational details are intentionally not printed.
    }
  }

  out <- list(
    margin_fit = margin_fit,
    Y = Y,
    X = X_list,
    Z = Z_list,
    selected_margins = M,
    n = n,
    k = k,
    d0 = d0,
    marginal_dist = marginal.dist,
    zero_inflation = zero_inflation,
    zi_index = zi_index,
    n_zi = n_zi,
    zi_position = zi_position,
    I_y_list = I_y_list,
    I_0_list = I_0_list,
    pb_list = pb_list,
    y0_list = y0_list,
    I_reused_from_margin = isTRUE(cache_ok),
    call = match.call()
  )
  class(out) <- c("multi_zipm_I", class(out))
  .progress_msg("Marginal quantities are ready.")
  .progress_msg("============================================================\n")
  out
}


# ============================================================
# .subset_multi_zipm_I
#
# Internal helper for progressive warm-start stages. It subsets
# observation-level quantities already computed by prepare_multi_zipm_I(),
# so marginal interval probabilities I are NEVER recomputed.
# ============================================================
.subset_multi_zipm_I <- function(prepared, idx) {
  if (!inherits(prepared, "multi_zipm_I")) {
    stop("`prepared` must be a multi_zipm_I object.")
  }

  idx <- sort(unique(as.integer(idx)))
  if (length(idx) < 1L || any(!is.finite(idx)) ||
      min(idx) < 1L || max(idx) > prepared$n) {
    stop("Invalid observation indices in `.subset_multi_zipm_I()`.")
  }

  # Build only the fields required by the dependence backend.  Do not start
  # from `out <- prepared`, because that keeps the full-data object graph alive
  # during progressive fitting.
  out <- list(
    margin_fit = prepared$margin_fit,
    Y = prepared$Y[idx, , drop = FALSE],
    X = lapply(prepared$X, function(x) {
      if (is.null(x)) return(NULL)
      as.matrix(x)[idx, , drop = FALSE]
    }),
    Z = lapply(prepared$Z, function(z) {
      if (is.null(z)) return(NULL)
      as.matrix(z)[idx, , drop = FALSE]
    }),
    selected_margins = prepared$selected_margins,
    n = length(idx),
    k = prepared$k,
    d0 = prepared$d0,
    marginal_dist = prepared$marginal_dist,
    zero_inflation = prepared$zero_inflation,
    zi_index = prepared$zi_index,
    n_zi = prepared$n_zi,
    zi_position = prepared$zi_position,
    I_y_list = lapply(prepared$I_y_list, function(x) x[idx, , drop = FALSE]),
    I_0_list = lapply(prepared$I_0_list, function(x) x[idx, , drop = FALSE]),
    pb_list = lapply(prepared$pb_list, function(x) x[idx]),
    y0_list = lapply(prepared$y0_list, function(x) {
      if (is.null(x)) NULL else x[idx]
    }),
    I_reused_from_margin = prepared$I_reused_from_margin,
    sample_index = idx,
    call = prepared$call
  )
  class(out) <- class(prepared)
  out
}

.fit_multizimp_copula <- function(
    prepared,
    copula.block_size = 500L,
    d0 = NULL,
    J = 1,
    L = 0,

    latent.dependence = c("gaussian", "biexp"),
    latent.factor.skewness = 0.5,
    latent.loglike = NULL,
    gh_order_N = 9,
    gh.threads = 1L,
    approx.hermite.gh_order = 21L,
    approx.hermite.joint.gh_order = 7L,
    approx.hermite.joint.order = 4L,
    approx.hermite.threads = 1L,

    # Hermite L=0 optimizer.  "polar" alternates V|r and r|V;
    # "joint" retains the original direct-Psi nlminb optimizer.
    hermite.optimizer = c("polar", "joint"),
    hermite.polar.cycles = 6L,
    hermite.polar.v.iter = 8L,
    hermite.polar.r.iter = 4L,
    hermite.polar.tol = 1e-5,
    hermite.polar.r.max = 0.995,
    hermite.polar.joint.refine.iter = 5L,

    qmc.Q = 2^12,
    qmc.n_rep = 4L,
    qmc.seed = 1L,
    qmc.chunk_size = 4096L,
    final.simulation = c("prior_qmc", "genz_conditional_qmc", "truncated_antithetic_qmc", "tilted_group_qmc", "both"),
    genz.conditional.Q = 256L,
    genz.conditional.antithetic = TRUE,
    simulation.rotation.energy = 0.99,
    simulation.eigen.tol = 1e-4,
    simulation.min.rank = 1L,
    simulation.max.rank = NULL,
    tilted.groups = 8L,
    tilted.cov.scale = 1.25,
    tilted.between.scale = 1.0,
    tilted.min.eigen = 1e-5,

    zero.dependence = c("gaussian", "biexp"),
    zero.factor.skewness = 0.5,

    gh_order_B = 9,
    gl_order_N = 15,
    gl_order_B = 15,
    biexp.root_tol = 1e-10,
    psi0_start = NULL,
    rho0_start = NULL,
    gh.progressive = TRUE,
    progressive.sample.frac = NULL,
    progressive.sample.n = NULL,
    progressive.sample.seed = 1L,
    progressive.target.iter = 10L,
    progressive.target.ladder = TRUE,
    progressive.target.ladder.frac = c(0.25, 0.50),
    progressive.target.ladder.iter = 5L,
    .gh_progressive_internal = FALSE,
    .gh_node_progressive_internal = FALSE,
    .skip_final_qmc = FALSE,
    use.gradient = TRUE,
    progress = TRUE,
    trace = 0,
    eval.max = 1000,
    iter.max = 100,
    tiny = 1e-300,
    big_penalty = 1e12,
    rel.tol = 1e-6
) {

  .progress_msg <- function(...) {
    if (isTRUE(progress)) {
      cat(..., "\n", sep = "")
      flush.console()
    }
  }

  # Implementation-level diagnostics are intentionally hidden from users.
  .internal_msg <- function(...) invisible(NULL)

  ## Clean branch policy: Hermite has one optimizer (polar) and one final
  ## likelihood evaluator (Genz conditional RQMC). Legacy alternatives remain
  ## private below only to avoid destabilizing validated L=0 code paths.
  final.simulation <- "genz_conditional_qmc"
  hermite.optimizer <- "polar"
  if (!is.finite(simulation.rotation.energy) || simulation.rotation.energy <= 0 || simulation.rotation.energy > 1)
    stop("`simulation.rotation.energy` must lie in (0,1].")
  if (!is.finite(simulation.eigen.tol) || simulation.eigen.tol < 0)
    stop("`simulation.eigen.tol` must be non-negative.")
  simulation.min.rank <- max(1L, as.integer(round(simulation.min.rank)))
  if (!is.null(simulation.max.rank))
    simulation.max.rank <- max(1L, as.integer(round(simulation.max.rank)))
  hermite.optimizer <- match.arg(hermite.optimizer)
  hermite.polar.cycles <- as.integer(round(hermite.polar.cycles))
  hermite.polar.v.iter <- as.integer(round(hermite.polar.v.iter))
  hermite.polar.r.iter <- as.integer(round(hermite.polar.r.iter))
  hermite.polar.joint.refine.iter <- 0L
  if (any(c(hermite.polar.cycles, hermite.polar.v.iter, hermite.polar.r.iter,
            hermite.polar.joint.refine.iter) < 0L)) {
    stop("Hermite polar iteration controls must be non-negative integers.")
  }
  if (!is.finite(hermite.polar.tol) || hermite.polar.tol <= 0) {
    stop("`hermite.polar.tol` must be positive.")
  }
  if (!is.finite(hermite.polar.r.max) || hermite.polar.r.max <= 0 ||
      hermite.polar.r.max >= 1) {
    stop("`hermite.polar.r.max` must lie strictly in (0,1).")
  }
  tilted.groups <- max(1L, as.integer(round(tilted.groups)))
  if (!is.finite(tilted.cov.scale) || tilted.cov.scale <= 0) {
    stop("`tilted.cov.scale` must be positive.")
  }
  if (!is.finite(tilted.between.scale) || tilted.between.scale < 0) {
    stop("`tilted.between.scale` must be non-negative.")
  }
  if (!is.finite(tilted.min.eigen) || tilted.min.eigen <= 0) {
    stop("`tilted.min.eigen` must be positive.")
  }

  ## =========================================================
  ## 1. Reuse precomputed marginal quantities
  ## =========================================================

  if (!inherits(prepared, "multi_zipm_I")) {
    stop("`prepared` must be an object returned by prepare_multi_zipm_I().")
  }

  margin_fit <- prepared$margin_fit
  Y <- prepared$Y
  X_list <- prepared$X
  Z_list <- prepared$Z
  M <- prepared$selected_margins
  n <- prepared$n
  k <- prepared$k

  if (!is.null(d0) && as.integer(d0) != prepared$d0) {
    stop("`d0` must match the value used in prepare_multi_zipm_I().")
  }
  d0 <- prepared$d0

  .internal_msg("\n============================================================")

  if (is.null(copula.block_size)) {
    copula.block_size <- n
  }

  if (!is.numeric(copula.block_size) ||
      length(copula.block_size) != 1L ||
      !is.finite(copula.block_size) ||
      copula.block_size < 1L ||
      abs(copula.block_size - round(copula.block_size)) > 1e-12) {
    stop("`copula.block_size` must be NULL or a positive integer.")
  }

  copula.block_size <- min(
    n,
    as.integer(round(copula.block_size))
  )

  if (k < 2L) stop(".fit_multizimp_copula() requires k >= 2.")
  if (!is.list(X_list) || length(X_list) != k) stop("margin$X must be a list of length k.")
  if (!is.list(Z_list) || length(Z_list) != k) stop("margin$Z must be a list of length k.")

  d0 <- as.integer(d0)
  J <- as.integer(J)

  if (d0 < 2L) stop("d0 must be >= 2.")
  if (J < 1L || J >= k) stop("J must satisfy 1 <= J < k.")

  if (!is.numeric(L) ||
      length(L) != 1L ||
      !is.finite(L) ||
      L < 0L ||
      abs(L - round(L)) > 1e-12) {
    stop("L must be a nonnegative integer.")
  }

  L <- as.integer(round(L))

  # The structural-zero factor is intentionally restricted to L = 0 or 1.
  # Count-factor dimension J remains unrestricted by this choice.
  if (!(L %in% c(0L, 1L))) {
    stop("L must be 0 or 1. The zero layer is implemented as independent (L=0) or one-factor (L=1).")
  }

  latent.dependence <- match.arg(latent.dependence)
  zero.dependence <- match.arg(zero.dependence)

  ## Gaussian count-factor likelihood method.
  ## Only two methods are supported:
  ##   GH      : raw tensor Gauss-Hermite quadrature
  ##   hermite : Hermite-polynomial approximation
  ## Default is selected from J, but expert users may override it.
  if (latent.dependence == "gaussian") {
    if (is.null(latent.loglike)) {
      latent.loglike <- if (J <= 2L) "GH" else "hermite"
    } else {
      if (!is.character(latent.loglike) || length(latent.loglike) != 1L ||
          is.na(latent.loglike)) {
        stop("`latent.loglike` must be NULL, 'GH', or 'hermite'.")
      }
      if (tolower(latent.loglike) == "gh") latent.loglike <- "GH"
      if (tolower(latent.loglike) %in% c("hermite", "hermite polynomial", "hermite polynom")) {
        latent.loglike <- "hermite"
      }
      if (!(latent.loglike %in% c("GH", "hermite"))) {
        stop("`latent.loglike` must be one of 'GH' or 'hermite'.")
      }
    }
  } else {
    if (!is.null(latent.loglike)) {
      warning("`latent.loglike` is ignored when `latent.dependence != 'gaussian'`.")
    }
    latent.loglike <- "quadrature"
  }

  ## GH runs directly without interactive warning/confirmation.

  gh_order_N <- as.integer(round(gh_order_N))
  if (!is.finite(gh_order_N) || gh_order_N < 2L) {
    stop("`gh_order_N` must be an integer >= 2.")
  }

  ## Exact Gaussian GH with one structural-zero factor uses the SAME
  ## one-dimensional quadrature order for the count and zero layers.
  ## gh_order_B is deliberately ignored in this branch.
  if (identical(latent.dependence, "gaussian") &&
      identical(latent.loglike, "GH") &&
      L == 1L &&
      identical(zero.dependence, "gaussian")) {
    gh_order_B <- gh_order_N
  }

  approx.hermite.joint.gh_order <- as.integer(round(approx.hermite.joint.gh_order))
  approx.hermite.joint.order <- as.integer(round(approx.hermite.joint.order))
  if (!is.finite(approx.hermite.joint.gh_order) || approx.hermite.joint.gh_order < 2L) {
    stop("`approx.hermite.joint.gh_order` must be an integer >= 2.")
  }
  if (!approx.hermite.joint.order %in% c(2L, 4L)) {
    stop("`approx.hermite.joint.order` must currently be 2 or 4.")
  }

  qmc.Q <- as.integer(round(qmc.Q))
  qmc.n_rep <- as.integer(round(qmc.n_rep))
  qmc.seed <- as.integer(round(qmc.seed))
  qmc.chunk_size <- as.integer(round(qmc.chunk_size))
  if (qmc.Q < 2L || qmc.n_rep < 1L || qmc.chunk_size < 1L) {
    stop("`qmc.Q`, `qmc.n_rep`, and `qmc.chunk_size` must be positive integers (qmc.Q >= 2).")
  }

  if (!is.numeric(latent.factor.skewness) ||
      length(latent.factor.skewness) != 1L ||
      !is.finite(latent.factor.skewness) ||
      latent.factor.skewness <= -2 ||
      latent.factor.skewness >= 2) {
    stop("latent.factor.skewness must lie strictly between -2 and 2.")
  }

  if (!is.numeric(zero.factor.skewness) ||
      length(zero.factor.skewness) != 1L ||
      !is.finite(zero.factor.skewness) ||
      zero.factor.skewness <= -2 ||
      zero.factor.skewness >= 2) {
    stop("zero.factor.skewness must lie strictly between -2 and 2.")
  }

  latent.factor.lambda <- if (latent.dependence == "biexp") {
    .biexp_lambda_from_skewness(latent.factor.skewness)
  } else {
    NA_real_
  }

  zero.factor.lambda <- if (
    zero.dependence == "biexp" &&
    L > 0L
  ) {
    .biexp_lambda_from_skewness(zero.factor.skewness)
  } else {
    NA_real_
  }

  marginal.dist <- vapply(
    M,
    function(x) as.character(x$distribution)[1L],
    character(1)
  )

  zero_inflation <- vapply(
    M,
    function(x) isTRUE(x$zero_inflation),
    logical(1)
  )

  zi_index <- which(zero_inflation)
  n_zi <- length(zi_index)

  # L = 0:
  #   structural-zero indicators are independent.
  #
  # L >= 1:
  #   estimate an L-factor Gaussian copula for the zero layer,
  #   provided at least two ZIPM margins exist.
  estimate_B_factor <-
    (L > 0L) &&
    n_zi >= 2L

  npar_psi <- k * J - J * (J - 1L) / 2L

  if (estimate_B_factor) {

    if (L >= n_zi) {
      stop(
        "When L > 0 and zero dependence is estimated, ",
        "L must be smaller than the number of zero-inflated margins; ",
        "n_zi = ", n_zi, "."
      )
    }

    npar_rho <- n_zi * L - L * (L - 1L) / 2L

  } else {

    npar_rho <- 0L
  }

  if (!isTRUE(.gh_progressive_internal) && !isTRUE(.gh_node_progressive_internal)) {
    .progress_msg("\nStarting factor-copula estimation")
    .progress_msg("  Observations: ", format(n, big.mark = ",", scientific = FALSE))
    .progress_msg("  Count variables: ", k)
    .progress_msg("  Mixed-Poisson margins: ", k - n_zi)
    .progress_msg("  Zero-inflated mixed-Poisson margins: ", n_zi)
    .progress_msg("  Count dependence: ", J, " latent factor", if (J == 1L) "" else "s",
                  " (", if (latent.dependence == "gaussian") "Gaussian factor copula" else paste0(latent.dependence, " factor copula"), ")")
    if (n_zi == 0L) {
      .progress_msg("  Structural-zero layer: not present")
    } else if (!estimate_B_factor) {
      .progress_msg("  Structural-zero dependence: independent")
    } else {
      .progress_msg("  Structural-zero dependence: ", L, " latent factor", if (L == 1L) "" else "s",
                    " (", if (zero.dependence == "gaussian") "Gaussian factor copula" else paste0(zero.dependence, " factor copula"), ")")
    }
    if (identical(latent.dependence, "gaussian") && identical(latent.loglike, "hermite")) {
      .progress_msg("  Likelihood used for fitting: Hermite lower-bound approximation")
    } else if (identical(latent.dependence, "gaussian") && identical(latent.loglike, "GH")) {
      .progress_msg("  Likelihood used for fitting: Gauss-Hermite quadrature")
    } else {
      .progress_msg("  Likelihood used for fitting: numerical quadrature")
    }
    .progress_msg("  Parallel workers/threads: ", gh.threads)
    .progress_msg("Constructing initial parameters...")
  }

  progressive_used <- FALSE
  progressive_history <- list()
  progressive_path <- J
  progressive_sample_index <- NULL

  ## Progressive factor-dimension initialization is only a warm start.
  ## Every predecessor J (1, 2, 4, ...) gets at most 10 optimizer iterations;
  ## only the requested final J is optimized with the user's full iter.max.
  progressive.warm.iter <- 10L

  progressive.target.iter <- as.integer(round(progressive.target.iter))
  if (!is.finite(progressive.target.iter) || progressive.target.iter < 0L) {
    stop("`progressive.target.iter` must be a non-negative integer.")
  }

  progressive.target.ladder.iter <- as.integer(round(progressive.target.ladder.iter))
  if (!is.finite(progressive.target.ladder.iter) || progressive.target.ladder.iter < 0L) {
    stop("`progressive.target.ladder.iter` must be a non-negative integer.")
  }
  if (!is.logical(progressive.target.ladder) || length(progressive.target.ladder) != 1L ||
      is.na(progressive.target.ladder)) {
    stop("`progressive.target.ladder` must be TRUE or FALSE.")
  }
  progressive.target.ladder.frac <- as.numeric(progressive.target.ladder.frac)
  # An empty ladder is valid when progressive.target.ladder = FALSE.
  # This is used by fit_multizimp(), where every GH rank stage is
  # fitted directly on the full sample (no observation-size ladder).
  if (isTRUE(progressive.target.ladder)) {
    if (length(progressive.target.ladder.frac) < 1L ||
        any(!is.finite(progressive.target.ladder.frac)) ||
        any(progressive.target.ladder.frac <= 0 | progressive.target.ladder.frac >= 1)) {
      stop("`progressive.target.ladder.frac` must contain fractions strictly in (0,1).")
    }
    progressive.target.ladder.frac <- sort(unique(progressive.target.ladder.frac))
  } else {
    progressive.target.ladder.frac <- numeric(0)
  }

  ## =========================================================
  ## Progressive warm start for Gaussian count factors
  ##
  ## GH paths (raw tensor GH; J <= 5 only):
  ##   J=1 : 1
  ##   J=2 : 1 -> 2
  ##   J=3 : 1 -> 3
  ##   J=4 : 1 -> 2 -> 4
  ##   J=5 : 1 -> 3 -> 5
  ##
  ## Hermite paths use powers of two, with the exact target appended:
  ##   1, 2, 4, 8, ... , J
  ## e.g. J=3: 1->2->3; J=5: 1->2->4->5; J=10: 1->2->4->8->10.
  ##
  ## All predecessor stages may use ONE fixed subsample. The final
  ## target-J stage always uses 100% of the prepared observations.
  ## No I matrix is recomputed: subsampling slices the precomputed I.
  ## =========================================================

  make_progressive_path <- function(J, method) {
    J <- as.integer(J)
    if (method == "GH") {
      if (J > 5L) stop("Progressive GH supports J <= 5 only.")
      return(switch(
        as.character(J),
        `1` = 1L,
        `2` = c(1L, 2L),
        `3` = c(1L, 3L),
        `4` = c(1L, 2L, 4L),
        `5` = c(1L, 3L, 5L)
      ))
    }

    # Hermite: powers of two below J, then exact target J.
    pp <- 1L
    while (tail(pp, 1L) * 2L < J) {
      pp <- c(pp, tail(pp, 1L) * 2L)
    }
    if (tail(pp, 1L) != J) pp <- c(pp, J)
    unique(pp)
  }

  make_progressive_sample <- function(n) {
    if (!is.null(progressive.sample.n) && !is.null(progressive.sample.frac)) {
      stop("Specify at most one of `progressive.sample.n` and `progressive.sample.frac`.")
    }

    n_stage <- n
    if (!is.null(progressive.sample.n)) {
      if (!is.numeric(progressive.sample.n) || length(progressive.sample.n) != 1L ||
          !is.finite(progressive.sample.n) || progressive.sample.n < 1L) {
        stop("`progressive.sample.n` must be NULL or a positive integer.")
      }
      n_stage <- min(n, as.integer(round(progressive.sample.n)))
    } else if (!is.null(progressive.sample.frac)) {
      if (!is.numeric(progressive.sample.frac) || length(progressive.sample.frac) != 1L ||
          !is.finite(progressive.sample.frac) || progressive.sample.frac <= 0 ||
          progressive.sample.frac > 1) {
        stop("`progressive.sample.frac` must lie in (0,1].")
      }
      n_stage <- max(1L, min(n, as.integer(round(n * progressive.sample.frac))))
    }

    if (n_stage >= n) return(seq_len(n))
    set.seed(as.integer(progressive.sample.seed))
    sort(sample.int(n, size = n_stage, replace = FALSE))
  }

  if (latent.dependence == "gaussian" &&
      isTRUE(gh.progressive) &&
      !isTRUE(.gh_progressive_internal) &&
      is.null(psi0_start) &&
      J > 1L) {

    progressive_path <- make_progressive_path(J, latent.loglike)
    predecessor_J <- head(progressive_path, -1L)

    if (length(predecessor_J) > 0L) {
      progressive_used <- TRUE
      progressive_sample_index <- make_progressive_sample(n)
      stage_prepared <- if (length(progressive_sample_index) < n) {
        .subset_multi_zipm_I(prepared, progressive_sample_index)
      } else {
        prepared
      }

      .internal_msg(
        "[dependence] progressive J: ", paste(progressive_path, collapse = " -> "),
        " | full-data warm n=", stage_prepared$n, "/", n
      )

      prev_fit <- NULL

      for (J_stage in predecessor_J) {
        stage_psi0 <- if (is.null(prev_fit)) {
          NULL
        } else {
          upgrade_factor_start(
            Psi_old = prev_fit$count_layer$psi,
            q_new = J_stage,
            loading0 = 0.05,
            offdiag0 = 0.02
          )
        }

        stage_rho0 <- rho0_start
        if (!is.null(prev_fit) &&
            isTRUE(prev_fit$zero_layer$dependence_estimated) &&
            !is.null(prev_fit$zero_layer$rho0)) {
          stage_rho0 <- prev_fit$zero_layer$rho0
        }



        prev_fit <- .fit_multizimp_copula(
          prepared = stage_prepared,
          copula.block_size = min(copula.block_size, stage_prepared$n),
          d0 = d0,
          J = J_stage,
          L = L,
          latent.dependence = latent.dependence,
          latent.factor.skewness = latent.factor.skewness,
          latent.loglike = latent.loglike,
          gh_order_N = gh_order_N,
          gh.threads = gh.threads,
          ## Hermite predecessor-J stages are only warm starts:
          ## use the cheapest Q=2 grid here.  The requested final J
          ## subsequently follows the progressive-Q ladder below.
          approx.hermite.gh_order =
            if (identical(latent.loglike, "hermite")) 2L else approx.hermite.gh_order,
          approx.hermite.joint.gh_order =
            if (identical(latent.loglike, "hermite")) 2L else approx.hermite.joint.gh_order,
          approx.hermite.joint.order = approx.hermite.joint.order,
          approx.hermite.threads = approx.hermite.threads,
          qmc.Q = qmc.Q,
          qmc.n_rep = qmc.n_rep,
          qmc.seed = qmc.seed,
          qmc.chunk_size = qmc.chunk_size,
          zero.dependence = zero.dependence,
          zero.factor.skewness = zero.factor.skewness,
          gh_order_B = gh_order_B,
          gl_order_N = gl_order_N,
          gl_order_B = gl_order_B,
          biexp.root_tol = biexp.root_tol,
          psi0_start = stage_psi0,
          rho0_start = stage_rho0,
          ## Warm predecessor-J fits are deliberately short.  Do not start a
          ## second progressive quadrature ladder inside a warm J stage: use
          ## the requested quadrature order directly and spend at most
          ## `progressive.warm.iter` optimizer iterations at this J.
          gh.progressive = FALSE,
          progressive.sample.frac = NULL,
          progressive.sample.n = NULL,
          progressive.sample.seed = progressive.sample.seed,
          progressive.target.iter = progressive.target.iter,
          .gh_progressive_internal = TRUE,
          .skip_final_qmc = TRUE,
          use.gradient = use.gradient,
          progress = FALSE,
          trace = trace,
          eval.max = max(30L, 5L * progressive.warm.iter),
          iter.max = progressive.warm.iter,
          tiny = tiny,
          big_penalty = big_penalty,
          rel.tol = rel.tol
        )

        .internal_msg(
          "[dependence] warm J=", J_stage,
          " | n=", stage_prepared$n,
          " | logLik=", signif(prev_fit$copula_fit$log_lik_value, 10),
          " | iter=", if (!is.null(prev_fit$copula_fit$optim$iterations)) prev_fit$copula_fit$optim$iterations else NA_integer_,
          " | elapsed=", signif(prev_fit$copula_fit$elapsed_sec, 5), " sec"
        )

        progressive_history[[length(progressive_history) + 1L]] <- list(
          J = J_stage,
          method = latent.loglike,
          n = stage_prepared$n,
          logLik = prev_fit$copula_fit$log_lik_value,
          convergence = prev_fit$copula_fit$convergence,
          elapsed_sec = prev_fit$copula_fit$elapsed_sec
        )
      }

      psi0_start <- upgrade_factor_start(
        Psi_old = prev_fit$count_layer$psi,
        q_new = J,
        loading0 = 0.05,
        offdiag0 = 0.02
      )

      if (isTRUE(prev_fit$zero_layer$dependence_estimated) &&
          !is.null(prev_fit$zero_layer$rho0)) {
        rho0_start <- prev_fit$zero_layer$rho0
      }

      ## Before switching to the full sample, allow a few iterations at
      ## the FINAL J on the same progressive subsample.  This avoids an
      ## abrupt jump from (J_previous, subsample) to (J_target, full data).
      if (stage_prepared$n < n && progressive.target.iter > 0L) {
        target_sample_fit <- .fit_multizimp_copula(
          prepared = stage_prepared,
          copula.block_size = min(copula.block_size, stage_prepared$n),
          d0 = d0, J = J, L = L,
          latent.dependence = latent.dependence,
          latent.factor.skewness = latent.factor.skewness,
          latent.loglike = latent.loglike,
          gh_order_N = gh_order_N,
          gh.threads = gh.threads,
          approx.hermite.gh_order = approx.hermite.gh_order,
          approx.hermite.joint.gh_order = approx.hermite.joint.gh_order,
          approx.hermite.joint.order = approx.hermite.joint.order,
          approx.hermite.threads = approx.hermite.threads,
          qmc.Q = qmc.Q, qmc.n_rep = qmc.n_rep, qmc.seed = qmc.seed,
          qmc.chunk_size = qmc.chunk_size,
          zero.dependence = zero.dependence,
          zero.factor.skewness = zero.factor.skewness,
          gh_order_B = gh_order_B, gl_order_N = gl_order_N, gl_order_B = gl_order_B,
          biexp.root_tol = biexp.root_tol,
          psi0_start = psi0_start, rho0_start = rho0_start,
          gh.progressive = FALSE,
          progressive.sample.frac = NULL, progressive.sample.n = NULL,
          progressive.sample.seed = progressive.sample.seed,
          progressive.target.iter = 0L,
          progressive.target.ladder = FALSE,
          progressive.target.ladder.frac = progressive.target.ladder.frac,
          progressive.target.ladder.iter = 0L,
          .gh_progressive_internal = TRUE,
          .skip_final_qmc = TRUE,
          use.gradient = use.gradient, progress = FALSE, trace = 0,
          eval.max = max(50L, 5L * progressive.target.iter),
          iter.max = progressive.target.iter,
          tiny = tiny, big_penalty = big_penalty, rel.tol = rel.tol
        )

        # target_sample_fit is ALREADY at the final dimension J.
        # We only need to re-encode its fitted Psi as raw starting
        # parameters for the subsequent full-data fit; do not call
        # upgrade_factor_start(), which is only for q_new > q_old.
        psi0_start <- factor_start_from_Psi(
          target_sample_fit$count_layer$psi
        )
        if (isTRUE(target_sample_fit$zero_layer$dependence_estimated) &&
            !is.null(target_sample_fit$zero_layer$rho0)) {
          rho0_start <- target_sample_fit$zero_layer$rho0
        }
        progressive_history[[length(progressive_history) + 1L]] <- list(
          J = J, method = latent.loglike, n = stage_prepared$n,
          target_sample_warm = TRUE,
          logLik = target_sample_fit$copula_fit$log_lik_value,
          convergence = target_sample_fit$copula_fit$convergence,
          elapsed_sec = target_sample_fit$copula_fit$elapsed_sec
        )
        .internal_msg(
          "[dependence] warm target J=", J,
          " | n=", stage_prepared$n,
          " | iterations<=", progressive.target.iter,
          " | logLik=", signif(target_sample_fit$copula_fit$log_lik_value, 10),
          " | elapsed=", signif(target_sample_fit$copula_fit$elapsed_sec, 5), " sec"
        )
      }

      ## Optional nested sample ladder at the FINAL J.  This makes the
      ## transition to the full sample gradual (e.g. 10% -> 25% -> 50%
      ## -> 100%) and usually reduces the number of expensive full-data
      ## nlminb iterations.  I matrices are only subset, never recomputed.
      if (isTRUE(progressive.target.ladder) &&
          progressive.target.ladder.iter > 0L &&
          length(progressive_sample_index) < n) {

        current_idx <- progressive_sample_index
        current_frac <- length(current_idx) / n
        ladder_frac <- progressive.target.ladder.frac[
          progressive.target.ladder.frac > current_frac + 1e-12 &
          progressive.target.ladder.frac < 1 - 1e-12
        ]

        if (length(ladder_frac) > 0L) {
          for (ff in ladder_frac) {
            n_ff <- min(n - 1L, max(length(current_idx), as.integer(round(n * ff))))
            if (n_ff <= length(current_idx)) next

            remain <- setdiff(seq_len(n), current_idx)
            set.seed(as.integer(progressive.sample.seed) + n_ff + J)
            add_idx <- sample(remain, size = n_ff - length(current_idx), replace = FALSE)
            current_idx <- sort(c(current_idx, add_idx))
            ladder_prepared <- .subset_multi_zipm_I(prepared, current_idx)

            ladder_fit <- .fit_multizimp_copula(
              prepared = ladder_prepared,
              copula.block_size = min(copula.block_size, ladder_prepared$n),
              d0 = d0, J = J, L = L,
              latent.dependence = latent.dependence,
              latent.factor.skewness = latent.factor.skewness,
              latent.loglike = latent.loglike,
              gh_order_N = gh_order_N,
              gh.threads = gh.threads,
              approx.hermite.gh_order = approx.hermite.gh_order,
              approx.hermite.threads = approx.hermite.threads,
              qmc.Q = qmc.Q, qmc.n_rep = qmc.n_rep, qmc.seed = qmc.seed,
              qmc.chunk_size = qmc.chunk_size,
              zero.dependence = zero.dependence,
              zero.factor.skewness = zero.factor.skewness,
              gh_order_B = gh_order_B, gl_order_N = gl_order_N, gl_order_B = gl_order_B,
              biexp.root_tol = biexp.root_tol,
              psi0_start = psi0_start, rho0_start = rho0_start,
              gh.progressive = FALSE,
              progressive.sample.frac = NULL, progressive.sample.n = NULL,
              progressive.sample.seed = progressive.sample.seed,
              progressive.target.iter = 0L,
              progressive.target.ladder = FALSE,
              progressive.target.ladder.frac = progressive.target.ladder.frac,
              progressive.target.ladder.iter = 0L,
              .gh_progressive_internal = TRUE,
              .skip_final_qmc = TRUE,
              use.gradient = use.gradient, progress = FALSE, trace = 0,
              eval.max = max(30L, 5L * progressive.target.ladder.iter),
              iter.max = progressive.target.ladder.iter,
              tiny = tiny, big_penalty = big_penalty, rel.tol = rel.tol
            )

            psi0_start <- factor_start_from_Psi(ladder_fit$count_layer$psi)
            if (isTRUE(ladder_fit$zero_layer$dependence_estimated) &&
                !is.null(ladder_fit$zero_layer$rho0)) {
              rho0_start <- ladder_fit$zero_layer$rho0
            }

            progressive_history[[length(progressive_history) + 1L]] <- list(
              J = J, method = latent.loglike, n = ladder_prepared$n,
              target_ladder_warm = TRUE,
              logLik = ladder_fit$copula_fit$log_lik_value,
              convergence = ladder_fit$copula_fit$convergence,
              elapsed_sec = ladder_fit$copula_fit$elapsed_sec
            )

            .internal_msg(
              "[dependence] warm target J=", J,
              " | n=", ladder_prepared$n, "/", n,
              " | iterations<=", progressive.target.ladder.iter,
              " | logLik=", signif(ladder_fit$copula_fit$log_lik_value, 10),
              " | elapsed=", signif(ladder_fit$copula_fit$elapsed_sec, 5), " sec"
            )
          }
        }
      }

      .internal_msg("[dependence] final J=", J, " | n=", n, "/", n)
    }
  }

  ## =========================================================
  ## Progressive GH quadrature order for Psi/Rho estimation
  ##
  ## For exact Gaussian GH, estimate the SAME J-dimensional model
  ## successively with per-dimension orders
  ##     2 -> 4 -> 8 -> ... -> gh_order_N,
  ## appending the exact requested target when it is not a power of two.
  ## Each stage warm-starts BOTH Psi and (when L=1) Rho.  For L=1 the
  ## structural-zero GH order is automatically identical to the current
  ## count GH order, so there is no separate B-node tuning parameter.
  ## =========================================================
  ## Apply node-order warm-up to BOTH:
  ##   L = 0: N nodes 2 -> 4 -> 8 -> ... -> target
  ##   L = 1: N and B nodes jointly (2,2) -> (4,4) -> ... -> (target,target)
  if (identical(latent.dependence, "gaussian") &&
      identical(latent.loglike, "GH") &&
      as.integer(L) %in% c(0L, 1L) &&
      isTRUE(gh.progressive) &&
      !isTRUE(.gh_node_progressive_internal) &&
      gh_order_N > 2L) {

    gh_node_path <- 2L
    while (tail(gh_node_path, 1L) * 2L < gh_order_N) {
      gh_node_path <- c(gh_node_path, tail(gh_node_path, 1L) * 2L)
    }
    if (tail(gh_node_path, 1L) != gh_order_N) {
      gh_node_path <- c(gh_node_path, gh_order_N)
    }
    gh_node_warm <- head(gh_node_path, -1L)

    if (length(gh_node_warm) > 0L) {
      .internal_msg(
        "[dependence] progressive GH nodes: ",
        paste(gh_node_path, collapse = " -> "),
        if (estimate_B_factor && L == 1L && identical(zero.dependence, "gaussian"))
          " | B nodes = N nodes" else ""
      )

      node_psi0 <- psi0_start
      node_rho0 <- rho0_start

      for (gh_stage in gh_node_warm) {
        node_fit <- .fit_multizimp_copula(
          prepared = prepared,
          copula.block_size = copula.block_size,
          d0 = d0, J = J, L = L,
          latent.dependence = latent.dependence,
          latent.factor.skewness = latent.factor.skewness,
          latent.loglike = latent.loglike,
          gh_order_N = gh_stage,
          gh.threads = gh.threads,
          approx.hermite.gh_order = approx.hermite.gh_order,
          approx.hermite.joint.gh_order = approx.hermite.joint.gh_order,
          approx.hermite.joint.order = approx.hermite.joint.order,
          approx.hermite.threads = approx.hermite.threads,
          qmc.Q = qmc.Q, qmc.n_rep = qmc.n_rep, qmc.seed = qmc.seed,
          qmc.chunk_size = qmc.chunk_size,
          zero.dependence = zero.dependence,
          zero.factor.skewness = zero.factor.skewness,
          gh_order_B = gh_stage,
          gl_order_N = gl_order_N, gl_order_B = gl_order_B,
          biexp.root_tol = biexp.root_tol,
          psi0_start = node_psi0, rho0_start = node_rho0,
          gh.progressive = FALSE,
          progressive.sample.frac = NULL, progressive.sample.n = NULL,
          progressive.sample.seed = progressive.sample.seed,
          progressive.target.iter = 0L,
          progressive.target.ladder = FALSE,
          progressive.target.ladder.frac = numeric(0),
          progressive.target.ladder.iter = 0L,
          .gh_progressive_internal = TRUE,
          .gh_node_progressive_internal = TRUE,
          .skip_final_qmc = TRUE,
          use.gradient = use.gradient, progress = FALSE, trace = 0,
          eval.max = eval.max, iter.max = iter.max,
          tiny = tiny, big_penalty = big_penalty, rel.tol = rel.tol
        )

        node_psi0 <- factor_start_from_Psi(node_fit$count_layer$psi)
        if (isTRUE(node_fit$zero_layer$dependence_estimated) &&
            !is.null(node_fit$zero_layer$rho0)) {
          node_rho0 <- node_fit$zero_layer$rho0
        }

        progressive_history[[length(progressive_history) + 1L]] <- list(
          J = J, method = "GH", gh_order = gh_stage, n = n,
          node_warm = TRUE,
          logLik = node_fit$copula_fit$log_lik_value,
          convergence = node_fit$copula_fit$convergence,
          elapsed_sec = node_fit$copula_fit$elapsed_sec
        )

        .internal_msg(
          "[dependence] warm GH nodes=", gh_stage,
          if (estimate_B_factor && L == 1L && identical(zero.dependence, "gaussian"))
            paste0(" | B nodes=", gh_stage) else "",
          " | logLik=", signif(node_fit$copula_fit$log_lik_value, 10),
          " | iter=", if (!is.null(node_fit$copula_fit$optim$iterations)) node_fit$copula_fit$optim$iterations else NA_integer_,
          " | elapsed=", signif(node_fit$copula_fit$elapsed_sec, 5), " sec"
        )
      }

      psi0_start <- node_psi0
      rho0_start <- node_rho0
      .internal_msg(
        "[dependence] final GH nodes=", gh_order_N,
        if (estimate_B_factor && L == 1L && identical(zero.dependence, "gaussian"))
          paste0(" | B nodes=", gh_order_N) else ""
      )
    }
  }

  ## =========================================================
  ## Progressive Hermite quadrature order
  ##
  ## At the requested FINAL J, warm-start the SAME model through
  ##
  ##     Q = 2 -> 4 -> 8 -> ... -> Q_full,
  ##
  ## appending the exact requested Q_full when it is not a power of two.
  ## Every intermediate Q is limited to at most 10 optimizer iterations.
  ## The final Q_full is NOT fitted here: it falls through to the ordinary
  ## optimization below and therefore uses the user's full iter.max/eval.max.
  ##
  ## This is deliberately separate from progressive J:
  ## predecessor J stages above use Q=2 for Hermite, while this ladder is
  ## run only after the requested final J has been reached.
  ## =========================================================
  if (identical(latent.dependence, "gaussian") &&
      identical(latent.loglike, "hermite") &&
      isTRUE(gh.progressive) &&
      !isTRUE(.gh_node_progressive_internal) &&
      approx.hermite.joint.gh_order > 2L) {

    hermite_Q_full <- as.integer(approx.hermite.joint.gh_order)

    hermite_Q_path <- 2L
    while (tail(hermite_Q_path, 1L) * 2L < hermite_Q_full) {
      hermite_Q_path <- c(
        hermite_Q_path,
        tail(hermite_Q_path, 1L) * 2L
      )
    }
    if (tail(hermite_Q_path, 1L) != hermite_Q_full) {
      hermite_Q_path <- c(hermite_Q_path, hermite_Q_full)
    }

    ## Only intermediate Q values are warm stages.
    ## Q_full is optimized later with the full user controls.
    hermite_Q_warm <- head(hermite_Q_path, -1L)

    if (length(hermite_Q_warm) > 0L) {

      .internal_msg(
        "[dependence] progressive Hermite Q: ",
        paste(hermite_Q_path, collapse = " -> "),
        " | warm iterations<=10; final Q=", hermite_Q_full,
        " uses full optimization"
      )

      Q_psi0 <- psi0_start
      Q_rho0 <- rho0_start
      hermite_Q_warm_iter <- min(10L, as.integer(iter.max))

      for (Q_stage in hermite_Q_warm) {

        Q_fit <- .fit_multizimp_copula(
          prepared = prepared,
          copula.block_size = copula.block_size,
          d0 = d0,
          J = J,
          L = L,
          latent.dependence = latent.dependence,
          latent.factor.skewness = latent.factor.skewness,
          latent.loglike = latent.loglike,

          ## Keep all public quadrature-order aliases synchronized.
          gh_order_N = Q_stage,
          gh.threads = gh.threads,
          approx.hermite.gh_order = Q_stage,
          approx.hermite.joint.gh_order = Q_stage,
          approx.hermite.joint.order = approx.hermite.joint.order,
          approx.hermite.threads = approx.hermite.threads,

          qmc.Q = qmc.Q,
          qmc.n_rep = qmc.n_rep,
          qmc.seed = qmc.seed,
          qmc.chunk_size = qmc.chunk_size,

          zero.dependence = zero.dependence,
          zero.factor.skewness = zero.factor.skewness,
          gh_order_B = Q_stage,
          gl_order_N = Q_stage,
          gl_order_B = Q_stage,
          biexp.root_tol = biexp.root_tol,

          psi0_start = Q_psi0,
          rho0_start = Q_rho0,

          ## Prevent recursive progressive-J / progressive-Q ladders.
          gh.progressive = FALSE,
          progressive.sample.frac = NULL,
          progressive.sample.n = NULL,
          progressive.sample.seed = progressive.sample.seed,
          progressive.target.iter = 0L,
          progressive.target.ladder = FALSE,
          progressive.target.ladder.frac = numeric(0),
          progressive.target.ladder.iter = 0L,
          .gh_progressive_internal = TRUE,
          .gh_node_progressive_internal = TRUE,
          .skip_final_qmc = TRUE,

          use.gradient = use.gradient,
          progress = FALSE,
          trace = 0,

          ## Hard cap: no intermediate Q stage can exceed 10 iterations.
          eval.max = max(30L, 5L * hermite_Q_warm_iter),
          iter.max = hermite_Q_warm_iter,

          tiny = tiny,
          big_penalty = big_penalty,
          rel.tol = rel.tol
        )

        Q_psi0 <- factor_start_from_Psi(Q_fit$count_layer$psi)

        if (isTRUE(Q_fit$zero_layer$dependence_estimated) &&
            !is.null(Q_fit$zero_layer$rho0)) {
          Q_rho0 <- Q_fit$zero_layer$rho0
        }

        progressive_history[[length(progressive_history) + 1L]] <- list(
          J = J,
          method = "hermite",
          hermite_Q = Q_stage,
          n = n,
          Q_warm = TRUE,
          max_iterations = hermite_Q_warm_iter,
          logLik = Q_fit$copula_fit$log_lik_value,
          convergence = Q_fit$copula_fit$convergence,
          elapsed_sec = Q_fit$copula_fit$elapsed_sec
        )

        .internal_msg(
          "[dependence] warm Hermite Q=", Q_stage,
          " | iterations<=", hermite_Q_warm_iter,
          " | logLik=", signif(Q_fit$copula_fit$log_lik_value, 10),
          " | iter=",
          if (!is.null(Q_fit$copula_fit$optim$iterations))
            Q_fit$copula_fit$optim$iterations else NA_integer_,
          " | elapsed=", signif(Q_fit$copula_fit$elapsed_sec, 5), " sec"
        )
      }

      ## Warm-start the ordinary final-Q fit below.
      psi0_start <- Q_psi0
      rho0_start <- Q_rho0

      .internal_msg(
        "[dependence] final Hermite Q=", hermite_Q_full,
        " | full optimization"
      )
    }
  }

  ## =========================================================
  ## 3. Dependence-specific quadrature
  ##
  ## Counting and zero layers are allowed to use different
  ## factor copula families.
  ## =========================================================

  make_gh_grid <- function(q, gh_order) {

    if (!requireNamespace("statmod", quietly = TRUE)) {
      stop("Package 'statmod' is required.")
    }

    gh <- statmod::gauss.quad(
      gh_order,
      kind = "hermite"
    )

    idx <- expand.grid(
      rep(list(seq_len(gh_order)), q)
    )

    Zgh <- sqrt(2) * as.matrix(
      sapply(
        seq_len(q),
        function(r) gh$nodes[idx[[r]]]
      )
    )

    Wgh <- apply(
      as.matrix(
        sapply(
          seq_len(q),
          function(r) gh$weights[idx[[r]]]
        )
      ),
      1,
      prod
    ) / pi^(q / 2)

    list(z = Zgh, w = Wgh)
  }


  ## ---------------------------------------------------------
  ## Counting / latent-intensity factor layer
  ## ---------------------------------------------------------

  if (latent.dependence == "gaussian") {

    t_inner <- qnorm(
      seq_len(d0 - 1L) / d0
    )
    biexp_N <- NULL

    if (latent.loglike == "GH") {
      gh_N <- make_gh_grid(J, gh_order_N)
      Znode_N <- gh_N$z
      Wnode_N <- gh_N$w
      n_node_N <- length(Wnode_N)

    } else {
      ## Hermite polynomial does not construct a J-dimensional
      ## tensor Gauss-Hermite grid during optimization.
      gh_N <- NULL
      Znode_N <- matrix(0, nrow = 1L, ncol = J)
      Wnode_N <- 1
      n_node_N <- 1L

    }

  } else {

    gh_N <- NULL
    t_inner <- NULL

    biexp_N <- .make_biexp_factor_backend(
      lambda = latent.factor.lambda,
      q = J,
      d0 = d0,
      gl_order = gl_order_N,
      root_tol = biexp.root_tol
    )

    Znode_N <- biexp_N$Znode
    Wnode_N <- biexp_N$Wnode
    n_node_N <- biexp_N$n_node

  }


  ## ---------------------------------------------------------
  ## Structural-zero factor layer
  ##
  ## If L = 0, or fewer than two ZIPM margins exist, this
  ## layer is absent and zero.dependence is ignored.
  ## ---------------------------------------------------------

  if (!estimate_B_factor) {

    gh_B <- NULL
    biexp_B <- NULL

    Znode_B <- NULL
    Wnode_B <- 1
    n_node_B <- 1L

  } else if (zero.dependence == "gaussian") {

    ## The joint-Hermite backend integrates the Gaussian zero factors
    ## analytically inside the (J+L)-dimensional quadratic core.  Do not
    ## build the exponentially growing gh_order_B^L tensor grid.
    # For the nested Hermite architecture with L = 1 we MUST keep an
    # explicit one-dimensional GH grid for the structural-zero factor.
    # Only the legacy joint-(J+L) Hermite backend (reserved for L >= 2)
    # integrates the zero factors inside its Gaussian core and therefore
    # does not need Znode_B/Wnode_B.
    if (latent.dependence == "gaussian" && latent.loglike == "hermite" && L >= 2L) {
      gh_B <- NULL
      biexp_B <- NULL
      Znode_B <- NULL
      Wnode_B <- 1
      n_node_B <- 1L
    } else {
      gh_B <- make_gh_grid(L, gh_order_B)
      biexp_B <- NULL
      Znode_B <- gh_B$z
      Wnode_B <- gh_B$w
      n_node_B <- length(Wnode_B)
    }


  } else {

    gh_B <- NULL

    biexp_B <- .make_biexp_factor_backend(
      lambda = zero.factor.lambda,
      q = L,
      d0 = d0,
      gl_order = gl_order_B,
      root_tol = biexp.root_tol
    )

    Znode_B <- biexp_B$Znode
    Wnode_B <- biexp_B$Wnode
    n_node_B <- biexp_B$n_node

  }

  ## =========================================================
  ## 4. Fixed marginal quantities (precomputed once)
  ## =========================================================

  I_y_list <- prepared$I_y_list
  I_0_list <- prepared$I_0_list
  pb_list <- prepared$pb_list
  y0_list <- prepared$y0_list
  zi_position <- prepared$zi_position


  ## =========================================================
  ## Fixed row scaling for Hermite GEMMs
  ##
  ## I_y and I_0 are already precomputed in `prepared` and do not
  ## depend on Psi.  Hermite evaluates matrix products most stably
  ## after a per-observation row rescaling.  The split-I refactor
  ## must therefore rebuild ONLY these cheap scaled views here; it
  ## must NOT recompute any marginal interval probabilities.
  ## =========================================================

  I_y_scaled_list <- vector("list", k)
  log_I_y_scale_list <- vector("list", k)
  I_0_scaled_list <- vector("list", k)
  log_I_0_scale_list <- vector("list", k)

  if (latent.dependence == "gaussian" && latent.loglike == "hermite") {
    for (j in seq_len(k)) {
      Iy <- I_y_list[[j]]
      sy <- apply(Iy, 1L, max)
      ay <- is.finite(sy) & sy > 0

      Iys <- matrix(0, nrow = n, ncol = d0)
      if (any(ay)) {
        Iys[ay, ] <- sweep(Iy[ay, , drop = FALSE], 1L, sy[ay], "/")
      }

      lsy <- rep(log(tiny), n)
      lsy[ay] <- log(sy[ay])

      I_y_scaled_list[[j]] <- Iys
      log_I_y_scale_list[[j]] <- lsy

      if (zero_inflation[j]) {
        I0 <- I_0_list[[j]]
        s0 <- apply(I0, 1L, max)
        a0 <- is.finite(s0) & s0 > 0

        I0s <- matrix(0, nrow = n, ncol = d0)
        if (any(a0)) {
          I0s[a0, ] <- sweep(I0[a0, , drop = FALSE], 1L, s0[a0], "/")
        }

        ls0 <- rep(log(tiny), n)
        ls0[a0] <- log(s0[a0])

        I_0_scaled_list[[j]] <- I0s
        log_I_0_scale_list[[j]] <- ls0
      }
    }
  }

  ## Fixed Gaussian zero-layer quantiles.  These depend only on the
  ## already-fitted marginal zero probabilities and can be cached once.
  qpb_list <- vector("list", k)
  if (estimate_B_factor && zero.dependence == "gaussian") {
    for (j in zi_index) {
      qpb_list[[j]] <- qnorm(
        pmin(1 - 1e-12, pmax(1e-12, pb_list[[j]]))
      )
    }
  }

  ## =========================================================
  ## 5. Chain rule helper
  ## =========================================================

  chain_grad_loading_to_raw <- function(
      grad_Load,
      Araw,
      map,
      npar,
      q
  ) {

    grad_raw <- numeric(npar)

    for (id in seq_len(npar)) {

      j <- map[id, 1L]
      r <- map[id, 2L]

      a_j <- Araw[j, ]
      c_j <- sqrt(1 + sum(a_j^2))

      Jac <- diag(1 / c_j, q) -
        tcrossprod(a_j, a_j) / c_j^3

      dA_draw <- if (j == r) Araw[j, r] else 1

      grad_raw[id] <- sum(
        grad_Load[j, ] * Jac[, r] * dA_draw
      )
    }

    grad_raw
  }

  ## =========================================================
  ## 6. Count-factor kernel: node x checkerboard-cell only
  ##
  ## MEMORY DESIGN
  ##
  ## The old implementation stored
  ##
  ##   dPN_j : n x n_node_N x J
  ##
  ## for every margin j.
  ##
  ## This is unnecessary.  The derivative of the conditional
  ## checkerboard cell probabilities does NOT depend on
  ## observation i.  We therefore store only
  ##
  ##   Omega_j        : n_node_N x d0
  ##   dOmega_j/dr    : n_node_N x d0
  ##
  ## and form
  ##
  ##   PN_j      = d0 * I_j %*% t(Omega_j)
  ##   dPN_j/dr  = d0 * I_j %*% t(dOmega_j/dr)
  ##
  ## only for the current observation block / margin / factor.
  ## =========================================================

  compute_N_kernel <- function(
      psi0,
      need_grad = FALSE
  ) {

    tmp <- make_factor_loadings(
      par = psi0,
      k = k,
      q = J,
      return_all = TRUE
    )

    Psi <- tmp$Psi
    Araw <- tmp$A
    map <- tmp$map

    norm2 <- rowSums(Psi^2)
    svec <- sqrt(pmax(0, 1 - norm2))

    if (any(!is.finite(svec)) ||
        any(svec <= 0)) {
      stop("Invalid count-factor loading matrix.")
    }

    Omega_list <- vector("list", k)

    dOmega_list <- if (need_grad) {
      vector("list", k)
    } else {
      NULL
    }

    if (latent.dependence == "gaussian") {

      eta <- Znode_N %*% t(Psi)

      for (j in seq_len(k)) {

        A_inner <- sweep(
          matrix(
            t_inner,
            nrow = n_node_N,
            ncol = d0 - 1L,
            byrow = TRUE
          ),
          1,
          eta[, j],
          "-"
        )

        aval <- A_inner / svec[j]
        q_inner <- pnorm(aval)

        q_full <- cbind(
          0,
          q_inner,
          1
        )

        Omega_j <-
          q_full[, 2:(d0 + 1L), drop = FALSE] -
          q_full[, 1:d0, drop = FALSE]

        Omega_list[[j]] <- Omega_j

        if (need_grad) {

          phi_aval <- dnorm(aval)

          dOmega_j <- array(
            0,
            dim = c(
              n_node_N,
              d0,
              J
            )
          )

          for (r in seq_len(J)) {

            da <- sweep(
              A_inner * Psi[j, r] / svec[j]^3,
              1,
              Znode_N[, r] / svec[j],
              "-"
            )

            dq <- phi_aval * da

            dq_full <- cbind(
              0,
              dq,
              0
            )

            dOmega_j[, , r] <-
              dq_full[, 2:(d0 + 1L), drop = FALSE] -
              dq_full[, 1:d0, drop = FALSE]
          }

          dOmega_list[[j]] <- dOmega_j
        }
      }

    } else {

      for (j in seq_len(k)) {

        psi_j <- Psi[j, ]

        Qobj <- biexp_N$compute_Q_and_dQ(
          psi_j = psi_j
        )

        Qj <- Qobj$Q
        dQj <- Qobj$dQ
        sigma_j <- Qobj$sigma
        eta_j <- Qobj$eta

        Q_inner <- Qj[2:d0]

        A_inner <-
          matrix(
            Q_inner,
            nrow = n_node_N,
            ncol = d0 - 1L,
            byrow = TRUE
          ) -
          matrix(
            eta_j,
            nrow = n_node_N,
            ncol = d0 - 1L
          )

        aval <- A_inner / sigma_j
        q_inner <- biexp_N$F(aval)

        q_full <- cbind(
          0,
          q_inner,
          1
        )

        Omega_j <-
          q_full[, 2:(d0 + 1L), drop = FALSE] -
          q_full[, 1:d0, drop = FALSE]

        Omega_list[[j]] <- Omega_j

        if (need_grad) {

          fe_inner <- biexp_N$f(aval)

          dOmega_j <- array(
            0,
            dim = c(
              n_node_N,
              d0,
              J
            )
          )

          for (r in seq_len(J)) {

            dQ_r <- dQj[2:d0, r]

            dQ_mat <- matrix(
              dQ_r,
              nrow = n_node_N,
              ncol = d0 - 1L,
              byrow = TRUE
            )

            Zr_mat <- matrix(
              Znode_N[, r],
              nrow = n_node_N,
              ncol = d0 - 1L
            )

            da <-
              (dQ_mat - Zr_mat) / sigma_j +
              A_inner * psi_j[r] / sigma_j^3

            dq <- fe_inner * da

            dq_full <- cbind(
              0,
              dq,
              0
            )

            dOmega_j[, , r] <-
              dq_full[, 2:(d0 + 1L), drop = FALSE] -
              dq_full[, 1:d0, drop = FALSE]
          }

          dOmega_list[[j]] <- dOmega_j
        }
      }
    }

    # These transposes are reused for every observation block.
    # Building them once avoids repeated t() calls in make_PN_block().
    OmegaT_list <- lapply(Omega_list, t)

    # For the gradient, stack all J derivative directions into one
    # d0 x (n_node_N * J) matrix per margin.  This lets each block
    # compute all dPN/dr with a single BLAS matrix multiplication.
    dOmegaT_stack_list <- if (need_grad) {
      lapply(
        dOmega_list,
        function(dOm) {
          do.call(
            cbind,
            lapply(seq_len(J), function(r) t(dOm[, , r]))
          )
        }
      )
    } else {
      NULL
    }

    list(
      Omega_list = Omega_list,
      OmegaT_list = OmegaT_list,
      dOmegaT_stack_list = dOmegaT_stack_list,
      Psi = Psi,
      A = Araw,
      map = map
    )
  }


  ## =========================================================
  ## 7. Zero-factor kernel
  ##
  ## No n x n_node_B x L derivative arrays are stored.
  ##
  ## Gaussian:
  ##   PB and dPB are generated only for the requested
  ##   observation block / node / loading.
  ##
  ## BiExp:
  ##   the expensive marginal quantiles Q_S(p_B) and dQ/dRho
  ##   are retained, but these require only O(n L) memory per
  ##   ZIPM margin rather than O(n n_node_B L).
  ## =========================================================

  compute_B_kernel <- function(
      rho0,
      need_grad = FALSE
  ) {

    if (!estimate_B_factor) {
      stop(
        "compute_B_kernel() is only available when zero dependence is estimated."
      )
    }

    tmp <- make_factor_loadings(
      par = rho0,
      k = n_zi,
      q = L,
      return_all = TRUE
    )

    Rho <- tmp$Psi
    Araw <- tmp$A
    map <- tmp$map

    svec <- sqrt(
      pmax(
        0,
        1 - rowSums(Rho^2)
      )
    )

    if (any(!is.finite(svec)) ||
        any(svec <= 0)) {
      stop("Invalid zero-factor loading matrix.")
    }

    state <- NULL

    if (zero.dependence == "biexp") {

      state <- vector(
        "list",
        n_zi
      )

      for (s in seq_len(n_zi)) {

        j <- zi_index[s]

        pb_j <- pmin(
          1 - 1e-12,
          pmax(
            1e-12,
            pb_list[[j]]
          )
        )

        rho_s <- Rho[s, ]

        Qobj <- biexp_B$compute_prob_Q_and_dQ(
          loading = rho_s,
          prob = pb_j
        )

        state[[s]] <- list(
          Q = Qobj$Q,
          dQ = if (need_grad) {
            Qobj$dQ
          } else {
            NULL
          },
          sigma = Qobj$sigma,
          eta = as.vector(
            Znode_B %*% rho_s
          )
        )
      }
    }

    list(
      Rho = Rho,
      svec = svec,
      A = Araw,
      map = map,
      state = state
    )
  }


  ## ---------------------------------------------------------
  ## Conditional structural-zero active probability for one
  ## observation block and one zero-factor quadrature node.
  ##
  ## If need_grad = TRUE, only ONE derivative direction l is
  ## returned; no 3-D derivative array is created.
  ## ---------------------------------------------------------

  get_PB_block <- function(
      j,
      b,
      idx,
      obj_B = NULL,
      need_grad = FALSE,
      l = NULL
  ) {

    if (!zero_inflation[j]) {

      return(
        list(
          PB = rep(1, length(idx)),
          dPB = NULL
        )
      )
    }

    if (!estimate_B_factor) {

      return(
        list(
          PB = pb_list[[j]][idx],
          dPB = NULL
        )
      )
    }

    s <- zi_position[j]

    if (s <= 0L) {
      stop("Internal ZIPM-margin indexing error.")
    }

    if (need_grad) {

      if (is.null(l) ||
          length(l) != 1L ||
          l < 1L ||
          l > L) {
        stop("A valid zero-factor derivative index `l` is required.")
      }
    }

    rho_s <- obj_B$Rho[s, ]
    sigma_s <- obj_B$svec[s]

    if (zero.dependence == "gaussian") {

      qpb <- qpb_list[[j]][idx]

      eta_b <- sum(
        Znode_B[b, ] *
          rho_s
      )

      A_pb <- qpb - eta_b
      aval <- A_pb / sigma_s
      PB <- pnorm(aval)

      dPB <- NULL

      if (need_grad) {

        da <-
          A_pb *
          rho_s[l] /
          sigma_s^3 -
          Znode_B[b, l] /
          sigma_s

        dPB <- dnorm(aval) * da
      }

    } else {

      st <- obj_B$state[[s]]

      Qpb <- st$Q[idx]
      eta_b <- st$eta[b]

      A_pb <- Qpb - eta_b
      aval <- A_pb / sigma_s

      PB <- biexp_B$F(aval)

      dPB <- NULL

      if (need_grad) {

        dQ_l <- st$dQ[idx, l]

        da <-
          (dQ_l - Znode_B[b, l]) /
          sigma_s +
          A_pb *
          rho_s[l] /
          sigma_s^3

        dPB <- biexp_B$f(aval) * da
      }
    }

    list(
      PB = PB,
      dPB = dPB
    )
  }


  # Derivative-only version used by the zero-factor gradient.
  # It avoids recomputing PB itself when only dPB/d loading is needed.
  get_dPB_block <- function(
      j,
      b,
      idx,
      obj_B,
      l
  ) {

    s <- zi_position[j]
    if (s <= 0L) {
      stop("Internal ZIPM-margin indexing error.")
    }

    rho_s <- obj_B$Rho[s, ]
    sigma_s <- obj_B$svec[s]

    if (zero.dependence == "gaussian") {

      qpb <- qpb_list[[j]][idx]
      eta_b <- sum(Znode_B[b, ] * rho_s)
      A_pb <- qpb - eta_b
      aval <- A_pb / sigma_s

      da <-
        A_pb * rho_s[l] / sigma_s^3 -
        Znode_B[b, l] / sigma_s

      return(dnorm(aval) * da)
    }

    st <- obj_B$state[[s]]
    Qpb <- st$Q[idx]
    eta_b <- st$eta[b]
    A_pb <- Qpb - eta_b
    aval <- A_pb / sigma_s
    dQ_l <- st$dQ[idx, l]

    da <-
      (dQ_l - Znode_B[b, l]) / sigma_s +
      A_pb * rho_s[l] / sigma_s^3

    biexp_B$f(aval) * da
  }

  ## =========================================================
  ## 8. Small block helpers
  ## =========================================================

  ## Fast cache for the only estimated zero-factor case used by this
  ## implementation: L = 1.  For a fixed parameter vector and
  ## observation block, PB and dPB/dRho depend on (ZIPM margin, B-node)
  ## but NOT on the count-factor node.  The old code recomputed pnorm(),
  ## dnorm(), eta and the same vectors repeatedly in likelihood and
  ## gradient passes.  Cache them once per block.
  make_PB_cache_L1 <- function(idx, obj_B, need_grad = FALSE) {

    if (!estimate_B_factor) return(NULL)
    if (L != 1L) stop("Internal error: make_PB_cache_L1() requires L = 1.")

    nb <- length(idx)
    PB <- vector("list", n_zi)
    dPB <- if (isTRUE(need_grad)) vector("list", n_zi) else NULL

    for (s in seq_len(n_zi)) {
      j <- zi_index[s]
      rho <- obj_B$Rho[s, 1L]
      sig <- obj_B$svec[s]

      if (zero.dependence == "gaussian") {
        qpb <- qpb_list[[j]][idx]
        eta <- rho * Znode_B[, 1L]
        A <- outer(qpb, eta, "-")
        aval <- A / sig
        P <- pnorm(aval)
        PB[[s]] <- P

        if (isTRUE(need_grad)) {
          # da/drho = A*rho/sigma^3 - z/sigma
          DA <- A * (rho / sig^3) -
            matrix(Znode_B[, 1L] / sig, nrow = nb, ncol = n_node_B, byrow = TRUE)
          dPB[[s]] <- dnorm(aval) * DA
        }
      } else {
        st <- obj_B$state[[s]]
        Qpb <- st$Q[idx]
        eta <- st$eta
        A <- outer(Qpb, eta, "-")
        aval <- A / sig
        PB[[s]] <- matrix(biexp_B$F(as.vector(aval)), nrow = nb, ncol = n_node_B)

        if (isTRUE(need_grad)) {
          dQ <- st$dQ[idx, 1L]
          DA <- outer(dQ, Znode_B[, 1L], "-") / sig + A * (rho / sig^3)
          dPB[[s]] <- matrix(biexp_B$f(as.vector(aval)), nrow = nb, ncol = n_node_B) * DA
        }
      }
    }

    list(PB = PB, dPB = dPB)
  }

  get_PB_cached <- function(j, b, idx, obj_B, cache = NULL) {
    if (!zero_inflation[j]) return(rep(1, length(idx)))
    if (!estimate_B_factor) return(pb_list[[j]][idx])
    if (!is.null(cache)) return(cache$PB[[zi_position[j]]][, b])
    get_PB_block(j = j, b = b, idx = idx, obj_B = obj_B, need_grad = FALSE)$PB
  }

  get_dPB_cached <- function(j, b, idx, obj_B, cache = NULL) {
    if (!is.null(cache) && !is.null(cache$dPB)) {
      return(cache$dPB[[zi_position[j]]][, b])
    }
    get_dPB_block(j = j, b = b, idx = idx, obj_B = obj_B, l = 1L)
  }

  make_PN_block <- function(
      idx,
      obj_N
  ) {

    nb <- length(idx)

    PN_y <- vector(
      "list",
      k
    )

    PN_0 <- vector(
      "list",
      k
    )

    for (j in seq_len(k)) {

      OmegaT_j <- obj_N$OmegaT_list[[j]]

      PN_y[[j]] <-
        d0 *
        (
          I_y_list[[j]][idx, , drop = FALSE] %*%
            OmegaT_j
        )

      # PN_0 is used only by ZIPM margins.
      if (zero_inflation[j]) {

        PN_0[[j]] <-
          d0 *
          (
            I_0_list[[j]][idx, , drop = FALSE] %*%
              OmegaT_j
          )

      } else {

        PN_0[[j]] <- NULL
      }

      if (any(!is.finite(PN_y[[j]]))) {
        stop("Non-finite PN_y values.")
      }

      if (zero_inflation[j] &&
          any(!is.finite(PN_0[[j]]))) {
        stop("Non-finite PN_0 values.")
      }
    }

    list(
      PN_y = PN_y,
      PN_0 = PN_0,
      n = nb
    )
  }


  make_T_block <- function(
      j,
      b,
      idx,
      PN,
      obj_B,
      PB_b = NULL
  ) {

    if (!zero_inflation[j]) {

      return(
        PN$PN_y[[j]]
      )
    }

    if (is.null(PB_b)) {
      PB_b <- get_PB_block(
        j = j,
        b = b,
        idx = idx,
        obj_B = obj_B,
        need_grad = FALSE
      )$PB
    }

    y0 <- y0_list[[j]][idx]
    ypos <- !y0

    Tj <- matrix(
      0,
      nrow = length(idx),
      ncol = n_node_N
    )

    if (any(y0)) {

      Tj[y0, ] <-
        (1 - PB_b[y0]) +
        sweep(
          PN$PN_0[[j]][y0, , drop = FALSE],
          1,
          PB_b[y0],
          "*"
        )
    }

    if (any(ypos)) {

      Tj[ypos, ] <-
        sweep(
          PN$PN_y[[j]][ypos, , drop = FALSE],
          1,
          PB_b[ypos],
          "*"
        )
    }

    Tj
  }


  make_block_index <- function() {

    split(
      seq_len(n),
      ceiling(
        seq_len(n) /
          copula.block_size
      )
    )
  }

  obs_blocks <- make_block_index()


  ## =========================================================
  ## 8a. Hermite polynomial Gaussian-factor approximation
  ## =========================================================

  approx.hermite.gh_order <- as.integer(round(approx.hermite.gh_order))
  gh.threads <- as.integer(round(gh.threads))
  if (!is.finite(gh.threads) || gh.threads < 1L) {
    stop("`gh.threads` must be an integer >= 1.")
  }

  approx.hermite.threads <- as.integer(round(approx.hermite.threads))

  if (approx.hermite.gh_order < 2L) {
    stop("`approx.hermite.gh_order` must be an integer >= 2.")
  }
  if (!is.finite(approx.hermite.threads) || approx.hermite.threads < 1L) {
    stop("`approx.hermite.threads` must be an integer >= 1.")
  }

  ## Hermite polynomial uses one fixed standard-normal GH grid.
  ## No J-dimensional tensor grid is constructed.
  if (latent.dependence == "gaussian" && latent.loglike == "hermite") {
    herm_order_use <- if (isTRUE(estimate_B_factor) && L >= 1L &&
                          identical(zero.dependence, "gaussian")) {
      approx.hermite.joint.gh_order
    } else {
      approx.hermite.gh_order
    }
    gh_herm1 <- make_gh_grid(1L, herm_order_use)
    herm_x <- as.vector(gh_herm1$z)
    herm_w <- as.vector(gh_herm1$w)
    if (!isTRUE(.gh_progressive_internal) && L == 0L) {
      .internal_msg(
        "[dependence] Hermite L=0 kernels: C++ interval + grid + core | threads=",
        approx.hermite.threads, " | projection Q=", length(herm_x)
      )
    }
  } else {
    herm_x <- herm_w <- NULL
  }


  ## ---------------------------------------------------------
  ## Conditional checkerboard probability and first two
  ## derivatives with respect to eta = psi_j' z.
  ##
  ## For q_r(eta) = Phi((t_r-eta)/s):
  ##
  ##   q'_r  = -phi(a_r)/s
  ##   q''_r = -a_r phi(a_r)/s^2.
  ##
  ## The code returns H0,H1,H2 = h, dh/deta, d2h/deta2 as
  ## n_b x k matrices.  Observations are vectorized inside each
  ## margin; there is no observation loop.
  ## ---------------------------------------------------------

  ## Faster special case at eta=0.  For each margin the
  ## checkerboard derivative vectors are scalar d0-vectors, so the
  ## n_b observations are evaluated by BLAS matrix-vector products.
  ## ---------------------------------------------------------
  ## Exact Gaussian integral of an exponential quadratic:
  ##
  ## E exp{ const + a'Z + 1/2 Z' B Z }
  ## =
  ## exp(const) |I-B|^{-1/2}
  ## exp{ 1/2 a'(I-B)^{-1}a }.
  ##
  ## QROW contains row-wise coefficients q_ij such that
  ## B_i = bscale * Psi' diag(q_i) Psi.
  ##
  ## The expensive construction of all B_i is vectorized through
  ## row-wise Kronecker features.  Only the small J x J Cholesky
  ## factorization remains observation-specific.
  ## ---------------------------------------------------------

  ## ---------------------------------------------------------
  ## Degree-2 Hermite projection of log h_j.
  ##
  ## Let X ~ N(0,1), eta_j = ||psi_j|| X.  We project
  ##
  ## log h_j(||psi_j|| X)
  ##   ~= c0 + c1 H1(X) + c2 H2(X)
  ##
  ## with H1=X, H2=X^2-1.  This becomes a quadratic form in z,
  ## which is integrated analytically by expquad_logint_batch().
  ##
  ## Only a one-dimensional GH rule is used for the projection;
  ## there is no tensor product in J.
  ## ---------------------------------------------------------

  ## ---------------------------------------------------------
  ## Hermite polynomial projection weights for the quadratic core.
  ## Only H0, H1, H2 are required; the remaining residual is handled
  ## by the first-cumulant correction under the tilted Gaussian.
  ## ---------------------------------------------------------

  if (latent.dependence == "gaussian") {
    herm_proj_W <- cbind(
      herm_w,
      herm_w * herm_x,
      0.5 * herm_w * (herm_x^2 - 1)
    )
  } else {
    herm_proj_W <- NULL
  }


  ## ============================================================
  ## Fast Hermite kernel.
  ##
  ## For fixed Psi, Omega_H,j and dOmega_H,j / d r_j do NOT depend
  ## on the observation.  Construct them once per objective/gradient
  ## evaluation, then reuse them over all observation blocks.
  ##
  ## r_j = ||psi_j||,
  ## s_j = sqrt(1-r_j^2),
  ## eta_h = r_j x_h,
  ##
  ## a_hm = (t_m - r_j x_h) / s_j.
  ##
  ## da/dr =
  ##   -x_h/s_j + (t_m-r_j x_h) r_j / s_j^3.
  ## ============================================================

  compute_H_kernel <- function(Psi) {
    if (!exists("hermite_interval_kernel_cpp", mode = "function", inherits = TRUE)) {
      stop(
        paste0(
          "Hermite count-cell kernel is not loaded. Compile hermite_count_kernels.cpp first with ",
          "Rcpp::sourceCpp('01_hermite_interval.cpp')."
        )
      )
    }
    hermite_interval_kernel_cpp(
      Psi = as.matrix(Psi),
      herm_x = herm_x,
      t_inner = t_inner,
      d0 = as.integer(d0),
      n_threads = as.integer(approx.hermite.threads)
    )
  }


  ## ============================================================
  ## Hermite polynomial: quadratic core + first-cumulant correction on a FIXED
  ## standard-normal Hermite grid.
  ##
  ## Let f_ij(x) = log h_ij(r_j x).  The quadratic Hermite core
  ## (m=0,1,2) defines the tilted Gaussian marginal
  ##
  ##     X_ij ~ N(MUX_ij, VARX_ij).
  ##
  ## Instead of moving the GH nodes observation by observation, use
  ## the density ratio
  ##
  ##     W_ij(x) = phi(x; MUX_ij, VARX_ij) / phi(x; 0, 1)
  ##
  ## so that
  ##
  ##     E_q[f_ij(X)] = E_phi[f_ij(X) W_ij(X)].
  ##
  ## Hence the SAME fixed nodes `herm_x` and the SAME matrices LH and
  ## dLH used by the finite-M Hermite path are reused for every
  ## observation.  No observation-specific checkerboard kernel is
  ## evaluated.
  ##
  ## The only approximation beyond the fixed 1-D quadrature is still
  ##
  ##     log E_q exp(R) ~= E_q[R].
  ##
  ## All gradients are analytic.  Differentiating the density ratio:
  ##
  ##   d log W / d mu = (x-mu)/v,
  ##   d log W / d v  = ((x-mu)^2-v)/(2 v^2).
  ## ============================================================


  ## -----------------------------------------------------------------
  ## L=1 Hermite acceleration: count-layer quantities do not depend on
  ## the zero-factor quadrature node b.  The old implementation
  ## recomputed I %*% Omega (and its radial derivative) Q_B times per
  ## observation block.  Cache them ONCE per block / parameter value.
  ## -----------------------------------------------------------------
  make_hermite_count_cache_block <- function(idx, Hker, need_grad = FALSE) {
    if (!exists("hermite_count_cache_cpp", mode = "function", inherits = TRUE)) {
      stop(
        paste0(
          "Hermite count-cache kernel is not loaded. Compile hermite_count_kernels.cpp first with ",
          "Rcpp::sourceCpp('01_hermite_interval.cpp')."
        )
      )
    }
    hermite_count_cache_cpp(
      Iy_list = I_y_scaled_list,
      log_sy_list = log_I_y_scale_list,
      I0_list = I_0_scaled_list,
      log_s0_list = log_I_0_scale_list,
      OmegaT_list = Hker$OmegaT_list,
      dOmegaT_dr_list = if (isTRUE(need_grad)) Hker$dOmegaT_dr_list else vector("list", k),
      dOmegaT_dx_list = if (isTRUE(need_grad)) Hker$dOmegaT_dx_list else vector("list", k),
      zero_inflation = zero_inflation,
      idx1 = as.integer(idx),
      d0 = d0,
      need_grad = isTRUE(need_grad),
      tiny = tiny,
      n_threads = as.integer(approx.hermite.threads)
    )
  }


  ## Build f_ij(x_q), its exact radial derivative, the quadratic
  ## coefficients c0,c1,c2, and (when needed) the derivative with
  ## respect to the zero-layer probability PB, all on the fixed
  ## standard-normal Hermite grid.
  make_hermite_grid_block <- function(
      idx,
      b,
      obj_B,
      Hker,
      need_grad = FALSE,
      PB_cache = NULL,
      count_cache = NULL
  ) {

    nb <- length(idx)
    QH <- length(herm_x)

    ## Fast L=0 path: all margin-wise work is fused in C++/OpenMP.
    if (!isTRUE(estimate_B_factor) && identical(as.integer(L), 0L)) {
      if (!exists("hermite_grid_L0_cpp", mode = "function", inherits = TRUE)) {
        stop(
          paste0(
            "Hermite L=0 grid kernel is not loaded. Compile hermite_count_kernels.cpp first with ",
            "Rcpp::sourceCpp('01_hermite_interval.cpp')."
          )
        )
      }
      return(
        hermite_grid_L0_cpp(
          Iy_list = I_y_scaled_list,
          log_sy_list = log_I_y_scale_list,
          I0_list = I_0_scaled_list,
          log_s0_list = log_I_0_scale_list,
          OmegaT_list = Hker$OmegaT_list,
          dOmegaT_dr_list = if (isTRUE(need_grad)) Hker$dOmegaT_dr_list else vector("list", k),
          pb_list = pb_list,
          y0_list = y0_list,
          zero_inflation = zero_inflation,
          idx1 = as.integer(idx),
          Wproj = herm_proj_W[, 1:3, drop = FALSE],
          d0 = d0,
          need_grad = isTRUE(need_grad),
          tiny = tiny,
          n_threads = as.integer(approx.hermite.threads)
        )
      )
    }

    C0 <- matrix(0, nrow = nb, ncol = k)
    C1 <- matrix(0, nrow = nb, ncol = k)
    C2 <- matrix(0, nrow = nb, ncol = k)

    D0 <- matrix(0, nrow = nb, ncol = k)
    D1 <- matrix(0, nrow = nb, ncol = k)
    D2 <- matrix(0, nrow = nb, ncol = k)

    LH_list <- vector("list", k)
    dLH_list <- if (isTRUE(need_grad)) vector("list", k) else NULL
    FPB_list <- if (isTRUE(need_grad) && estimate_B_factor) vector("list", k) else NULL

    dC1_dPB <- matrix(0, nrow = nb, ncol = k)
    dC2_dPB <- matrix(0, nrow = nb, ncol = k)

    Wproj <- herm_proj_W[, 1:3, drop = FALSE]

    log_count_scaled <- function(Is, logsc, OmT) {

      Ps <- d0 * (Is %*% OmT)

      logP <- matrix(
        log(tiny),
        nrow = nrow(Is),
        ncol = ncol(Ps)
      )

      goodP <- is.finite(Ps) & Ps > 0
      if (any(goodP)) logP[goodP] <- log(Ps[goodP])

      sweep(logP, 1L, logsc, "+")
    }

    count_value_score_scaled <- function(Is, logsc, OmT, dOmT) {

      Ps  <- d0 * (Is %*% OmT)
      dPs <- d0 * (Is %*% dOmT)

      logP <- matrix(
        log(tiny),
        nrow = nrow(Is),
        ncol = ncol(Ps)
      )

      dlogP <- matrix(0, nrow = nrow(Is), ncol = ncol(Ps))

      goodP <- is.finite(Ps) & Ps > 0
      if (any(goodP)) {
        logP[goodP] <- log(Ps[goodP])
        dlogP[goodP] <- dPs[goodP] / Ps[goodP]
      }

      logP <- sweep(logP, 1L, logsc, "+")

      list(logP = logP, dlogP = dlogP)
    }


    for (j in seq_len(k)) {

      OmT <- Hker$OmegaT_list[[j]]

      if (!is.null(count_cache)) {
        logY <- count_cache$logY[[j]]
        dlogY <- if (isTRUE(need_grad)) count_cache$dlogY[[j]] else NULL
      } else {
        Iys <- I_y_scaled_list[[j]][idx, , drop = FALSE]
        lsy <- log_I_y_scale_list[[j]][idx]

        if (isTRUE(need_grad)) {
          YY <- count_value_score_scaled(
            Is = Iys,
            logsc = lsy,
            OmT = OmT,
            dOmT = Hker$dOmegaT_dr_list[[j]]
          )
          logY <- YY$logP
          dlogY <- YY$dlogP
        } else {
          logY <- log_count_scaled(Iys, lsy, OmT)
          dlogY <- NULL
        }
      }

      FPB <- if (isTRUE(need_grad) && estimate_B_factor && zero_inflation[j]) {
        matrix(0, nrow = nb, ncol = QH)
      } else {
        NULL
      }

      if (!zero_inflation[j]) {

        LH <- logY
        dLH <- if (isTRUE(need_grad)) dlogY else NULL

      } else {

        PB <- get_PB_cached(
          j = j, b = b, idx = idx, obj_B = obj_B, cache = PB_cache
        )

        PB <- pmin(pmax(PB, 1e-12), 1 - 1e-12)

        y0 <- y0_list[[j]][idx]
        yp <- !y0

        LH <- matrix(log(tiny), nrow = nb, ncol = QH)
        dLH <- if (isTRUE(need_grad)) matrix(0, nrow = nb, ncol = QH) else NULL

        if (any(yp)) {

          LH[yp, ] <-
            sweep(
              logY[yp, , drop = FALSE],
              1L,
              log(PB[yp]),
              "+"
            )

          if (isTRUE(need_grad)) {
            dLH[yp, ] <- dlogY[yp, , drop = FALSE]
          }

          if (!is.null(FPB)) {
            FPB[yp, ] <- matrix(
              1 / PB[yp],
              nrow = sum(yp),
              ncol = QH
            )
          }
        }

        if (any(y0)) {

          if (!is.null(count_cache)) {
            log0 <- count_cache$log0[[j]][y0, , drop = FALSE]
            ZZ <- if (isTRUE(need_grad)) {
              list(dlogP = count_cache$dlog0[[j]][y0, , drop = FALSE])
            } else NULL
          } else {
            I0s <- I_0_scaled_list[[j]][idx[y0], , drop = FALSE]
            ls0 <- log_I_0_scale_list[[j]][idx[y0]]

            if (isTRUE(need_grad)) {
              ZZ <- count_value_score_scaled(
                Is = I0s,
                logsc = ls0,
                OmT = OmT,
                dOmT = Hker$dOmegaT_dr_list[[j]]
              )
              log0 <- ZZ$logP
            } else {
              log0 <- log_count_scaled(I0s, ls0, OmT)
              ZZ <- NULL
            }
          }

          pb0 <- PB[y0]

          logA <- log1p(-pb0)
          logB <- sweep(log0, 1L, log(pb0), "+")

          logA_mat <- matrix(
            logA,
            nrow = length(logA),
            ncol = QH
          )

          mx <- pmax(logA_mat, logB)
          logTot <-
            mx +
            log(
              exp(logA_mat - mx) +
              exp(logB - mx)
            )

          LH[y0, ] <- logTot

          if (isTRUE(need_grad)) {
            wCount <- exp(logB - logTot)
            dLH[y0, ] <- wCount * ZZ$dlogP
          }

          if (!is.null(FPB)) {
            P0 <- exp(pmin(log0, 0))
            Tot <- pmax(exp(logTot), tiny)
            FPB[y0, ] <- (P0 - 1) / Tot
          }
        }
      }

      ## c0,c1,c2 from the same fixed-grid LH matrix.
      CMAT <- LH %*% Wproj
      C0[, j] <- CMAT[, 1L]
      C1[, j] <- CMAT[, 2L]
      C2[, j] <- CMAT[, 3L]

      if (isTRUE(need_grad)) {

        DCMAT <- dLH %*% Wproj
        D0[, j] <- DCMAT[, 1L]
        D1[, j] <- DCMAT[, 2L]
        D2[, j] <- DCMAT[, 3L]

        dLH_list[[j]] <- dLH

        if (!is.null(FPB)) {
          FPB_list[[j]] <- FPB
          dC1_dPB[, j] <- as.vector(FPB %*% Wproj[, 2L])
          dC2_dPB[, j] <- as.vector(FPB %*% Wproj[, 3L])
        }
      }

      LH_list[[j]] <- LH
    }

    list(
      C0 = C0,
      C1 = C1,
      C2 = C2,
      D0 = D0,
      D1 = D1,
      D2 = D2,
      LH = LH_list,
      dLH = dLH_list,
      FPB = FPB_list,
      dC1_dPB = dC1_dPB,
      dC2_dPB = dC2_dPB
    )
  }


  f_int_hermite_fast_batch <- function(
      idx,
      b,
      obj_B,
      Hker,
      need_grad = FALSE,
      obs_weight = NULL,
      PB_cache = NULL,
      count_cache = NULL
  ) {

    if (!exists("hermite_interval_kernel_cpp", mode = "function", inherits = TRUE) ||
        !exists("hermite_grid_L0_cpp", mode = "function", inherits = TRUE) ||
        !exists("hermite_count_cache_cpp", mode = "function", inherits = TRUE) ||
        !exists("hermite_state_cpp", mode = "function", inherits = TRUE) ||
        !exists("hermite_reweight_cpp", mode = "function", inherits = TRUE) ||
        !exists("hermite_finish_cpp", mode = "function", inherits = TRUE)) {
      stop(
        paste0(
          "Hermite kernels are not fully loaded. Compile both helper files first: ",
          "Rcpp::sourceCpp('01_hermite_interval.cpp'); ",
          "Rcpp::sourceCpp('02_hermite_polar.cpp')."
        )
      )
    }

    nb <- length(idx)
    QH <- length(herm_x)

    cc <-
      make_hermite_grid_block(
        idx = idx,
        b = b,
        obj_B = obj_B,
        Hker = Hker,
        need_grad = isTRUE(need_grad),
        PB_cache = PB_cache,
        count_cache = count_cache
      )

    st <-
      hermite_state_cpp(
        Psi = Hker$Psi,
        C0 = cc$C0,
        C1 = cc$C1,
        C2 = cc$C2,
        n_threads = approx.hermite.threads
      )

    MUX <- st$mux
    VARX <- st$varx

    ## ---------------------------------------------------------------
    ## Scale-optimized Hermite polynomial correction.
    ## All normalized density-ratio reweighting and its analytic
    ## derivatives are computed in one C++ kernel.  No n_b x Q
    ## temporary weight/score matrices are allocated in R.
    ## ---------------------------------------------------------------
    rw <- hermite_reweight_cpp(
      LH_list = cc$LH,
      dLH_list = if (isTRUE(need_grad)) cc$dLH else vector("list", k),
      FPB_list = if (isTRUE(need_grad) && estimate_B_factor) cc$FPB else vector("list", k),
      MUX = MUX,
      VARX = VARX,
      logbase = st$logbase,
      C0 = cc$C0,
      C1 = cc$C1,
      C2 = cc$C2,
      D0 = cc$D0,
      D1 = cc$D1,
      D2 = cc$D2,
      herm_x = herm_x,
      herm_w = herm_w,
      need_grad = isTRUE(need_grad),
      estimate_B_factor = isTRUE(estimate_B_factor),
      n_threads = as.integer(approx.hermite.threads)
    )

    CORR <- rw$CORR
    CORR_MU <- rw$CORR_MU
    CORR_VAR <- rw$CORR_VAR
    CORR_R <- rw$CORR_R
    GPB_DIRECT <- rw$GPB_DIRECT

    if (is.null(obs_weight)) {
      obs_weight <- rep(1, nb)
    }

    ans <-
      hermite_finish_cpp(
        Psi = Hker$Psi,
        C0 = cc$C0,
        C1 = cc$C1,
        C2 = cc$C2,
        D0 = cc$D0,
        D1 = cc$D1,
        D2 = cc$D2,
        logbase = st$logbase,
        MUX = MUX,
        VARX = VARX,
        MU = st$mu,
        VPACK = st$Vpack,
        CORR = CORR,
        CORR_MU = CORR_MU,
        CORR_VAR = CORR_VAR,
        CORR_R = CORR_R,
        obs_weight = as.numeric(obs_weight),
        need_grad = isTRUE(need_grad),
        n_threads = approx.hermite.threads
      )

    ## Raw-mass diagnostics are produced inside the C++ reweighting kernel.
    ans$h2_logmass_max_abs <- rw$h2_logmass_max_abs
    ans$h2_logmass_mean_abs <- rw$h2_logmass_mean_abs
    ans$h2_mass_max_error <- rw$h2_mass_max_error
    ans$h2_mass_mean_error <- rw$h2_mass_mean_error

    if (!isTRUE(need_grad) || !estimate_B_factor) {
      ans$grad_PB <- NULL
      return(ans)
    }

    ## Analytic derivative wrt observation-specific PB.  c0 cancels
    ## exactly between the Gaussian core and residual correction.
    ans$grad_PB <-
      GPB_DIRECT +
      ans$adj_c1 * cc$dC1_dPB +
      ans$adj_c2 * cc$dC2_dPB

    ans
  }


  ## ============================================================
  ## Joint Hermite lower bound for arbitrary Gaussian zero-factor
  ## dimension L >= 1.
  ##
  ## For ZIPM margins define x = u_j' Z and s = v_j' W.  Even when
  ## W is L-dimensional, each margin depends on W only through the
  ## scalar projection s.  Zero observations therefore require only
  ## a bivariate Hermite projection in (x,s).  The quadratic core is
  ## assembled in the joint S=(Z,W) space of dimension J+L and is
  ## integrated analytically.  The residual first-cumulant correction
  ## remains at most 2-D margin by margin.
  ##
  ## This backend is used only for:
  ##   latent.dependence = "gaussian",
  ##   zero.dependence   = "gaussian",
  ##   latent.loglike    = "hermite",
  ##   L >= 1.
  ## The existing L=0 code path is untouched.  BiExp zero factors keep
  ## the legacy quadrature path because their joint prior is non-Gaussian.
  ## ============================================================

  .ensure_joint_hermite_cpp <- function() {
    if (!exists("joint_hermite_cpp", mode = "function", inherits = TRUE)) {
      stop("Hermite C++ kernel is not loaded. Run Rcpp::sourceCpp('joint_hermite.cpp') first.")
    }
    invisible(TRUE)
  }

  ## ============================================================
  ## L=1 architecture (nested by construction)
  ##
  ## IMPORTANT: for a single structural-zero factor we deliberately DO NOT
  ## use the joint (J+L)-dimensional Hermite surrogate.  Instead L=1 uses
  ## the same count-layer Hermite approximation as L=0 and integrates only
  ## the one-dimensional zero factor outside that count approximation.
  ##
  ## Consequently, at Rho = 0 the conditional zero probabilities are
  ## constant over the zero-factor nodes and the quadrature weights sum to
  ## one, so the L=1 objective reduces numerically to the L=0 objective for
  ## the same Psi (up to floating-point/quadrature normalization error).
  ## This makes L=0 genuinely nested in L=1.
  ##
  ## The joint-Hermite backend is retained only for a possible future L>=2
  ## implementation, where direct zero-factor quadrature would become
  ## multidimensional.
  ## ============================================================
  ## One joint-Hermite C++ route for every Gaussian structural-zero factor model.
  use_nested_hermite_L1 <- FALSE

  use_joint_hermite_allL <-
    isTRUE(estimate_B_factor) &&
    L >= 1L &&
    identical(latent.dependence, "gaussian") &&
    identical(zero.dependence, "gaussian") &&
    identical(latent.loglike, "hermite")


  .joint_hermite_block <- function(idx, pobj, obj_B, Hker, need_grad = FALSE, count_cache = NULL) {
    ## One stable C++ entry point; no version-dependent dispatch in R.
    .ensure_joint_hermite_cpp()
    cc <- if (is.null(count_cache)) {
      make_hermite_count_cache_block(idx, Hker, need_grad = isTRUE(need_grad))
    } else {
      count_cache
    }

    ## The C++ kernel expects one list element per margin.  NULL elements
    ## are retained for PM margins and are never dereferenced there.
    qpb_b <- vector("list", k)
    y0_b <- vector("list", k)
    for (j in seq_len(k)) {
      if (zero_inflation[j]) {
        qpb_b[[j]] <- qpb_list[[j]][idx]
        y0_b[[j]] <- y0_list[[j]][idx]
      } else {
        qpb_b[[j]] <- numeric(0)
        y0_b[[j]] <- logical(0)
      }
    }

    joint_hermite_cpp(
      logY_list = cc$logY,
      dlogY_dr_list = if (isTRUE(need_grad)) cc$dlogY else vector("list", k),
      dlogY_dx_list = if (isTRUE(need_grad)) cc$dlogY_dx else vector("list", k),
      log0_list = cc$log0,
      dlog0_dr_list = if (isTRUE(need_grad)) cc$dlog0 else vector("list", k),
      dlog0_dx_list = if (isTRUE(need_grad)) cc$dlog0_dx else vector("list", k),
      y0_list = y0_b,
      qpb_list = qpb_b,
      Psi = pobj$Psi,
      Rho = obj_B$Rho,
      zi_position = as.integer(zi_position),
      herm_x = herm_x,
      herm_w = herm_w,
      need_grad = isTRUE(need_grad),
      tiny = tiny,
      n_threads = as.integer(approx.hermite.threads)
    )
  }


  .joint_hermite_loadings_valid <- function(Psi, Rho, eps = 1e-10) {
    rn_psi <- sqrt(rowSums(Psi^2))
    rn_rho <- sqrt(rowSums(Rho^2))

    ## Zero loading is a valid independence boundary.
    ## Only reject non-finite values or norms at/above one.
    all(is.finite(rn_psi)) &&
      all(is.finite(rn_rho)) &&
      all(rn_psi >= 0) &&
      all(rn_rho >= 0) &&
      all(rn_psi < 1 - eps) &&
      all(rn_rho < 1 - eps)
  }

  negloglik_copula_hermite_joint_allL <- function(par) {
    rho0 <- par[seq_len(npar_rho)]
    psi0 <- par[npar_rho + seq_len(npar_psi)]

    pobj <- try(make_factor_loadings(
      par = psi0, k = k, q = J, return_all = TRUE
    ), silent = TRUE)
    if (inherits(pobj, "try-error")) return(big_penalty)

    obj_B <- try(compute_B_kernel(rho0, need_grad = FALSE), silent = TRUE)
    if (inherits(obj_B, "try-error")) return(big_penalty)

    if (!.joint_hermite_loadings_valid(pobj$Psi, obj_B$Rho)) {
      return(big_penalty)
    }

    Hker <- try(compute_H_kernel(pobj$Psi), silent = TRUE)
    if (inherits(Hker, "try-error")) return(big_penalty)

    nll <- 0
    for (idx in obs_blocks) {
      z <- try(.joint_hermite_block(
        idx = idx, pobj = pobj, obj_B = obj_B, Hker = Hker, need_grad = FALSE
      ), silent = TRUE)
      if (inherits(z, "try-error")) return(big_penalty)
      lv <- z$log_value

      ## v16: `lv` is a Hermite surrogate, not an exact log probability.
      ## Therefore do NOT impose lv <= 0 observation by observation.
      ## Such a constraint creates artificial barriers and can freeze nlminb,
      ## especially when J or L increases.  Only numerical failures are handled.

      if (any(!is.finite(lv))) {
        ## A few observations can still fail only in the Hermite surrogate.
        ## Do NOT assign an artificial -1e6 contribution: that creates
        ## million-scale jumps in the objective.  Replace only those failed
        ## observation contributions by a typical finite contribution from
        ## the same block.  Final model likelihood is evaluated separately
        ## by RQMC from the original model.
        ok <- is.finite(lv)
        if (!any(ok)) return(big_penalty)
        lv_fill <- stats::median(lv[ok])
        if (!is.finite(lv_fill)) lv_fill <- mean(lv[ok])
        lv[!ok] <- lv_fill
      }
      nll <- nll - sum(lv)
    }
    if (is.finite(nll)) nll else big_penalty
  }

  grad_negloglik_copula_hermite_joint_allL <- function(par) {
    rho0 <- par[seq_len(npar_rho)]
    psi0 <- par[npar_rho + seq_len(npar_psi)]

    pobj <- try(
      make_factor_loadings(
        par = psi0, k = k, q = J, return_all = TRUE
      ),
      silent = TRUE
    )
    if (inherits(pobj, "try-error")) {
      return(rep(0, length(par)))
    }

    obj_B <- try(
      compute_B_kernel(rho0, need_grad = FALSE),
      silent = TRUE
    )
    if (inherits(obj_B, "try-error")) {
      return(rep(0, length(par)))
    }

    ## Important: the joint-Hermite implementation uses
    ## rho_s = t_s v_s, t_s = ||rho_s||.  At t_s = 0 the
    ## representation v_s = rho_s/t_s is undefined even though
    ## the original Gaussian factor model itself is valid.
    ##
    ## nlminb may evaluate the gradient at rejected line-search
    ## trial points.  Do not let such a trial point abort the R
    ## session; the objective returns big_penalty there.
    if (!.joint_hermite_loadings_valid(pobj$Psi, obj_B$Rho)) {
      return(rep(0, length(par)))
    }

    Hker <- try(compute_H_kernel(pobj$Psi), silent = TRUE)
    if (inherits(Hker, "try-error")) {
      return(rep(0, length(par)))
    }

    grad_Psi <- matrix(0, nrow = k, ncol = J)
    grad_Rho <- matrix(0, nrow = n_zi, ncol = L)

    for (idx in obs_blocks) {
      z <- try(
        .joint_hermite_block(
          idx = idx, pobj = pobj, obj_B = obj_B, Hker = Hker, need_grad = TRUE
        ),
        silent = TRUE
      )
      if (inherits(z, "try-error")) {
        return(rep(0, length(par)))
      }

      ## v16: never zero the ENTIRE gradient because one observation has a
      ## non-finite surrogate value.  The objective replaces only failed lv
      ## contributions by a typical finite block value; here we retain every
      ## finite analytic gradient component and zero only components that are
      ## themselves non-finite.  In particular, positive finite surrogate
      ## log-values are allowed during optimization.
      z$grad_Psi[!is.finite(z$grad_Psi)] <- 0
      z$grad_Rho[!is.finite(z$grad_Rho)] <- 0
      ## C++ returns gradient of the surrogate log likelihood.
      ## The optimizer minimizes negative log likelihood.
      grad_Psi <- grad_Psi - z$grad_Psi
      grad_Rho <- grad_Rho - z$grad_Rho
    }

    grad_psi_raw <- chain_grad_loading_to_raw(
      grad_Load = grad_Psi,
      Araw = pobj$A,
      map = pobj$map,
      npar = npar_psi,
      q = J
    )
    grad_rho_raw <- chain_grad_loading_to_raw(
      grad_Load = grad_Rho,
      Araw = obj_B$A,
      map = obj_B$map,
      npar = npar_rho,
      q = L
    )
    c(grad_rho_raw, grad_psi_raw)
  }


  negloglik_copula_hermite <- function(par) {

    if (estimate_B_factor) {

      rho0 <- par[seq_len(npar_rho)]
      psi0 <- par[npar_rho + seq_len(npar_psi)]

      obj_B <-
        try(
          compute_B_kernel(
            rho0,
            need_grad = FALSE
          ),
          silent = TRUE
        )

      if (inherits(obj_B, "try-error")) return(big_penalty)

      n_node_B_local <- n_node_B
      Wnode_B_local <- Wnode_B

    } else {

      obj_B <- NULL
      psi0 <- par[seq_len(npar_psi)]
      n_node_B_local <- 1L
      Wnode_B_local <- 1
    }

    pobj <-
      try(
        make_factor_loadings(
          par = psi0,
          k = k,
          q = J,
          return_all = TRUE
        ),
        silent = TRUE
      )

    if (inherits(pobj, "try-error")) return(big_penalty)

    Hker <- try(compute_H_kernel(pobj$Psi), silent = TRUE)
    if (inherits(Hker, "try-error")) return(big_penalty)

    nll <- 0

    for (idx in obs_blocks) {

      nb <- length(idx)
      PB_cache <- if (estimate_B_factor && L == 1L) make_PB_cache_L1(idx, obj_B, need_grad = FALSE) else NULL
      count_cache <- make_hermite_count_cache_block(idx, Hker, need_grad = FALSE)
      LOG <- matrix(-Inf, nrow = nb, ncol = n_node_B_local)

      for (b in seq_len(n_node_B_local)) {

        hh <-
          try(
            f_int_hermite_fast_batch(
              idx = idx,
              b = b,
              obj_B = obj_B,
              Hker = Hker,
              need_grad = FALSE,
              PB_cache = PB_cache,
              count_cache = count_cache
            ),
            silent = TRUE
          )

        if (inherits(hh, "try-error")) return(big_penalty)

        good <- is.finite(hh$log_value) & Wnode_B_local[b] > 0

        if (any(good)) {
          LOG[good, b] <-
            log(Wnode_B_local[b]) +
            hh$log_value[good]
        }
      }

      row_max <- apply(LOG, 1L, max)
      if (any(!is.finite(row_max))) return(big_penalty)

      logLi <- row_max + log(rowSums(exp(LOG - row_max)))
      if (any(!is.finite(logLi))) return(big_penalty)

      nll <- nll - sum(logLi)
    }

    if (is.finite(nll)) nll else big_penalty
  }


  grad_negloglik_copula_hermite <- function(par) {

    if (estimate_B_factor) {

      rho0 <- par[seq_len(npar_rho)]
      psi0 <- par[npar_rho + seq_len(npar_psi)]

      ## need_grad=TRUE is essential for the BiExp zero layer because
      ## dQ/dRho enters dPB/dRho.
      obj_B <-
        compute_B_kernel(
          rho0,
          need_grad = TRUE
        )

      n_node_B_local <- n_node_B
      Wnode_B_local <- Wnode_B

    } else {

      obj_B <- NULL
      psi0 <- par[seq_len(npar_psi)]
      n_node_B_local <- 1L
      Wnode_B_local <- 1
    }

    pobj <-
      make_factor_loadings(
        par = psi0,
        k = k,
        q = J,
        return_all = TRUE
      )

    Hker <- compute_H_kernel(pobj$Psi)

    grad_Psi <- matrix(0, nrow = k, ncol = J)

    grad_Rho <- if (estimate_B_factor) {
      matrix(0, nrow = n_zi, ncol = L)
    } else {
      NULL
    }

    for (idx in obs_blocks) {

      nb <- length(idx)
      PB_cache <- if (estimate_B_factor && L == 1L) make_PB_cache_L1(idx, obj_B, need_grad = TRUE) else NULL
      count_cache <- make_hermite_count_cache_block(idx, Hker, need_grad = TRUE)

      if (!estimate_B_factor) {

        hh <-
          f_int_hermite_fast_batch(
            idx = idx,
            b = 1L,
            obj_B = NULL,
            Hker = Hker,
            need_grad = TRUE,
            obs_weight = rep(1, nb),
            count_cache = count_cache
          )

        grad_Psi <- grad_Psi - hh$grad_Psi
        next
      }

      ## First pass: posterior weights over zero-factor nodes.
      LOG <- matrix(-Inf, nrow = nb, ncol = n_node_B_local)

      for (b in seq_len(n_node_B_local)) {

        hh <-
          f_int_hermite_fast_batch(
            idx = idx,
            b = b,
            obj_B = obj_B,
            Hker = Hker,
            need_grad = FALSE,
            PB_cache = PB_cache,
            count_cache = count_cache
          )

        good <- is.finite(hh$log_value) & Wnode_B_local[b] > 0

        if (any(good)) {
          LOG[good, b] <-
            log(Wnode_B_local[b]) +
            hh$log_value[good]
        }
      }

      row_max <- apply(LOG, 1L, max)
      if (any(!is.finite(row_max))) {
        return(rep(0, length(par)))
      }

      WW <- exp(LOG - row_max)
      den <- rowSums(WW)

      if (any(!is.finite(den)) || any(den <= 0)) {
        return(rep(0, length(par)))
      }

      WW <- WW / den

      ## Second pass: analytic Psi and Rho gradients.
      for (b in seq_len(n_node_B_local)) {

        wb <- WW[, b]

        if (!any(is.finite(wb) & wb != 0)) next

        hh <-
          f_int_hermite_fast_batch(
            idx = idx,
            b = b,
            obj_B = obj_B,
            Hker = Hker,
            need_grad = TRUE,
            obs_weight = wb,
            PB_cache = PB_cache,
            count_cache = count_cache
          )

        grad_Psi <- grad_Psi - hh$grad_Psi

        for (s in seq_len(n_zi)) {

          j <- zi_index[s]
          gpb <- hh$grad_PB[, j]

          if (!any(is.finite(gpb) & gpb != 0)) next

          for (l in seq_len(L)) {

            dPB <- if (L == 1L) {
              get_dPB_cached(j = j, b = b, idx = idx, obj_B = obj_B, cache = PB_cache)
            } else {
              get_dPB_block(j = j, b = b, idx = idx, obj_B = obj_B, l = l)
            }

            grad_Rho[s, l] <-
              grad_Rho[s, l] -
              sum(wb * gpb * dPB)
          }
        }
      }
    }

    grad_psi_raw <-
      chain_grad_loading_to_raw(
        grad_Load = grad_Psi,
        Araw = pobj$A,
        map = pobj$map,
        npar = npar_psi,
        q = J
      )

    if (!estimate_B_factor) {
      return(grad_psi_raw)
    }

    grad_rho_raw <-
      chain_grad_loading_to_raw(
        grad_Load = grad_Rho,
        Araw = obj_B$A,
        map = obj_B$map,
        npar = npar_rho,
        q = L
      )

    c(
      grad_rho_raw,
      grad_psi_raw
    )
  }




  ## Final RQMC evaluation for a fitted latent Gaussian factor.
  ## This is called only after Price optimization has converged.
  ## The zero layer remains integrated by its existing GH/GL rule.
  ## =========================================================


  ## ---------------------------------------------------------
  ## Final RQMC for Gaussian count + Gaussian zero factors.
  ## Integrates the exact likelihood jointly in dimension J+L;
  ## there is no GH tensor over the zero-factor layer.
  ## ---------------------------------------------------------
  qmc_loglik_joint_gaussian_allL <- function(psi0, rho0, Q_eval = qmc.Q, n_rep_eval = qmc.n_rep) {
    if (!requireNamespace("randtoolbox", quietly = TRUE)) {
      stop("Package 'randtoolbox' is required for final QMC evaluation.")
    }
    pobj <- make_factor_loadings(par = psi0, k = k, q = J, return_all = TRUE)
    Psi <- pobj$Psi
    sN <- sqrt(pmax(1e-14, 1 - rowSums(Psi^2)))
    bobj <- compute_B_kernel(rho0, need_grad = FALSE)
    Rho <- bobj$Rho
    sB <- sqrt(pmax(1e-14, 1 - rowSums(Rho^2)))
    Dq <- J + L

    U0 <- randtoolbox::sobol(
      n = Q_eval, dim = Dq, scrambling = 0,
      seed = qmc.seed, normal = FALSE
    )
    U0 <- as.matrix(U0)
    set.seed(qmc.seed)
    loglik_rep <- rep(NA_real_, n_rep_eval)
    Li_rep <- matrix(NA_real_, nrow = n, ncol = n_rep_eval)

    for (rr in seq_len(n_rep_eval)) {
      shift <- runif(Dq)
      U <- sweep(U0, 2L, shift, "+") %% 1
      U <- pmin(pmax(U, 1e-15), 1 - 1e-15)
      Li <- numeric(n)

      for (qs in seq(1L, Q_eval, by = qmc.chunk_size)) {
        qe <- min(qs + qmc.chunk_size - 1L, Q_eval)
        S <- qnorm(U[qs:qe, , drop = FALSE])
        Zq <- S[, seq_len(J), drop = FALSE]
        Wq <- S[, J + seq_len(L), drop = FALSE]
        nq <- nrow(S)

        etaN <- Zq %*% t(Psi)
        etaB <- Wq %*% t(Rho)
        OmegaT_list_q <- vector("list", k)
        for (j in seq_len(k)) {
          A_inner <- sweep(
            matrix(t_inner, nrow = nq, ncol = d0 - 1L, byrow = TRUE),
            1L, etaN[, j], "-"
          )
          q_inner <- pnorm(A_inner / sN[j])
          q_full <- cbind(0, q_inner, 1)
          Omega <- q_full[, 2:(d0 + 1L), drop = FALSE] -
            q_full[, 1:d0, drop = FALSE]
          OmegaT_list_q[[j]] <- t(Omega)
        }

        for (idx in obs_blocks) {
          nb <- length(idx)
          Tprod <- matrix(1, nrow = nb, ncol = nq)
          for (j in seq_len(k)) {
            PNy <- d0 * (I_y_list[[j]][idx, , drop = FALSE] %*% OmegaT_list_q[[j]])
            if (!zero_inflation[j]) {
              Tj <- PNy
            } else {
              s <- zi_position[j]
              qpb <- qpb_list[[j]][idx]
              Aval <- sweep(
                matrix(-etaB[, s], nrow = nb, ncol = nq, byrow = TRUE),
                1L, qpb, "+"
              ) / sB[s]
              PB <- pnorm(Aval)
              y0 <- y0_list[[j]][idx]
              Tj <- matrix(0, nrow = nb, ncol = nq)
              if (any(y0)) {
                PN0 <- d0 * (I_0_list[[j]][idx[y0], , drop = FALSE] %*% OmegaT_list_q[[j]])
                Tj[y0, ] <- (1 - PB[y0, , drop = FALSE]) +
                  PB[y0, , drop = FALSE] * PN0
              }
              if (any(!y0)) {
                Tj[!y0, ] <- PB[!y0, , drop = FALSE] * PNy[!y0, , drop = FALSE]
              }
            }
            Tprod <- Tprod * pmax(Tj, tiny)
          }
          Li[idx] <- Li[idx] + rowSums(Tprod) / Q_eval
        }
      }

      if (all(is.finite(Li)) && all(Li > tiny)) {
        Li_rep[, rr] <- Li
        loglik_rep[rr] <- sum(log(Li))
      }
    }

    good <- is.finite(loglik_rep)
    if (!any(good)) stop("All final joint QMC likelihood replicates failed.")
    Li_mean <- rowMeans(Li_rep[, good, drop = FALSE])
    if (any(!is.finite(Li_mean)) || any(Li_mean <= tiny)) {
      stop("Invalid observation-level mean joint QMC likelihood.")
    }
    list(
      log_lik_value = sum(log(Li_mean)),
      log_lik_se = if (sum(good) > 1L) sd(loglik_rep[good]) / sqrt(sum(good)) else NA_real_,
      replicate_log_lik = loglik_rep,
      Q = Q_eval,
      n_rep = n_rep_eval,
      seed = qmc.seed,
      dimension = Dq,
      method = "joint Gaussian RQMC"
    )
  }

  ## =========================================================
  ## Genz-style sequential conditional RQMC for Gaussian factor
  ## PM / ZIPM likelihoods.
  ##
  ## New C++ backend: genz_conditional_zipm_qmc_cpp().
  ##   * L = 0: PM and independent-zero ZIPM margins.
  ##   * L = 1: Gaussian zero-factor ZIPM margins.
  ##
  ## Parallelization is performed inside C++ over observations.
  ## U_count drives the sequential count-factor truncation and
  ## U_zero drives the structural-zero branch / zero-factor update.
  ## =========================================================
  genz_conditional_qmc_loglik_at_fit <- function(psi0, rho0 = numeric(0)) {
    if (!requireNamespace("randtoolbox", quietly = TRUE)) {
      stop("Package 'randtoolbox' is required.")
    }
    if (!exists("genz_conditional_zipm_qmc_cpp", mode = "function", inherits = TRUE)) {
      stop(
        paste0(
          "C++ helper `genz_conditional_zipm_qmc_cpp()` is not loaded. ",
          "Compile 05_genz_final.cpp (or run Rcpp::compileAttributes() + devtools::load_all())."
        )
      )
    }
    if (!identical(latent.dependence, "gaussian")) {
      stop("Genz conditional RQMC requires Gaussian count dependence.")
    }
    if (estimate_B_factor && !identical(zero.dependence, "gaussian")) {
      stop("The new Genz zero-factor path currently requires zero.dependence = 'gaussian'.")
    }
    if (estimate_B_factor && L != 1L) {
      stop("The new Genz zero-factor path currently supports exactly L = 1.")
    }

    pobj <- make_factor_loadings(par = psi0, k = k, q = J, return_all = TRUE)
    Psi <- pobj$Psi

    if (estimate_B_factor) {
      bobj <- make_factor_loadings(par = rho0, k = n_zi, q = L, return_all = TRUE)
      Rho <- bobj$Psi
    } else {
      ## A 0-column matrix explicitly selects the L=0 branch in C++.
      Rho <- matrix(numeric(0), nrow = n_zi, ncol = 0L)
    }

    ## The C++ interface requires one element per margin.  Empty vectors
    ## are used where a quantity is irrelevant and are never dereferenced.
    I0_genz <- vector("list", k)
    pb_genz <- vector("list", k)
    qpb_genz <- vector("list", k)
    y0_genz <- vector("list", k)
    for (j in seq_len(k)) {
      if (zero_inflation[j]) {
        I0_genz[[j]] <- I_0_list[[j]]
        pb_genz[[j]] <- pb_list[[j]]
        qpb_genz[[j]] <- if (estimate_B_factor) qpb_list[[j]] else numeric(0)
        y0_genz[[j]] <- y0_list[[j]]
      } else {
        I0_genz[[j]] <- matrix(numeric(0), nrow = 0L, ncol = 0L)
        pb_genz[[j]] <- numeric(0)
        qpb_genz[[j]] <- numeric(0)
        y0_genz[[j]] <- logical(0)
      }
    }

    Quse <- max(2L, as.integer(genz.conditional.Q))
    ant <- isTRUE(genz.conditional.antithetic)
    Qbase <- if (ant) as.integer(ceiling(Quse / 2)) else Quse

    ## Use one 2k-dimensional Sobol net and split it into independent roles.
    ## This is preferable to reusing the same coordinates for count and zero
    ## branch decisions while retaining low-discrepancy structure jointly.
    U0 <- randtoolbox::sobol(
      n = Qbase, dim = 2L * k, scrambling = 0,
      seed = qmc.seed, normal = FALSE
    )
    U0 <- as.matrix(U0)

    set.seed(qmc.seed)
    loglik_rep <- rep(NA_real_, qmc.n_rep)
    logLi_rep <- matrix(NA_real_, nrow = n, ncol = qmc.n_rep)
    nth <- max(1L, as.integer(approx.hermite.threads))

    .internal_msg(paste0(
      "[dependence] Genz conditional RQMC | k=", k,
      " | count factors J=", J,
      " | zero factors L=", if (estimate_B_factor) L else 0L,
      " | Q=", Quse,
      " | antithetic=", ant,
      " | rep=", qmc.n_rep,
      " | threads=", nth
    ))

    for (rr in seq_len(qmc.n_rep)) {
      shift <- runif(2L * k)
      Uall <- sweep(U0, 2L, shift, "+") %% 1
      Uall <- pmin(pmax(Uall, 1e-15), 1 - 1e-15)

      if (ant) {
        Uall <- rbind(Uall, 1 - Uall)
        if (nrow(Uall) > Quse) {
          Uall <- Uall[seq_len(Quse), , drop = FALSE]
        }
      }

      U_count <- Uall[, seq_len(k), drop = FALSE]
      U_zero <- Uall[, k + seq_len(k), drop = FALSE]

      li <- genz_conditional_zipm_qmc_cpp(
        U_count = U_count,
        U_zero = U_zero,
        Psi = Psi,
        Rho = Rho,
        zi_position = as.integer(zi_position),
        t_inner = t_inner,
        I_y_list = I_y_list,
        I_0_list = I0_genz,
        zero_inflation = zero_inflation,
        pb_list = pb_genz,
        qpb_list = qpb_genz,
        y0_list = y0_genz,
        tiny = tiny,
        n_threads = nth
      )

      if (length(li) == n && all(is.finite(li))) {
        logLi_rep[, rr] <- li
        loglik_rep[rr] <- sum(li)
      }
    }

    good <- which(is.finite(loglik_rep))
    if (!length(good)) stop("All Genz-conditional RQMC replicates failed.")

    ## Average replicate likelihoods observation-wise on the log scale.
    logLi_mean <- numeric(n)
    for (i in seq_len(n)) {
      v <- logLi_rep[i, good]
      vm <- max(v)
      logLi_mean[i] <- vm + log(sum(exp(v - vm))) - log(length(good))
    }

    list(
      log_lik_value = sum(logLi_mean),
      log_lik_se = if (length(good) > 1L) sd(loglik_rep[good]) / sqrt(length(good)) else NA_real_,
      replicate_log_lik = loglik_rep,
      Q = Quse,
      n_rep = qmc.n_rep,
      antithetic = ant,
      sequence_dimension = 2L * k,
      count_factor_dimension = J,
      zero_factor_dimension = if (estimate_B_factor) L else 0L,
      qmc_threads = nth,
      parallel_dimension = "observations",
      method = if (estimate_B_factor) {
        "Genz sequential conditional RQMC with Gaussian zero factor (L=1)"
      } else if (any(zero_inflation)) {
        "Genz sequential conditional RQMC with independent zero inflation (L=0)"
      } else {
        "Genz sequential conditional RQMC for PM (L=0)"
      }
    )
  }

  ## =========================================================
  ## Rotated + truncated + antithetic RQMC for L = 0.
  ##
  ## 1) Orthogonally rotate Psi by eigenvectors of Psi'Psi.
  ## 2) Retain directions explaining `simulation.rotation.energy`, while
  ##    also dropping directions with lambda/lambda_max < simulation.eigen.tol.
  ## 3) Fold discarded factor variance into the idiosyncratic Gaussian term.
  ##    Thus every latent margin remains N(0,1); only weak cross-margin
  ##    dependence directions are removed.
  ## 4) Use antithetic Sobol pairs Z and -Z in the retained dimension.
  ## =========================================================
  ## =========================================================
  ## 8a-bis. Grouped tilted RQMC importance sampling for L = 0
  ##
  ## The Hermite optimizer/lower bound is UNCHANGED.  This routine uses
  ## only the quadratic Hermite core after fitting to construct Gaussian
  ## proposals.  Observations with similar Hermite tilt means are grouped,
  ## so every group shares one proposal q_g = N(mu_g, V_g).  This preserves
  ## the efficient Q-parallel likelihood kernel while moving simulation
  ## points toward high-posterior-mass regions when J is large.
  ##
  ## Exact importance correction is applied node by node:
  ##   log w_q = log phi_J(z_q) - log q_g(z_q).
  ## Hence the target remains the ORIGINAL model likelihood; Hermite is used
  ## only to design the proposal.
  ## =========================================================
  ## =========================================================
  ## 8b. Specialized exact GH backend for L = 1
  ##
  ## IMPORTANT: this fast path is completely isolated from L = 0.
  ## It keeps count GH (q) and zero GH (b) as separate axes and
  ## contracts them in C++.  No nb x (Q_N * Q_B) R matrices are
  ## materialized.  gh_order_B is unchanged.
  ## =========================================================

  use_L1_joint_GH <- isTRUE(estimate_B_factor) &&
    identical(as.integer(L), 1L) &&
    latent.dependence == "gaussian" &&
    zero.dependence == "gaussian" &&
    latent.loglike == "GH"

  ## L=1 exact GH uses the package-level C++ backend for both serial and
  ## parallel execution.  gh_L1_eval_cpp_cached() accepts n_threads = 1,
  ## so the former runtime-compiled legacy fallback is unnecessary.
  if (use_L1_joint_GH) {
    if (!exists("gh_L1_eval_cpp_cached", mode = "function", inherits = TRUE)) {
      stop(
        "C++ L=1 GH backend `gh_L1_eval_cpp_cached()` is not loaded. ",
        "Run Rcpp::compileAttributes() and devtools::load_all()."
      )
    }
    .internal_msg(
      "[dependence] L=1 GH backend: package-level C++ observation blocks | threads=",
      as.integer(gh.threads),
      " | count nodes=", length(Wnode_N),
      " | zero nodes=", length(Wnode_B)
    )
  }

  negloglik_copula_GH_L1_joint <- function(par) {
    rho0 <- par[seq_len(npar_rho)]
    psi0 <- par[npar_rho + seq_len(npar_psi)]

    obj_B <- try(compute_B_kernel(rho0, need_grad = FALSE), silent = TRUE)
    if (inherits(obj_B, "try-error")) return(big_penalty)
    pobj <- try(make_factor_loadings(
      par = psi0, k = k, q = J, return_all = TRUE
    ), silent = TRUE)
    if (inherits(pobj, "try-error")) return(big_penalty)

    ans <- try(gh_L1_eval_cpp_cached(
      I_y_list = I_y_list,
      I_0_list = I_0_list,
      qpb_list = qpb_list[zi_index],
      y0_list = y0_list,
      zi_pos = as.integer(zi_position),
      Psi = pobj$Psi,
      Rho = as.vector(obj_B$Rho[, 1L]),
      Znode_N = Znode_N,
      Wnode_N = Wnode_N,
      Znode_B = as.vector(Znode_B[, 1L]),
      Wnode_B = Wnode_B,
      t_inner = t_inner,
      d0 = as.integer(d0),
      block_size = as.integer(copula.block_size),
      need_grad = FALSE,
      tiny = tiny,
      n_threads = as.integer(gh.threads)
    ), silent = TRUE)

    if (inherits(ans, "try-error") || !is.finite(ans$nll)) return(big_penalty)
    ans$nll
  }

  grad_negloglik_copula_GH_L1_joint <- function(par) {
    rho0 <- par[seq_len(npar_rho)]
    psi0 <- par[npar_rho + seq_len(npar_psi)]

    obj_B <- compute_B_kernel(rho0, need_grad = FALSE)
    pobj <- make_factor_loadings(
      par = psi0, k = k, q = J, return_all = TRUE
    )

    ans <- gh_L1_eval_cpp_cached(
      I_y_list = I_y_list,
      I_0_list = I_0_list,
      qpb_list = qpb_list[zi_index],
      y0_list = y0_list,
      zi_pos = as.integer(zi_position),
      Psi = pobj$Psi,
      Rho = as.vector(obj_B$Rho[, 1L]),
      Znode_N = Znode_N,
      Wnode_N = Wnode_N,
      Znode_B = as.vector(Znode_B[, 1L]),
      Wnode_B = Wnode_B,
      t_inner = t_inner,
      d0 = as.integer(d0),
      block_size = as.integer(copula.block_size),
      need_grad = TRUE,
      tiny = tiny,
      n_threads = as.integer(gh.threads)
    )

    grad_psi0 <- chain_grad_loading_to_raw(
      grad_Load = ans$grad_Psi, Araw = pobj$A, map = pobj$map,
      npar = npar_psi, q = J)
    grad_rho0 <- chain_grad_loading_to_raw(
      grad_Load = matrix(ans$grad_Rho, ncol = 1L),
      Araw = obj_B$A, map = obj_B$map,
      npar = npar_rho, q = 1L)
    c(grad_rho0, grad_psi0)
  }

  ## =========================================================
  ## Fast exact GH backend for Gaussian L=0 + ZIP margins.
  ##
  ## Structural-zero indicators are independent when L=0, hence for a
  ## ZIP margin j the node-wise contribution is
  ##
  ##   y > 0 : pb_ij * PN_y(i,q)
  ##   y = 0 : (1-pb_ij) + pb_ij * PN_0(i,q).
  ##
  ## The expensive product over margins / GH nodes and its analytic
  ## loading gradient are fused in C++ and parallelized over observations.
  ## This backend is valid for any J; in particular J = 1,2,3.
  ## =========================================================
  ## Package-level C++ backend: gh_L0_ZIP_block_cpp()
  ## is compiled from src/gh_L0_ZIP.cpp and exposed through RcppExports.

  use_fast_GH_L0_ZIP <-
    identical(latent.dependence, "gaussian") &&
    identical(latent.loglike, "GH") &&
    identical(as.integer(L), 0L) &&
    !estimate_B_factor &&
    any(zero_inflation)

  .fast_gh_L0_ZIP_eval <- function(psi0, need_grad = FALSE) {
    pobj <- make_factor_loadings(par = psi0, k = k, q = J, return_all = TRUE)
    obj_N <- compute_N_kernel(psi0, need_grad = isTRUE(need_grad))
    nll <- 0
    grad_Psi <- if (isTRUE(need_grad)) matrix(0, nrow = k, ncol = J) else NULL
    threads_used <- 1L

    dPN_cols_all <- if (isTRUE(need_grad)) seq_len(n_node_N * J) else integer(0)

    for (idx in obs_blocks) {
      PN <- make_PN_block(idx, obj_N)
      nb <- length(idx)
      pb_b <- lapply(seq_len(k), function(j) {
        if (zero_inflation[j]) as.numeric(pb_list[[j]][idx]) else numeric(0)
      })
      y0_b <- lapply(seq_len(k), function(j) {
        if (zero_inflation[j]) as.logical(y0_list[[j]][idx]) else logical(0)
      })

      if (isTRUE(need_grad)) {
        dPNy <- vector("list", k)
        dPN0 <- vector("list", k)
        for (j in seq_len(k)) {
          Iy <- I_y_list[[j]][idx, , drop = FALSE]
          dOy <- obj_N$dOmegaT_stack_list[[j]]
          dPNy[[j]] <- d0 * (Iy %*% dOy)
          if (zero_inflation[j]) {
            I0 <- I_0_list[[j]][idx, , drop = FALSE]
            dPN0[[j]] <- d0 * (I0 %*% dOy)
          } else {
            dPN0[[j]] <- matrix(0, nrow = nb, ncol = length(dPN_cols_all))
          }
        }
      } else {
        dPNy <- vector("list", k)
        dPN0 <- vector("list", k)
      }

      ans <- gh_L0_ZIP_block_cpp(
        PNy_list = PN$PN_y,
        PN0_list = PN$PN_0,
        dPNy_list = dPNy,
        dPN0_list = dPN0,
        pb_list = pb_b,
        y0_list = y0_b,
        zero_inflation = zero_inflation,
        WN = Wnode_N,
        J = as.integer(J),
        need_grad = isTRUE(need_grad),
        tiny = tiny,
        n_threads = as.integer(gh.threads)
      )
      if (!is.finite(ans$nll)) return(list(nll = Inf, grad_raw = NULL))
      nll <- nll + ans$nll
      threads_used <- ans$threads_used
      if (isTRUE(need_grad)) grad_Psi <- grad_Psi + ans$grad_Psi
    }

    grad_raw <- if (isTRUE(need_grad)) {
      chain_grad_loading_to_raw(
        grad_Load = grad_Psi, Araw = pobj$A, map = pobj$map,
        npar = npar_psi, q = J
      )
    } else NULL

    list(nll = nll, grad_Psi = grad_Psi, grad_raw = grad_raw,
         threads_used = threads_used)
  }

  ## =========================================================
  ## Fast exact GH backend for Gaussian L=0 + all-PM margins.
  ## Uses a C++ two-pass log-domain evaluator and analytic gradient.
  ## Other cases fall back to the legacy R implementation below.
  ## =========================================================
  use_fast_GH_PM_L0 <-
    identical(latent.dependence, "gaussian") &&
    identical(latent.loglike, "GH") &&
    !estimate_B_factor &&
    all(!zero_inflation) &&
    exists("gh_pm_L0_eval_cpp_cached", mode = "function", inherits = TRUE)

  .fast_gh_eval <- function(psi0, need_grad = FALSE) {
    pobj <- make_factor_loadings(
      par = psi0, k = k, q = J, return_all = TRUE
    )
    ans <- gh_pm_L0_eval_cpp_cached(
      I_y_list = I_y_list,
      Psi = pobj$Psi,
      Znode = Znode_N,
      Wnode = Wnode_N,
      t_inner = t_inner,
      d0 = as.integer(d0),
      block_size = as.integer(copula.block_size),
      need_grad = isTRUE(need_grad),
      tiny = tiny,
      n_threads = as.integer(gh.threads)
    )
    if (isTRUE(need_grad)) {
      ans$grad_raw <- chain_grad_loading_to_raw(
        grad_Load = ans$grad_Psi,
        Araw = pobj$A,
        map = pobj$map,
        npar = npar_psi,
        q = J
      )
    } else {
      ans$grad_raw <- NULL
    }
    ans
  }

  if (!isTRUE(.gh_progressive_internal) &&
      identical(latent.dependence, "gaussian") &&
      identical(latent.loglike, "GH")) {
    if (isTRUE(use_fast_GH_L0_ZIP)) {
      .internal_msg(
        "[dependence] L=0 ZIP GH backend: C++ observation-parallel fused contraction | threads=",
        gh.threads, " | J=", J, " | count nodes=", n_node_N
      )
    } else if (isTRUE(use_fast_GH_PM_L0)) {
      .internal_msg(
        "[dependence] fast GH backend: C++ cached-PN scaled-product | PM-only L=0 | threads=",
        gh.threads, " | nodes=", n_node_N
      )
    } else {
      .internal_msg("[dependence] GH backend: legacy R path (fast L0 kernel not applicable or not loaded)")
    }
  }

  ## =========================================================
  ## 9. Negative log-likelihood -- block streaming
  ##
  ## No n x n_node_N matrices are retained across the full
  ## sample.  At most `copula.block_size` observations are
  ## active in memory.
  ## =========================================================

  negloglik_copula <- function(par) {

    if (use_joint_hermite_allL) {
      return(negloglik_copula_hermite_joint_allL(par))
    }
    if (latent.dependence == "gaussian" && latent.loglike == "hermite") {
      return(negloglik_copula_hermite(par))
    }
    if (use_L1_joint_GH) {
      return(negloglik_copula_GH_L1_joint(par))
    }
    if (isTRUE(use_fast_GH_L0_ZIP)) {
      z <- try(.fast_gh_L0_ZIP_eval(par[seq_len(npar_psi)], need_grad = FALSE), silent = TRUE)
      if (inherits(z, "try-error") || !is.finite(z$nll)) return(big_penalty)
      return(z$nll)
    }
    if (isTRUE(use_fast_GH_PM_L0)) {
      z <- try(.fast_gh_eval(par[seq_len(npar_psi)], need_grad = FALSE), silent = TRUE)
      if (inherits(z, "try-error") || !is.finite(z$nll)) return(big_penalty)
      return(z$nll)
    }

    if (estimate_B_factor) {

      rho0 <- par[
        seq_len(npar_rho)
      ]

      psi0 <- par[
        npar_rho +
          seq_len(npar_psi)
      ]

      obj_B <- try(
        compute_B_kernel(
          rho0,
          need_grad = FALSE
        ),
        silent = TRUE
      )

      if (inherits(
        obj_B,
        "try-error"
      )) {
        return(big_penalty)
      }

      n_node_B_local <- n_node_B
      Wnode_B_local <- Wnode_B

    } else {

      obj_B <- NULL

      psi0 <- par[
        seq_len(npar_psi)
      ]

      n_node_B_local <- 1L
      Wnode_B_local <- 1
    }

    obj_N <- try(
      compute_N_kernel(
        psi0,
        need_grad = FALSE
      ),
      silent = TRUE
    )

    if (inherits(
      obj_N,
      "try-error"
    )) {
      return(big_penalty)
    }

    nll <- 0

    for (idx in obs_blocks) {

      PN <- try(
        make_PN_block(
          idx,
          obj_N
        ),
        silent = TRUE
      )

      if (inherits(
        PN,
        "try-error"
      )) {
        return(big_penalty)
      }

      nb <- length(idx)
      Li <- numeric(nb)
      PB_cache <- if (estimate_B_factor && L == 1L) make_PB_cache_L1(idx, obj_B, need_grad = FALSE) else NULL

      for (b in seq_len(
        n_node_B_local
      )) {

        Tprod <- matrix(
          1,
          nrow = nb,
          ncol = n_node_N
        )

        for (j in seq_len(k)) {

          PB_b <- if (zero_inflation[j]) {
            get_PB_cached(j = j, b = b, idx = idx, obj_B = obj_B, cache = PB_cache)
          } else NULL

          Tj <- make_T_block(
            j = j, b = b, idx = idx, PN = PN, obj_B = obj_B, PB_b = PB_b
          )

          if (any(!is.finite(Tj))) {
            return(big_penalty)
          }

          Tprod <-
            Tprod *
            pmax(
              Tj,
              tiny
            )
        }

        Li <-
          Li +
          Wnode_B_local[b] *
          as.vector(
            Tprod %*%
              Wnode_N
          )
      }

      if (any(!is.finite(Li)) ||
          any(Li <= tiny)) {
        return(big_penalty)
      }

      nll <-
        nll -
        sum(
          log(Li)
        )

      # Explicitly release the large block objects before the
      # next block.
      rm(
        PN,
        Li
      )
    }

    if (!is.finite(nll)) {
      big_penalty
    } else {
      nll
    }
  }


  ## =========================================================
  ## 10. Analytic gradient -- block / derivative streaming
  ##
  ## Crucially, we DO NOT store:
  ##
  ##   dPN_y_list : k lists of n x node x J arrays
  ##   dPN_0_list : k lists of n x node x J arrays
  ##   dLi_dPsi   : n x k x J
  ##   dLi_dRho   : n x n_zi x L
  ##
  ## Only the final k x J and n_zi x L gradients are retained.
  ## =========================================================

  grad_negloglik_copula <- function(par) {

    if (use_joint_hermite_allL) {
      return(grad_negloglik_copula_hermite_joint_allL(par))
    }
    if (latent.dependence == "gaussian" && latent.loglike == "hermite") {
      return(grad_negloglik_copula_hermite(par))
    }
    if (use_L1_joint_GH) {
      return(grad_negloglik_copula_GH_L1_joint(par))
    }
    if (isTRUE(use_fast_GH_L0_ZIP)) {
      z <- try(.fast_gh_L0_ZIP_eval(par[seq_len(npar_psi)], need_grad = TRUE), silent = TRUE)
      if (inherits(z, "try-error") || is.null(z$grad_raw) || any(!is.finite(z$grad_raw))) {
        return(rep(0, length(par)))
      }
      return(z$grad_raw)
    }
    if (isTRUE(use_fast_GH_PM_L0)) {
      z <- try(.fast_gh_eval(par[seq_len(npar_psi)], need_grad = TRUE), silent = TRUE)
      if (inherits(z, "try-error") || is.null(z$grad_raw) || any(!is.finite(z$grad_raw))) {
        return(rep(0, length(par)))
      }
      return(z$grad_raw)
    }

    if (estimate_B_factor) {

      rho0 <- par[seq_len(npar_rho)]
      psi0 <- par[npar_rho + seq_len(npar_psi)]

      obj_B <- compute_B_kernel(
        rho0,
        need_grad = TRUE
      )

      n_node_B_local <- n_node_B
      Wnode_B_local <- Wnode_B

    } else {

      obj_B <- NULL
      psi0 <- par[seq_len(npar_psi)]
      n_node_B_local <- 1L
      Wnode_B_local <- 1
    }

    obj_N <- compute_N_kernel(
      psi0,
      need_grad = TRUE
    )

    grad_Psi <- matrix(0, nrow = k, ncol = J)

    grad_Rho <- if (estimate_B_factor) {
      matrix(0, nrow = n_zi, ncol = L)
    } else {
      NULL
    }

    # Column ranges for the stacked dPN matrices.  Direction r
    # occupies one contiguous n_node_N-column block.
    dPN_cols <- lapply(
      seq_len(J),
      function(r) {
        ((r - 1L) * n_node_N + 1L):(r * n_node_N)
      }
    )

    for (idx in obs_blocks) {

      nb <- length(idx)
      PB_cache <- if (estimate_B_factor && L == 1L) make_PB_cache_L1(idx, obj_B, need_grad = TRUE) else NULL

      PN <- make_PN_block(
        idx,
        obj_N
      )

      ## ------------------------------------------------------
      ## PASS 1: likelihood contribution for this block.
      ## ------------------------------------------------------

      Li <- numeric(nb)
      Tprod_list <- vector("list", n_node_B_local)

      for (b in seq_len(n_node_B_local)) {

        Tprod <- matrix(
          1,
          nrow = nb,
          ncol = n_node_N
        )

        for (j in seq_len(k)) {

          PB_b <- if (zero_inflation[j]) {
            get_PB_cached(j = j, b = b, idx = idx, obj_B = obj_B, cache = PB_cache)
          } else NULL

          Tj <- make_T_block(
            j = j, b = b, idx = idx, PN = PN, obj_B = obj_B, PB_b = PB_b
          )

          Tprod <- Tprod * pmax(Tj, tiny)
        }

        Tprod_list[[b]] <- Tprod

        Li <-
          Li +
          Wnode_B_local[b] *
          as.vector(Tprod %*% Wnode_N)
      }

      Li <- pmax(Li, tiny)

      ## ------------------------------------------------------
      ## PASS 2A: count-factor gradient.
      ##
      ## Main changes relative to the original implementation:
      ##   1) all J dPN/dr are produced by ONE GEMM per margin;
      ##   2) Tj and Tprod/Tj are computed once per (j,b), not
      ##      once again for every r;
      ##   3) no full dT matrix is allocated for ZIPM margins.
      ## ------------------------------------------------------

      for (j in seq_len(k)) {

        I_y_block <- I_y_list[[j]][idx, , drop = FALSE]

        I_0_block <- if (zero_inflation[j]) {
          I_0_list[[j]][idx, , drop = FALSE]
        } else {
          NULL
        }

        # d0 x (n_node_N * J), already assembled once in
        # compute_N_kernel().
        dOmegaT_stack <- obj_N$dOmegaT_stack_list[[j]]

        # One BLAS multiplication gives derivatives for every r.
        dPN_y_all <- d0 * (I_y_block %*% dOmegaT_stack)

        dPN_0_all <- if (zero_inflation[j]) {
          d0 * (I_0_block %*% dOmegaT_stack)
        } else {
          NULL
        }

        y0 <- if (zero_inflation[j]) {
          y0_list[[j]][idx]
        } else {
          NULL
        }
        ypos <- if (zero_inflation[j]) !y0 else NULL

        g_j <- numeric(J)

        for (b in seq_len(n_node_B_local)) {

          PB_b <- if (zero_inflation[j]) {
            get_PB_cached(j = j, b = b, idx = idx, obj_B = obj_B, cache = PB_cache)
          } else {
            NULL
          }

          Tj <- make_T_block(
            j = j,
            b = b,
            idx = idx,
            PN = PN,
            obj_B = obj_B,
            PB_b = PB_b
          )

          T_except_j <-
            Tprod_list[[b]] / pmax(Tj, tiny)

          for (r in seq_len(J)) {

            cols_r <- dPN_cols[[r]]
            dPN_y_r <- dPN_y_all[, cols_r, drop = FALSE]

            if (!zero_inflation[j]) {

              dLi <- as.vector(
                (T_except_j * dPN_y_r) %*% Wnode_N
              )

            } else {

              dPN_0_r <- dPN_0_all[, cols_r, drop = FALSE]
              dLi <- numeric(nb)

              # PB_b is constant over the count-factor GH nodes,
              # so apply it after the node integration instead of
              # constructing an n x n_node_N dT matrix.
              if (any(y0)) {
                dLi[y0] <-
                  PB_b[y0] *
                  as.vector(
                    (
                      T_except_j[y0, , drop = FALSE] *
                        dPN_0_r[y0, , drop = FALSE]
                    ) %*% Wnode_N
                  )
              }

              if (any(ypos)) {
                dLi[ypos] <-
                  PB_b[ypos] *
                  as.vector(
                    (
                      T_except_j[ypos, , drop = FALSE] *
                        dPN_y_r[ypos, , drop = FALSE]
                    ) %*% Wnode_N
                  )
              }
            }

            g_j[r] <-
              g_j[r] -
              sum(
                Wnode_B_local[b] *
                  dLi / Li
              )
          }
        }

        grad_Psi[j, ] <- grad_Psi[j, ] + g_j

        rm(dPN_y_all, dPN_0_all)
      }

      ## ------------------------------------------------------
      ## PASS 2B: zero-factor gradient.
      ##
      ## Tj and Tprod/Tj are now computed once per (margin,node)
      ## and reused for every zero-factor direction l.  Moreover,
      ## the expensive n x n_node_N part that is independent of l
      ## is integrated once; each l then only scales length-n
      ## vectors by dPB.
      ## ------------------------------------------------------

      if (estimate_B_factor) {

        for (s in seq_len(n_zi)) {

          j <- zi_index[s]
          y0 <- y0_list[[j]][idx]
          ypos <- !y0
          g_s <- numeric(L)

          for (b in seq_len(n_node_B_local)) {

            PB_b <- get_PB_cached(j = j, b = b, idx = idx, obj_B = obj_B, cache = PB_cache)

            Tj <- make_T_block(
              j = j,
              b = b,
              idx = idx,
              PN = PN,
              obj_B = obj_B,
              PB_b = PB_b
            )

            T_except_j <-
              Tprod_list[[b]] / pmax(Tj, tiny)

            # Integrate the part independent of l only once.
            base_zero <- if (any(y0)) {
              as.vector(
                (
                  T_except_j[y0, , drop = FALSE] *
                    (PN$PN_0[[j]][y0, , drop = FALSE] - 1)
                ) %*% Wnode_N
              )
            } else {
              numeric(0)
            }

            base_pos <- if (any(ypos)) {
              as.vector(
                (
                  T_except_j[ypos, , drop = FALSE] *
                    PN$PN_y[[j]][ypos, , drop = FALSE]
                ) %*% Wnode_N
              )
            } else {
              numeric(0)
            }

            for (l in seq_len(L)) {

              dPB_b <- get_dPB_cached(
                j = j, b = b, idx = idx, obj_B = obj_B, cache = PB_cache
              )

              dLi <- numeric(nb)

              if (any(y0)) {
                dLi[y0] <- dPB_b[y0] * base_zero
              }

              if (any(ypos)) {
                dLi[ypos] <- dPB_b[ypos] * base_pos
              }

              g_s[l] <-
                g_s[l] -
                sum(
                  Wnode_B_local[b] *
                    dLi / Li
                )
            }
          }

          grad_Rho[s, ] <- grad_Rho[s, ] + g_s
        }
      }

      rm(PN, Tprod_list, Li)
    }

    grad_psi0 <- chain_grad_loading_to_raw(
      grad_Load = grad_Psi,
      Araw = obj_N$A,
      map = obj_N$map,
      npar = npar_psi,
      q = J
    )

    if (estimate_B_factor) {

      grad_rho0 <- chain_grad_loading_to_raw(
        grad_Load = grad_Rho,
        Araw = obj_B$A,
        map = obj_B$map,
        npar = npar_rho,
        q = L
      )

      c(grad_rho0, grad_psi0)

    } else {

      grad_psi0
    }
  }

  ## =========================================================
  ## 10. Starting values
  ## =========================================================

  if (is.null(psi0_start)) {

    if (latent.dependence == "gaussian" && latent.loglike == "hermite") {

      psi0_start <- make_factor_start_nonzero(
        k = k,
        q = J,
        loading0 = 0.05,
        offdiag0 = 0.02
      )

    } else {

      psi0_start <- make_factor_start(
        k = k,
        q = J,
        loading0 = 0.05
      )
    }
  }

  if (length(psi0_start) != npar_psi) {
    stop("psi0_start must have length ", npar_psi, ".")
  }

  if (estimate_B_factor) {

    if (is.null(rho0_start)) {
      ## Joint-Hermite parameterization writes rho_s = t_s v_s and
      ## therefore requires every zero-loading row to have t_s > 0.
      ## The generic triangular initializer leaves all rows s > L
      ## exactly zero, which makes v_s = rho_s / ||rho_s|| undefined.
      ## Use the same non-zero initializer as for the count factor.
      ## Start close to the nested L=0 model.  A clearly non-zero radius
      ## is retained only because the polar joint-Hermite parameterization
      ## needs a defined direction at the starting point.
      rho0_start <- make_factor_start_nonzero(
        k = n_zi,
        q = L,
        loading0 = 0.005,
        offdiag0 = 0.002
      )
    }

    if (length(rho0_start) != npar_rho) {
      stop("rho0_start must have length ", npar_rho, ".")
    }

  } else {

    rho0_start <- numeric(0)
  }

  par_start <- if (estimate_B_factor) {
    c(rho0_start, psi0_start)
  } else {
    psi0_start
  }

  ## =========================================================
  ## 10b. Hermite L=0 polar optimizer helpers
  ##
  ## Psi[j,] = r_j v_j, ||v_j||=1.  For fixed r, the expensive
  ## one-dimensional Hermite count/grid quantities depend only on r and
  ## are cached across all V evaluations.  The lower bound and E(R)
  ## formula are unchanged.
  ## =========================================================

  .polar_active_dim <- function(j) min(j, J)

  .polar_from_Psi <- function(Psi) {
    r <- sqrt(rowSums(Psi^2))
    if (any(!is.finite(r)) || any(r <= 1e-10) || any(r >= 1)) {
      stop("Invalid Psi for Hermite polar parameterization.")
    }
    V <- Psi
    for (j in seq_len(k)) {
      p <- .polar_active_dim(j)
      V[j, ] <- 0
      V[j, seq_len(p)] <- Psi[j, seq_len(p)] / r[j]
    }
    list(r = r, V = V)
  }

  .polar_to_Psi <- function(r, V) {
    Psi <- V * r
    for (j in seq_len(k)) {
      p <- .polar_active_dim(j)
      if (p < J) Psi[j, (p + 1L):J] <- 0
    }
    Psi
  }

  # Fix the usual Gaussian-factor column-sign ambiguity so that the
  # identifying diagonal remains positive when Psi is converted back to
  # the original lower-triangular raw parameterization.
  .polar_fix_column_signs <- function(V) {
    for (a in seq_len(min(J, k))) {
      if (is.finite(V[a, a]) && V[a, a] < 0) V[, a] <- -V[, a]
    }
    V
  }

  .polar_tangent_bases <- function(V) {
    out <- vector("list", k)
    for (j in seq_len(k)) {
      p <- .polar_active_dim(j)
      if (p <= 1L) {
        out[[j]] <- matrix(0, nrow = p, ncol = 0L)
        next
      }
      v <- V[j, seq_len(p)]
      # null(t(v)) gives an orthonormal basis of the tangent space.
      out[[j]] <- qr.Q(qr(cbind(v, diag(p))), complete = TRUE)[, -1L, drop = FALSE]
    }
    out
  }

  .polar_pack_delta0 <- function() {
    numeric(sum(pmax(vapply(seq_len(k), .polar_active_dim, integer(1)) - 1L, 0L)))
  }

  .polar_delta_to_V <- function(delta, V0, Tb) {
    V <- matrix(0, nrow = k, ncol = J)
    pos <- 1L
    for (j in seq_len(k)) {
      p <- .polar_active_dim(j)
      v0 <- V0[j, seq_len(p)]
      if (p == 1L) {
        V[j, 1L] <- v0[1L]
        next
      }
      m <- p - 1L
      dd <- delta[pos:(pos + m - 1L)]
      pos <- pos + m
      y <- v0 + as.vector(Tb[[j]] %*% dd)
      sc <- sqrt(1 + sum(dd^2))
      V[j, seq_len(p)] <- y / sc
    }
    V
  }

  .polar_grad_delta <- function(delta, V0, Tb, r, grad_Psi) {
    g <- numeric(length(delta))
    pos <- 1L
    for (j in seq_len(k)) {
      p <- .polar_active_dim(j)
      if (p == 1L) next
      m <- p - 1L
      dd <- delta[pos:(pos + m - 1L)]
      sc <- sqrt(1 + sum(dd^2))
      v <- (V0[j, seq_len(p)] + as.vector(Tb[[j]] %*% dd)) / sc
      # dv/dd = T/sc - v dd'/sc^2
      Jac <- Tb[[j]] / sc - tcrossprod(v, dd) / (sc^2)
      g[pos:(pos + m - 1L)] <-
        r[j] * as.vector(crossprod(Jac, grad_Psi[j, seq_len(p)]))
      pos <- pos + m
    }
    g
  }

  .hermite_L0_cache_for_r <- function(r, Vref, need_grad = TRUE) {
    Psi_ref <- .polar_to_Psi(r, Vref)
    Hker <- compute_H_kernel(Psi_ref)
    block_cache <- lapply(
      obs_blocks,
      function(idx) make_hermite_count_cache_block(idx, Hker, need_grad = need_grad)
    )
    list(Hker = Hker, block_cache = block_cache, need_grad = need_grad)
  }

  .hermite_L0_eval_Psi <- function(Psi, need_grad = FALSE, fixed_cache = NULL) {
    if (estimate_B_factor || L != 0L) stop("Internal polar Hermite helper requires L=0.")

    if (is.null(fixed_cache)) {
      Hker <- compute_H_kernel(Psi)
      block_cache <- NULL
    } else {
      Hker <- fixed_cache$Hker
      # Only the row norms determine the interval/count cache.  Replace
      # the orientation used by the Gaussian state/finish kernels.
      Hker$Psi <- Psi
      block_cache <- fixed_cache$block_cache
    }

    nll <- 0
    grad_Psi <- if (need_grad) matrix(0, nrow = k, ncol = J) else NULL

    for (bb in seq_along(obs_blocks)) {
      idx <- obs_blocks[[bb]]
      cc <- if (is.null(block_cache)) {
        make_hermite_count_cache_block(idx, Hker, need_grad = need_grad)
      } else {
        block_cache[[bb]]
      }

      hh <- f_int_hermite_fast_batch(
        idx = idx,
        b = 1L,
        obj_B = NULL,
        Hker = Hker,
        need_grad = need_grad,
        obs_weight = if (need_grad) rep(1, length(idx)) else NULL,
        count_cache = cc
      )

      if (any(!is.finite(hh$log_value))) {
        if (need_grad) return(list(value = big_penalty, grad_Psi = matrix(0, k, J)))
        return(list(value = big_penalty, grad_Psi = NULL))
      }
      nll <- nll - sum(hh$log_value)
      if (need_grad) grad_Psi <- grad_Psi - hh$grad_Psi
    }

    list(value = if (is.finite(nll)) nll else big_penalty, grad_Psi = grad_Psi)
  }

  .hermite_polar_optimize_L0 <- function(psi_raw_start) {
    pobj0 <- make_factor_loadings(psi_raw_start, k = k, q = J, return_all = TRUE)
    pv <- .polar_from_Psi(pobj0$Psi)
    r <- pv$r
    V <- .polar_fix_column_signs(pv$V)

    # Keep radii safely inside the admissible interval.
    r <- pmin(r, hermite.polar.r.max * (1 - 1e-8))

    cur <- .hermite_L0_eval_Psi(.polar_to_Psi(r, V), need_grad = FALSE)$value
    total_iter <- 0L
    total_eval <- 1L
    conv <- 1L

    for (cy in seq_len(max(1L, hermite.polar.cycles))) {
      before <- cur

      ## ----------------------------------------------------------
      ## V | r : cache all r-dependent Hermite count/grid work.
      ## ----------------------------------------------------------
      if (hermite.polar.v.iter > 0L && J > 1L) {
        Tb <- .polar_tangent_bases(V)
        delta0 <- .polar_pack_delta0()
        rc <- .hermite_L0_cache_for_r(r, V, need_grad = TRUE)

        fV <- function(delta) {
          VV <- .polar_delta_to_V(delta, V, Tb)
          zz <- .hermite_L0_eval_Psi(.polar_to_Psi(r, VV), need_grad = FALSE, fixed_cache = rc)
          total_eval <<- total_eval + 1L
          zz$value
        }
        gV <- function(delta) {
          VV <- .polar_delta_to_V(delta, V, Tb)
          zz <- .hermite_L0_eval_Psi(.polar_to_Psi(r, VV), need_grad = TRUE, fixed_cache = rc)
          total_eval <<- total_eval + 1L
          .polar_grad_delta(delta, V, Tb, r, zz$grad_Psi)
        }

        oV <- nlminb(
          start = delta0,
          objective = fV,
          gradient = gV,
          control = list(
            trace = 0,
            iter.max = hermite.polar.v.iter,
            eval.max = max(20L, 5L * hermite.polar.v.iter),
            rel.tol = hermite.polar.tol
          )
        )
        Vcand <- .polar_delta_to_V(oV$par, V, Tb)
        vc <- .hermite_L0_eval_Psi(.polar_to_Psi(r, Vcand), need_grad = FALSE, fixed_cache = rc)$value
        total_eval <- total_eval + 1L
        if (is.finite(vc) && vc <= cur + 1e-8 * (1 + abs(cur))) {
          V <- .polar_fix_column_signs(Vcand)
          cur <- vc
        }
        total_iter <- total_iter + ifelse(is.null(oV$iterations), 0L, oV$iterations)
      }

      afterV <- cur

      ## ----------------------------------------------------------
      ## r | V : optimize all k radii together.  This is the expensive
      ## block because changing r invalidates the Hermite count/grid cache.
      ## ----------------------------------------------------------
      if (hermite.polar.r.iter > 0L) {
        rr <- pmin(pmax(r / hermite.polar.r.max, 1e-8), 1 - 1e-8)
        eta0 <- qlogis(rr)

        fR <- function(eta) {
          rt <- hermite.polar.r.max * plogis(eta)
          zz <- .hermite_L0_eval_Psi(.polar_to_Psi(rt, V), need_grad = FALSE)
          total_eval <<- total_eval + 1L
          zz$value
        }
        gR <- function(eta) {
          pp <- plogis(eta)
          rt <- hermite.polar.r.max * pp
          zz <- .hermite_L0_eval_Psi(.polar_to_Psi(rt, V), need_grad = TRUE)
          total_eval <<- total_eval + 1L
          gr <- rowSums(zz$grad_Psi * V)
          gr * (hermite.polar.r.max * pp * (1 - pp))
        }

        oR <- nlminb(
          start = eta0,
          objective = fR,
          gradient = gR,
          control = list(
            trace = 0,
            iter.max = hermite.polar.r.iter,
            eval.max = max(20L, 5L * hermite.polar.r.iter),
            rel.tol = hermite.polar.tol
          )
        )
        rcand <- hermite.polar.r.max * plogis(oR$par)
        rv <- .hermite_L0_eval_Psi(.polar_to_Psi(rcand, V), need_grad = FALSE)$value
        total_eval <- total_eval + 1L
        if (is.finite(rv) && rv <= cur + 1e-8 * (1 + abs(cur))) {
          r <- rcand
          cur <- rv
        }
        total_iter <- total_iter + ifelse(is.null(oR$iterations), 0L, oR$iterations)
      }

      improve <- before - cur
      if (!isTRUE(.gh_progressive_internal)) {
        .internal_msg(
          "[dependence][Hermite polar] cycle=", cy,
          " | logLik=", signif(-cur, 10),
          " | after V=", signif(-afterV, 10),
          " | improvement=", signif(improve, 6),
          " | r range=[", signif(min(r), 4), ", ", signif(max(r), 4), "]"
        )
      }

      if (improve <= hermite.polar.tol * (1 + abs(before))) {
        conv <- 0L
        break
      }
    }

    V <- .polar_fix_column_signs(V)
    Psi_final <- .polar_to_Psi(r, V)
    raw_final <- factor_start_from_Psi(Psi_final, max_row_norm = min(0.999999, hermite.polar.r.max))

    list(
      par = raw_final,
      objective = cur,
      convergence = conv,
      message = if (conv == 0L) "Hermite polar alternating converged" else "Hermite polar alternating iteration limit",
      iterations = total_iter,
      evaluations = c(function_eval = total_eval),
      polar_r = r,
      polar_V = V
    )
  }


  ## =========================================================
  ## Hermite polar optimizer for L = 1
  ##
  ## Alternating blocks:
  ##   V | r,rho  ->  rho | V,r  ->  r | V,rho
  ##
  ## Count Hermite interval/grid quantities depend on row norms r only.
  ## Hence V and rho updates reuse the same count cache.  Only the r block
  ## rebuilds it.  The L=0 optimizer above is intentionally left unchanged.
  ## =========================================================

  .hermite_L1_cache_for_r <- function(r, Vref, need_grad = TRUE) {
    Psi_ref <- .polar_to_Psi(r, Vref)
    Hker <- compute_H_kernel(Psi_ref)
    block_cache <- lapply(
      obs_blocks,
      function(idx) make_hermite_count_cache_block(idx, Hker, need_grad = need_grad)
    )
    list(Hker = Hker, block_cache = block_cache, need_grad = need_grad)
  }

  .hermite_L1_eval <- function(Psi, rho_vec, need_grad = FALSE, fixed_cache = NULL) {
    if (!estimate_B_factor || L != 1L) stop("Internal L1 polar Hermite helper requires L=1 with an estimated zero factor.")
    Rho <- matrix(rho_vec, ncol = 1L)
    if (!.joint_hermite_loadings_valid(Psi, Rho)) {
      return(list(value = big_penalty,
                  grad_Psi = if (need_grad) matrix(0, k, J) else NULL,
                  grad_Rho = if (need_grad) rep(0, n_zi) else NULL))
    }

    if (is.null(fixed_cache)) {
      Hker <- compute_H_kernel(Psi)
      block_cache <- NULL
    } else {
      Hker <- fixed_cache$Hker
      Hker$Psi <- Psi
      block_cache <- fixed_cache$block_cache
    }

    pobj <- list(Psi = Psi)
    obj_B <- list(Rho = Rho)
    nll <- 0
    grad_Psi <- if (need_grad) matrix(0, k, J) else NULL
    grad_Rho <- if (need_grad) rep(0, n_zi) else NULL

    for (bb in seq_along(obs_blocks)) {
      idx <- obs_blocks[[bb]]
      z <- try(.joint_hermite_block(
        idx = idx, pobj = pobj, obj_B = obj_B, Hker = Hker,
        need_grad = need_grad,
        count_cache = if (is.null(block_cache)) NULL else block_cache[[bb]]
      ), silent = TRUE)
      if (inherits(z, "try-error")) {
        return(list(value = big_penalty,
                    grad_Psi = if (need_grad) matrix(0, k, J) else NULL,
                    grad_Rho = if (need_grad) rep(0, n_zi) else NULL))
      }
      lv <- z$log_value
      ok <- is.finite(lv)
      if (!all(ok)) {
        if (!any(ok)) return(list(value = big_penalty,
                                  grad_Psi = if (need_grad) matrix(0, k, J) else NULL,
                                  grad_Rho = if (need_grad) rep(0, n_zi) else NULL))
        fill <- stats::median(lv[ok]); if (!is.finite(fill)) fill <- mean(lv[ok])
        lv[!ok] <- fill
      }
      nll <- nll - sum(lv)
      if (need_grad) {
        z$grad_Psi[!is.finite(z$grad_Psi)] <- 0
        z$grad_Rho[!is.finite(z$grad_Rho)] <- 0
        grad_Psi <- grad_Psi - z$grad_Psi
        grad_Rho <- grad_Rho - as.vector(z$grad_Rho[, 1L])
      }
    }
    list(value = if (is.finite(nll)) nll else big_penalty,
         grad_Psi = grad_Psi, grad_Rho = grad_Rho)
  }

  .hermite_polar_optimize_L1 <- function(rho_raw_start, psi_raw_start) {
    pobj0 <- make_factor_loadings(psi_raw_start, k = k, q = J, return_all = TRUE)
    pv <- .polar_from_Psi(pobj0$Psi)
    r <- pmin(pv$r, hermite.polar.r.max * (1 - 1e-8))
    V <- .polar_fix_column_signs(pv$V)

    rho_obj0 <- make_factor_loadings(rho_raw_start, k = n_zi, q = 1L, return_all = TRUE)
    rho <- as.vector(rho_obj0$Psi[, 1L])
    rho_max <- hermite.polar.r.max
    rho <- pmin(pmax(rho, -rho_max * (1 - 1e-8)), rho_max * (1 - 1e-8))

    cur <- .hermite_L1_eval(.polar_to_Psi(r, V), rho, need_grad = FALSE)$value
    total_iter <- 0L; total_eval <- 1L; conv <- 1L

    for (cy in seq_len(max(1L, hermite.polar.cycles))) {
      before <- cur
      rc <- .hermite_L1_cache_for_r(r, V, need_grad = TRUE)

      ## V | r,rho -- cheap block: reuse all r-dependent count work.
      if (hermite.polar.v.iter > 0L && J > 1L) {
        Tb <- .polar_tangent_bases(V); delta0 <- .polar_pack_delta0()
        fV <- function(delta) {
          VV <- .polar_delta_to_V(delta, V, Tb)
          zz <- .hermite_L1_eval(.polar_to_Psi(r, VV), rho, FALSE, rc)
          total_eval <<- total_eval + 1L; zz$value
        }
        gV <- function(delta) {
          VV <- .polar_delta_to_V(delta, V, Tb)
          zz <- .hermite_L1_eval(.polar_to_Psi(r, VV), rho, TRUE, rc)
          total_eval <<- total_eval + 1L
          .polar_grad_delta(delta, V, Tb, r, zz$grad_Psi)
        }
        oV <- nlminb(delta0, fV, gradient=gV,
                     control=list(trace=0, iter.max=hermite.polar.v.iter,
                                  eval.max=max(20L,5L*hermite.polar.v.iter), rel.tol=hermite.polar.tol))
        Vcand <- .polar_delta_to_V(oV$par, V, Tb)
        vc <- .hermite_L1_eval(.polar_to_Psi(r,Vcand), rho, FALSE, rc)$value
        total_eval <- total_eval + 1L
        if (is.finite(vc) && vc <= cur + 1e-8*(1+abs(cur))) { V <- .polar_fix_column_signs(Vcand); cur <- vc }
        total_iter <- total_iter + ifelse(is.null(oV$iterations),0L,oV$iterations)
      }
      afterV <- cur

      ## rho | V,r -- L=1 has no continuous direction parameter; optimize
      ## signed zero-factor loadings through rho=rho_max*tanh(alpha).
      alpha0 <- atanh(pmin(pmax(rho/rho_max, -1+1e-8), 1-1e-8))
      fB <- function(alpha) {
        rh <- rho_max * tanh(alpha)
        zz <- .hermite_L1_eval(.polar_to_Psi(r,V), rh, FALSE, rc)
        total_eval <<- total_eval + 1L; zz$value
      }
      gB <- function(alpha) {
        th <- tanh(alpha); rh <- rho_max * th
        zz <- .hermite_L1_eval(.polar_to_Psi(r,V), rh, TRUE, rc)
        total_eval <<- total_eval + 1L
        zz$grad_Rho * rho_max * (1-th^2)
      }
      oB <- nlminb(alpha0, fB, gradient=gB,
                   control=list(trace=0, iter.max=hermite.polar.r.iter,
                                eval.max=max(20L,5L*hermite.polar.r.iter), rel.tol=hermite.polar.tol))
      rhocand <- rho_max*tanh(oB$par)
      bv <- .hermite_L1_eval(.polar_to_Psi(r,V), rhocand, FALSE, rc)$value
      total_eval <- total_eval + 1L
      if (is.finite(bv) && bv <= cur + 1e-8*(1+abs(cur))) { rho <- rhocand; cur <- bv }
      total_iter <- total_iter + ifelse(is.null(oB$iterations),0L,oB$iterations)
      afterRho <- cur

      ## r | V,rho -- expensive block; changing r invalidates count cache.
      if (hermite.polar.r.iter > 0L) {
        eta0 <- qlogis(pmin(pmax(r/hermite.polar.r.max,1e-8),1-1e-8))
        fR <- function(eta) {
          rt <- hermite.polar.r.max*plogis(eta)
          zz <- .hermite_L1_eval(.polar_to_Psi(rt,V), rho, FALSE)
          total_eval <<- total_eval + 1L; zz$value
        }
        gR <- function(eta) {
          pp <- plogis(eta); rt <- hermite.polar.r.max*pp
          zz <- .hermite_L1_eval(.polar_to_Psi(rt,V), rho, TRUE)
          total_eval <<- total_eval + 1L
          gr <- rowSums(zz$grad_Psi*V)
          gr*(hermite.polar.r.max*pp*(1-pp))
        }
        oR <- nlminb(eta0, fR, gradient=gR,
                     control=list(trace=0, iter.max=hermite.polar.r.iter,
                                  eval.max=max(20L,5L*hermite.polar.r.iter), rel.tol=hermite.polar.tol))
        rcand <- hermite.polar.r.max*plogis(oR$par)
        rv <- .hermite_L1_eval(.polar_to_Psi(rcand,V), rho, FALSE)$value
        total_eval <- total_eval + 1L
        if (is.finite(rv) && rv <= cur + 1e-8*(1+abs(cur))) { r <- rcand; cur <- rv }
        total_iter <- total_iter + ifelse(is.null(oR$iterations),0L,oR$iterations)
      }

      improve <- before-cur
      .internal_msg("[dependence][Hermite polar L1] cycle=",cy,
                    " | logLik=",signif(-cur,10),
                    " | after V=",signif(-afterV,10),
                    " | after rho=",signif(-afterRho,10),
                    " | improvement=",signif(improve,6),
                    " | r=[",signif(min(r),4),", ",signif(max(r),4),"]",
                    " | rho=[",signif(min(rho),4),", ",signif(max(rho),4),"]")
      if (improve <= hermite.polar.tol*(1+abs(before))) { conv <- 0L; break }
    }

    Psi_final <- .polar_to_Psi(r, .polar_fix_column_signs(V))
    psi_raw_final <- factor_start_from_Psi(Psi_final, max_row_norm=min(0.999999,hermite.polar.r.max))
    rho_raw_final <- factor_start_from_Psi(matrix(rho,ncol=1L), max_row_norm=min(0.999999,rho_max))
    list(par=c(rho_raw_final,psi_raw_final), objective=cur, convergence=conv,
         message=if(conv==0L) "Hermite polar L1 alternating converged" else "Hermite polar L1 alternating iteration limit",
         iterations=total_iter, evaluations=c(function_eval=total_eval),
         polar_r=r, polar_V=V, polar_rho=rho)
  }

  ## =========================================================
  ## 11. Optimization
  ##
  ## Keep nlminb() unchanged.  A combined objective+gradient
  ## evaluator is intentionally NOT forced here: nlminb may perform
  ## objective-only line-search evaluations, so always calculating
  ## the gradient inside objective() can be slower overall.
  ## =========================================================

  if (latent.dependence == "gaussian" && latent.loglike == "hermite") {

    psi_start_check <- if (estimate_B_factor) {
      par_start[npar_rho + seq_len(npar_psi)]
    } else {
      par_start[seq_len(npar_psi)]
    }

    psi_start_obj <- make_factor_loadings(
      par = psi_start_check,
      k = k,
      q = J,
      return_all = TRUE
    )

    start_row_norm <- sqrt(rowSums(psi_start_obj$Psi^2))

    if (estimate_B_factor && use_joint_hermite_allL) {
      rho_start_check <- par_start[seq_len(npar_rho)]
      rho_start_obj <- make_factor_loadings(
        par = rho_start_check,
        k = n_zi,
        q = L,
        return_all = TRUE
      )
      rho_start_norm <- sqrt(rowSums(rho_start_obj$Psi^2))

      if (any(!is.finite(rho_start_norm)) ||
          any(rho_start_norm >= 1 - 1e-10)) {
        stop(
          "Invalid initial zero-factor loading row in joint Hermite backend. ",
          "Range of ||rho_j|| = [",
          signif(min(rho_start_norm), 6), ", ",
          signif(max(rho_start_norm), 6), "]."
        )
      }
    }

    if (any(start_row_norm <= 1e-8)) {
      stop(
        "Hermite-polynomial analytic gradient requires every starting loading row ",
        "to have positive norm. Supply a non-zero `psi0_start`."
      )
    }
  }


  if (!isTRUE(.gh_progressive_internal) && use_joint_hermite_allL) {
    .internal_msg(
      "[dependence] Hermite backend: joint Gaussian core dimension J+L=", J + L,
      " | residual dimension <= 2"
    )
  }
  if (!isTRUE(.gh_progressive_internal) && !isTRUE(.gh_node_progressive_internal)) {
    .progress_msg("Initial parameters are ready. Starting optimization...")
  }

  t0 <- proc.time()[3L]

  use_analytic_gradient_now <- isTRUE(use.gradient)

  ## Save a known-good fallback before the unrestricted final optimization.
  ## Progressive warm fitting has already brought par_start into the local
  ## region where the Hermite approximation is useful.
  start_objective_final <- try(negloglik_copula(par_start), silent = TRUE)
  if (inherits(start_objective_final, "try-error") ||
      !is.finite(start_objective_final)) {
    start_objective_final <- big_penalty
  }

  use_polar_hermite_now <-
    identical(hermite.optimizer, "polar") &&
    identical(latent.dependence, "gaussian") &&
    identical(latent.loglike, "hermite") &&
    ((!estimate_B_factor && L == 0L) || (estimate_B_factor && L == 1L && use_joint_hermite_allL)) &&
    !isTRUE(.gh_progressive_internal)

  if (isTRUE(use_polar_hermite_now)) {

    .internal_msg(
      if (estimate_B_factor && L == 1L)
        "[dependence] Hermite optimizer L1: alternating V|r,rho -> rho|V,r -> r|V,rho"
      else
        "[dependence] Hermite optimizer: alternating V|r -> r|V",
      " | cycles=", hermite.polar.cycles,
      " | V iter=", hermite.polar.v.iter,
      " | r iter=", hermite.polar.r.iter
    )

    opt <- if (estimate_B_factor && L == 1L) {
      .hermite_polar_optimize_L1(
        rho_raw_start = par_start[seq_len(npar_rho)],
        psi_raw_start = par_start[npar_rho + seq_len(npar_psi)]
      )
    } else {
      .hermite_polar_optimize_L0(par_start[seq_len(npar_psi)])
    }

    # Optional short joint-Psi refinement.  This keeps the speed benefit of
    # the polar phase but lets the original parameterization polish the final
    # stationary point.
    if (!estimate_B_factor && hermite.polar.joint.refine.iter > 0L && is.finite(opt$objective)) {
      oj <- if (use_analytic_gradient_now) {
        nlminb(
          start = opt$par, objective = negloglik_copula,
          gradient = grad_negloglik_copula,
          control = list(
            trace = 0,
            eval.max = max(20L, 5L * hermite.polar.joint.refine.iter),
            iter.max = hermite.polar.joint.refine.iter,
            rel.tol = rel.tol
          )
        )
      } else {
        nlminb(
          start = opt$par, objective = negloglik_copula,
          control = list(
            trace = 0,
            eval.max = max(20L, 5L * hermite.polar.joint.refine.iter),
            iter.max = hermite.polar.joint.refine.iter,
            rel.tol = rel.tol
          )
        )
      }
      if (is.finite(oj$objective) && oj$objective < opt$objective) {
        oj$message <- paste0("polar + joint refinement: ", oj$message)
        opt <- oj
      }
    }

  } else if (use_analytic_gradient_now) {

    opt <- nlminb(
      start = par_start,
      objective = negloglik_copula,
      gradient = grad_negloglik_copula,
      control = list(
        trace = trace,
        eval.max = eval.max,
        iter.max = iter.max, rel.tol = rel.tol
      )
    )

  } else {

    opt <- nlminb(
      start = par_start,
      objective = negloglik_copula,
      control = list(
        trace = trace,
        eval.max = eval.max,
        iter.max = iter.max,rel.tol = rel.tol
      )
    )
  }

  elapsed <- proc.time()[3L] - t0

  ## v16 final-solution sanity check -------------------------------------
  ## The Hermite objective is a surrogate, so its value is not constrained to
  ## be a genuine log probability at every trial point.  Reject only numerical
  ## failure / explicit penalty.  Exact model validity is checked by final RQMC.
  opt_invalid <-
    !is.finite(opt$objective) ||
    opt$objective >= 0.999 * big_penalty

  if (isTRUE(opt_invalid) &&
      is.finite(start_objective_final) &&
      start_objective_final < 0.999 * big_penalty) {
    if (!isTRUE(.gh_progressive_internal)) {
      .internal_msg(
        "[dependence] WARNING: invalid Hermite final optimum detected; ",
        "reverting to progressive warm start | bad approx logLik=",
        signif(-opt$objective, 10),
        " | warm approx logLik=", signif(-start_objective_final, 10)
      )
    }
    opt$par <- par_start
    opt$objective <- start_objective_final
    opt$convergence <- 90L
    opt$message <- paste0(
      "Hermite surrogate left valid probability domain; ",
      "reverted to progressive warm start. Original nlminb message: ",
      if (is.null(opt$message)) "" else opt$message
    )
  }

  if (!isTRUE(.gh_progressive_internal) && !isTRUE(.gh_node_progressive_internal)) {
    .progress_msg("Optimization completed.")
    if (identical(latent.dependence, "gaussian") && identical(latent.loglike, "hermite")) {
      .progress_msg("  Hermite lower-bound log-likelihood: ", signif(-opt$objective, 10))
    } else {
      .progress_msg("  Log-likelihood: ", signif(-opt$objective, 10))
    }
    .progress_msg("  Convergence code: ", opt$convergence)
    .progress_msg("  Optimization time: ", signif(elapsed, 5), " sec")
  }

  ## =========================================================
  ## 12. Recover loadings
  ## =========================================================

  if (estimate_B_factor) {

    rho0_hat <- opt$par[seq_len(npar_rho)]
    psi0_hat <- opt$par[npar_rho + seq_len(npar_psi)]

    rho_hat_obj <- make_factor_loadings(
  	par = rho0_hat,
  	k = n_zi,
  	q = L,
  	return_all = TRUE
    )
    R_B_hat <- rho_hat_obj$Psi %*% t(rho_hat_obj$Psi)
    diag(R_B_hat) <- 1

  } else {

    rho0_hat <- numeric(0)
    rho_hat_obj <- NULL
    R_B_hat <- NULL
    psi0_hat <- opt$par[seq_len(npar_psi)]
  }

  psi_hat_obj <- make_factor_loadings(
    par = psi0_hat,
    k = k,
    q = J,
    return_all = TRUE
  )

  
  R_N_hat <- psi_hat_obj$Psi %*% t(psi_hat_obj$Psi)
  diag(R_N_hat) <- 1

  ## v19 RQMC candidate safeguard ---------------------------------------
  ## The Hermite surrogate and the original likelihood need not rank
  ## parameter vectors identically.  For the joint L>=1 Gaussian-Hermite
  ## path, cheaply compare:
  ##   (a) the final Hermite optimum, and
  ##   (b) the progressive warm start
  ## using the SAME Sobol seed/shifts and a small Q.  If the warm start is
  ## clearly better under the original likelihood, retain it for the final
  ## high-accuracy RQMC evaluation.
  ##
  ## This is deliberately a candidate screen, not a noisy optimizer.
  if (use_joint_hermite_allL && !isTRUE(.skip_final_qmc)) {
    q_screen <- as.integer(min(qmc.Q, max(256L, min(512L, qmc.Q))))
    rep_screen <- as.integer(min(qmc.n_rep, 2L))

    get_joint_parts <- function(pp) {
      list(
        rho0 = pp[seq_len(npar_rho)],
        psi0 = pp[npar_rho + seq_len(npar_psi)]
      )
    }

    cand_opt <- get_joint_parts(opt$par)
    cand_warm <- get_joint_parts(par_start)

    screen_opt <- try(
      qmc_loglik_joint_gaussian_allL(
        psi0 = cand_opt$psi0, rho0 = cand_opt$rho0,
        Q_eval = q_screen, n_rep_eval = rep_screen
      ),
      silent = TRUE
    )
    screen_warm <- try(
      qmc_loglik_joint_gaussian_allL(
        psi0 = cand_warm$psi0, rho0 = cand_warm$rho0,
        Q_eval = q_screen, n_rep_eval = rep_screen
      ),
      silent = TRUE
    )

    if (!inherits(screen_opt, "try-error") &&
        !inherits(screen_warm, "try-error") &&
        is.finite(screen_opt$log_lik_value) &&
        is.finite(screen_warm$log_lik_value)) {

      .internal_msg(
        "[dependence] v19 RQMC screen | Hermite optimum=",
        signif(screen_opt$log_lik_value, 10),
        " | warm start=", signif(screen_warm$log_lik_value, 10),
        " | Q=", q_screen, " x ", rep_screen
      )

      ## Require a modest margin to avoid switching on screening noise.
      screen_margin <- max(
        2.0,
        3.0 * max(
          ifelse(is.finite(screen_opt$log_lik_se), screen_opt$log_lik_se, 0),
          ifelse(is.finite(screen_warm$log_lik_se), screen_warm$log_lik_se, 0)
        )
      )

      if (screen_warm$log_lik_value >
          screen_opt$log_lik_value + screen_margin) {
        .internal_msg(
          "[dependence] v19 RQMC safeguard: progressive warm start is better ",
          "for the original likelihood; using warm candidate for final RQMC."
        )
        opt$par <- par_start
        opt$objective <- start_objective_final

        rho0_hat <- opt$par[seq_len(npar_rho)]
        psi0_hat <- opt$par[npar_rho + seq_len(npar_psi)]
        rho_hat_obj <- make_factor_loadings(
          par = rho0_hat, k = n_zi, q = L, return_all = TRUE
        )
        R_B_hat <- rho_hat_obj$Psi %*% t(rho_hat_obj$Psi)
        diag(R_B_hat) <- 1
        psi_hat_obj <- make_factor_loadings(
          par = psi0_hat, k = k, q = J, return_all = TRUE
        )
        R_N_hat <- psi_hat_obj$Psi %*% t(psi_hat_obj$Psi)
        diag(R_N_hat) <- 1
      }
    }
  }

  ## Final likelihood. The Hermite optimizer/lower bound is unchanged.
  ## The validated public path uses Genz-style sequential conditional RQMC.
  final_genz_qmc <- NULL
  loglik_approx <- NULL
  loglik_genz_qmc <- NULL
  loglik_final <- -opt$objective

  if (latent.dependence == "gaussian" && latent.loglike == "hermite" && !isTRUE(.skip_final_qmc)) {
    if (!isTRUE(.gh_progressive_internal) && !isTRUE(.gh_node_progressive_internal)) {
      .progress_msg("Evaluating the fitted model likelihood by simulation...")
    }
    loglik_approx <- -opt$objective

    .internal_msg("[dependence] evaluating final Genz-style sequential conditional RQMC likelihood...")
    tg <- system.time(
      final_genz_qmc <- genz_conditional_qmc_loglik_at_fit(
        psi0 = psi0_hat,
        rho0 = rho0_hat
      )
    )["elapsed"]
    final_genz_qmc$elapsed_sec <- unname(tg)
    loglik_genz_qmc <- final_genz_qmc$log_lik_value
    loglik_final <- loglik_genz_qmc

    .internal_msg(paste0(
      "[dependence] Genz conditional RQMC | logLik=", signif(loglik_genz_qmc,10),
      " | SE=", signif(final_genz_qmc$log_lik_se,4),
      " | Q=", final_genz_qmc$Q,
      " | time=", signif(unname(tg),6), " sec"
    ))
  }

  if (latent.dependence == "gaussian" && latent.loglike == "hermite" && isTRUE(.skip_final_qmc)) {
    loglik_approx <- -opt$objective
    loglik_final <- loglik_approx
  }

  ## =========================================================
  ## 13. Zero-layer output
  ## =========================================================

  zero_layer <- if (n_zi == 0L) {

    list(
      type = "No structural-zero layer",
      zero_inflated_margins = integer(0),
      n_zero_inflated = 0L,
      dependence_estimated = FALSE,
      L = 0L,
      rho0 = NULL,
      rho = NULL,
      A = NULL,
      implied_correlation = NULL
    )

  } else if (n_zi == 1L) {

    list(
      type = "Single zero-inflated margin; no zero-dependence model",
      zero_inflated_margins = zi_index,
      n_zero_inflated = 1L,
      dependence_estimated = FALSE,
      L = 0L,
      rho0 = NULL,
      rho = NULL,
      A = NULL,
      implied_correlation = NULL
    )

  } else if (!estimate_B_factor) {

    list(
      type = "Independent Bernoulli structural-zero indicators",
      zero_inflated_margins = zi_index,
      n_zero_inflated = n_zi,
      dependence_estimated = FALSE,
      L = 0L,
      rho0 = NULL,
      rho = NULL,
      A = NULL,
      implied_correlation = diag(n_zi)
    )

  } else {

    list(
      type = paste0(
        if (zero.dependence == "gaussian") "Gaussian" else "BiExp",
        " L-factor copula for Bernoulli structural-zero indicators"
      ),
      zero_inflated_margins = zi_index,
      n_zero_inflated = n_zi,
      dependence_estimated = TRUE,
      dependence = zero.dependence,
      skewness = if (zero.dependence == "biexp") zero.factor.skewness else NULL,
      lambda = if (zero.dependence == "biexp") zero.factor.lambda else NULL,
      L = L,
      quadrature_order = if (use_joint_hermite_allL) NULL else if (zero.dependence == "gaussian") gh_order_B else gl_order_B,
      n_node = if (use_joint_hermite_allL) NULL else n_node_B,
      hermite_joint_core_dimension = if (use_joint_hermite_allL) J + L else NULL,
      hermite_residual_dimension = if (use_joint_hermite_allL) 2L else NULL,
      npar_rho = npar_rho,
      rho0 = rho0_hat,
      rho = rho_hat_obj$Psi,
      A = rho_hat_obj$A,
      implied_correlation = R_B_hat
    )
  }

  ## =========================================================
  ## 14. Return
  ## =========================================================

  # Large Hermite-scaled views are computation-only workspaces.  They are not
  # part of the fitted model and must not be retained in the returned object.
  .tmp_drop <- intersect(
    c("I_y_scaled_list", "log_I_y_scale_list",
      "I_0_scaled_list", "log_I_0_scale_list", "qpb_list"),
    ls(envir = environment())
  )
  if (length(.tmp_drop)) rm(list = .tmp_drop, envir = environment())
  rm(.tmp_drop)
  invisible(gc(FALSE))

  out <- list(
    method = "IFM",
    n = n,
    k = k,
    copula.block_size = copula.block_size,

    latent.dependence = latent.dependence,
    latent.factor.skewness =
      if (latent.dependence == "biexp") latent.factor.skewness else NULL,
    latent.factor.lambda =
      if (latent.dependence == "biexp") latent.factor.lambda else NULL,

    zero.dependence =
      if (estimate_B_factor) zero.dependence else "independent",
    zero.factor.skewness =
      if (estimate_B_factor && zero.dependence == "biexp") {
        zero.factor.skewness
      } else {
        NULL
      },
    zero.factor.lambda =
      if (estimate_B_factor && zero.dependence == "biexp") {
        zero.factor.lambda
      } else {
        NULL
      },

    # Keep only fitted marginal parameters/metadata needed to describe the
    # fitted multivariate model.  The original margin object (and its I cache,
    # Y/X/Z, selection workspaces, etc.) is deliberately not retained.
    marginal = lapply(seq_len(k), function(j) {
      mj <- M[[j]]
      list(
        margin = j,
        distribution = mj$distribution,
        var.power = if (!is.null(mj$var.power)) mj$var.power else mj$var_power,
        zero_inflation = isTRUE(mj$zero_inflation),
        be = mj$be,
        phi = mj$phi,
        gam = if (!is.null(mj$gam)) mj$gam else mj$al,
        log_lik_value = mj$log_lik_value
      )
    }),
    marginal.dist = marginal.dist,
    zero_inflation = zero_inflation,

    zero_layer = zero_layer,

    count_layer = list(
      type = paste0(
        if (latent.dependence == "gaussian") "Gaussian" else "BiExp",
        " J-factor checkerboard copula for mixed-Poisson latent variables"
      ),
      dependence = latent.dependence,
      skewness = if (latent.dependence == "biexp") latent.factor.skewness else NULL,
      lambda = if (latent.dependence == "biexp") latent.factor.lambda else NULL,
      J = J,
      progressive = progressive_used,
      progressive_method = if (progressive_used) latent.loglike else NULL,
      progressive_path = if (progressive_used) progressive_path else J,
      progressive_sample_n = if (progressive_used) length(progressive_sample_index) else n,
      progressive_sample_frac = if (progressive_used) length(progressive_sample_index) / n else 1,
      progressive_sample_index = if (progressive_used && length(progressive_sample_index) < n) {
        progressive_sample_index
      } else {
        NULL
      },
      progressive_history = if (progressive_used) {
        c(
          progressive_history,
          list(list(
            J = J,
            method = latent.loglike,
            n = n,
            logLik = loglik_final,
            convergence = opt$convergence,
            elapsed_sec = unname(elapsed),
            final_full_data = TRUE
          ))
        )
      } else {
        NULL
      },
      d0 = d0,
      loglike_method = if (latent.dependence == "gaussian") latent.loglike else "quadrature",
      quadrature_order = if (latent.dependence == "gaussian" && latent.loglike == "GH") {
        gh_order_N
      } else if (latent.dependence == "biexp") {
        gl_order_N
      } else {
        NULL
      },
      n_node = if (latent.dependence == "gaussian" && latent.loglike == "hermite") NULL else n_node_N,
      npar_psi = npar_psi,
      psi0 = psi0_hat,
      psi = psi_hat_obj$Psi,
      A = psi_hat_obj$A,
      implied_correlation = R_N_hat
    ),

    copula_fit = list(
      neg_log_lik_value = -loglik_final,
      optimization_neg_log_lik_value = opt$objective,
      log_lik_value = loglik_final,
      log_lik_approx = if (latent.dependence == "gaussian" && latent.loglike == "hermite") {
        loglik_approx
      } else {
        NULL
      },
      log_lik_qmc = if (latent.dependence == "gaussian" && latent.loglike == "hermite") {
        NULL
      } else {
        NULL
      },
      log_lik_tilted_qmc = if (latent.dependence == "gaussian" && latent.loglike == "hermite") {
        NULL
      } else {
        NULL
      },
      final_simulation = if (latent.dependence == "gaussian" && latent.loglike == "hermite") final.simulation else NULL,
      qmc = NULL,
      genz_conditional_qmc = final_genz_qmc,
      truncated_antithetic_qmc = NULL,
      tilted_qmc = NULL,
      optimization_loglike_method = if (latent.dependence == "gaussian") latent.loglike else "quadrature",
      convergence = opt$convergence,
      message = opt$message,
      elapsed_sec = unname(elapsed),
      optim = list(
        par = opt$par,
        objective = opt$objective,
        convergence = opt$convergence,
        message = opt$message,
        iterations = opt$iterations,
        evaluations = opt$evaluations
      )
    ),

    call = match.call()
  )

  # Compact user-facing summary. This is an ordinary list (no custom print
  # method) and contains the principal fitted-model results that would
  # otherwise have been reported through progress messages.
  out$summary <- list(
    method = out$method,
    n = n,
    k = k,
    n_mixed_poisson = k - n_zi,
    n_zero_inflated = n_zi,
    marginal = data.frame(
      margin = seq_len(k),
      distribution = marginal.dist,
      var.power = vapply(out$marginal, function(x) as.numeric(x$var.power), numeric(1)),
      zero_inflation = zero_inflation,
      logLik = vapply(out$marginal, function(x) as.numeric(x$log_lik_value), numeric(1)),
      stringsAsFactors = FALSE
    ),
    count_layer = list(
      dependence = out$count_layer$dependence,
      J = out$count_layer$J,
      loglike_method = out$count_layer$loglike_method,
      quadrature_order = out$count_layer$quadrature_order,
      psi = out$count_layer$psi,
      implied_correlation = out$count_layer$implied_correlation,
      progressive = out$count_layer$progressive,
      progressive_path = out$count_layer$progressive_path
    ),
    zero_layer = list(
      dependence = out$zero.dependence,
      L = out$zero_layer$L,
      dependence_estimated = out$zero_layer$dependence_estimated,
      rho = out$zero_layer$rho,
      implied_correlation = out$zero_layer$implied_correlation
    ),
    fit = list(
      logLik = out$copula_fit$log_lik_value,
      optimization_logLik = -out$copula_fit$optimization_neg_log_lik_value,
      convergence = out$copula_fit$convergence,
      message = out$copula_fit$message,
      elapsed_sec = out$copula_fit$elapsed_sec,
      final_simulation = out$copula_fit$final_simulation,
      genz_conditional_qmc = out$copula_fit$genz_conditional_qmc
    )
  )

  class(out) <- c("multi_zipm_fit", class(out))

  if (!isTRUE(.gh_progressive_internal) && !isTRUE(.gh_node_progressive_internal)) {
    if (identical(latent.dependence, "gaussian") && identical(latent.loglike, "hermite") &&
        !isTRUE(.skip_final_qmc) && is.finite(loglik_final)) {
      .progress_msg("  Simulation-evaluated log-likelihood: ", signif(loglik_final, 10))
    }
    .progress_msg("Factor-copula estimation finished.\n")
  }

  out
}



#' Fit a Multivariate Zero-Inflated Mixed Poisson Model
#'
#' @description
#' Fits a multivariate Mixed Poisson (MMP) or Zero-Inflated Mixed Poisson
#' (MZIMP) model with factor-copula dependence.
#'
#' Marginal models may either be fitted directly from `Y`, `X`, and
#' optionally `Z`, or supplied through an object returned by
#' [select_multi_margin()]. After the marginal models have been fitted or
#' supplied, `fit_multizimp()` estimates dependence among the latent
#' mixed-Poisson components and, when requested, dependence among
#' structural-zero indicators.
#'
#' The marginal and dependence parameters are estimated using an
#' inference-functions-for-margins (IFM) strategy.
#'
#' @details
#' For margin \eqn{j=1,\ldots,k}, the observed response is represented as
#'
#' \deqn{
#' Y_j = B_j N_j,
#' }
#'
#' where \eqn{N_j} follows a mixed Poisson model and \eqn{B_j} is an
#' activation indicator when zero inflation is present. For an ordinary
#' mixed Poisson margin, \eqn{B_j=1} with probability one.
#'
#' The model separates two sources of cross-margin dependence:
#'
#' \enumerate{
#'   \item dependence among the latent mixed-Poisson components
#'   underlying \eqn{N_1,\ldots,N_k}; and
#'   \item dependence among the activation indicators for zero-inflated
#'   margins.
#' }
#'
#' These two dependence structures are modeled separately through the
#' count layer and the structural-zero layer.
#'
#' @section Marginal models:
#'
#' Marginal models can be supplied in two ways.
#'
#' \strong{Direct specification.}
#'
#' If `margin = NULL`, `Y` and `X` must be supplied. The function fits
#' exactly one marginal model for each response. No marginal model
#' selection is performed.
#'
#' `mixing.distribution` specifies the mixing distribution and
#' `var.power` specifies its variance-power parameterization.
#'
#' `Z = NULL` gives ordinary mixed Poisson margins. Supplying `Z`
#' gives zero-inflated mixed Poisson margins.
#'
#' `X` and `Z` may be common matrices/data frames shared by all margins
#' or lists containing margin-specific designs. When lists are supplied,
#' their lengths must equal the number of response margins.
#'
#' \strong{Previously selected margins.}
#'
#' Alternatively, `margin` may be an object returned by
#' [select_multi_margin()]. In this case the selected marginal
#' distributions, variance powers, regression coefficients, dispersion
#' parameters, and zero-inflation specifications are reused without
#' refitting or reselecting the marginal models.
#'
#' If compatible checkerboard probability matrices were precomputed by
#' [select_multi_margin()], they are reused when possible.
#'
#' A `margin` object may contain a mixture of ordinary MP and ZIMP
#' margins.
#'
#' @section Count dependence:
#'
#' Dependence among the latent mixed-Poisson components is represented
#' by a `J`-factor copula. The number of factors must satisfy
#'
#' \deqn{
#' 1 \leq J < k.
#' }
#'
#' The factor family is controlled by `latent.dependence`:
#'
#' \describe{
#'   \item{`"gaussian"`}{
#'     Gaussian factor copula.
#'   }
#'   \item{`"biexp"`}{
#'     Factor copula based on a standardized bilateral-exponential
#'     distribution.
#'   }
#' }
#'
#' The fitted count-layer loading matrix is denoted by \eqn{\Psi} and has
#' dimension \eqn{k \times J}. The parameterization guarantees
#'
#' \deqn{
#' \sum_{r=1}^{J}\psi_{jr}^2 < 1
#' }
#'
#' for every margin \eqn{j}, leaving a positive idiosyncratic variance.
#'
#' For Gaussian factors, `latent.loglike` controls likelihood evaluation.
#' `"GH"` uses Gauss--Hermite quadrature and `"hermite"` uses the
#' Hermite-based approximation implemented by the package. If
#' `latent.loglike = NULL`, the method is selected automatically.
#'
#' @section Structural-zero dependence:
#'
#' Zero-inflated margins may additionally share a dependence structure
#' through their activation indicators.
#'
#' `L = 0` specifies no factor dependence in the structural-zero layer.
#' `L = 1` requests a one-factor zero-layer dependence model when at least
#' two zero-inflated margins are available.
#'
#' The zero-layer factor family is specified by `zero.dependence`, with
#' supported values `"gaussian"` and `"biexp"`.
#'
#' If no margin is zero inflated, or if only one margin is zero inflated,
#' no cross-margin zero-layer dependence can be estimated. Ordinary MP
#' margins do not participate in the structural-zero dependence model.
#'
#' Thus a multivariate model may contain any combination of MP and ZIMP
#' margins.
#'
#' @section Checkerboard probabilities:
#'
#' Dependence estimation uses checkerboard interval probabilities derived
#' from the fitted marginal distributions. The checkerboard resolution is
#' controlled by `d0`.
#'
#' These marginal quantities are computed before dependence optimization
#' and do not depend on the factor-loading parameters. When a compatible
#' cache is available in an object returned by [select_multi_margin()],
#' it is reused where possible.
#'
#' @section Estimation:
#'
#' The model is fitted using an IFM strategy. Marginal parameters are
#' estimated first, or taken from `margin`, and are subsequently held
#' fixed while the dependence parameters are estimated.
#'
#' For Gaussian count factors, `gh.progressive = TRUE` allows
#' lower-dimensional fitted factor models to be used as starting values
#' for higher-dimensional models.
#'
#' `parallel` controls supported parallel computations. `NULL` or `1`
#' requests serial computation; an integer greater than one requests the
#' corresponding level of parallelism where supported.
#'
#' @param Y Optional \eqn{n \times k} numeric matrix or data frame of
#'   multivariate count responses. Required when `margin = NULL`.
#'   Entries must be finite, non-negative integer counts and `k` must be
#'   at least two.
#'
#' @param X Count-component design used when `margin = NULL`. May be a
#'   matrix/data frame shared by all margins or a list of length `k`
#'   containing margin-specific designs. Intercepts must be supplied
#'   explicitly when required.
#'
#' @param Z Optional zero-component design used when `margin = NULL`.
#'   `NULL` specifies ordinary mixed Poisson margins. A matrix/data frame
#'   specifies a common zero-component design, while a list may provide
#'   margin-specific designs.
#'
#' @param margin Optional object returned by [select_multi_margin()].
#'   When supplied, the selected marginal models are reused and the
#'   direct marginal-fitting arguments `Y`, `X`, `Z`,
#'   `mixing.distribution`, `var.power`, and `margin.ini.par` are not used
#'   to reselect the margins.
#'
#' @param mixing.distribution Mixing distribution used for direct
#'   marginal fitting when `margin = NULL`. A single value may be shared
#'   by all margins or a vector may specify one distribution per margin.
#'   Supported distributions are `"gamma"`, `"lnorm"`, `"ig"`,
#'   `"weibull"`, `"loglogistic"`, and `"betaprime"`.
#'
#' @param var.power Optional variance-power specification used for direct
#'   marginal fitting.
#'
#' @param margin.ini.par Optional initial parameter specification for
#'   direct marginal fitting.
#'
#' @param parallel `NULL` or a positive integer controlling supported
#'   parallel computation. `NULL` or `1` requests serial execution.
#'
#' @param margin.block_size Optional positive integer controlling the
#'   block size used during parallel marginal fitting.
#'
#' @param copula.block_size Positive integer giving the number of
#'   observations processed in each dependence-likelihood block.
#'   The default is `500L`.
#'
#' @param d0 Integer greater than or equal to `2` specifying the
#'   checkerboard resolution used in the count dependence layer.
#'   The default is `10`.
#'
#' @param J Positive integer giving the number of factors in the
#'   count/intensity dependence layer. It must satisfy `1 <= J < k`.
#'   The default is `1L`.
#'
#' @param L Integer specifying the number of factors in the
#'   structural-zero dependence layer. The current implementation
#'   supports `0` and `1`. `L = 0` disables zero-layer factor dependence;
#'   `L = 1` requests a one-factor model when estimable.
#'
#' @param latent.dependence Character string specifying the count-layer
#'   factor family. Supported values are `"gaussian"` and `"biexp"`.
#'
#' @param latent.factor.skewness Numeric parameter controlling skewness
#'   of the standardized bilateral-exponential factor distribution when
#'   `latent.dependence = "biexp"`.
#'
#' @param latent.loglike Likelihood method for Gaussian count factors.
#'   May be `NULL`, `"GH"`, or `"hermite"`. If `NULL`, the method is
#'   selected automatically according to the factor dimension.
#'
#' @param gh_order Integer greater than or equal to `2` specifying the
#'   quadrature order used by dependence calculations where applicable.
#'   The default is `16L`.
#'
#' @param approx.hermite.joint.order Hermite approximation order used by
#'   the joint count/zero likelihood calculation where applicable.
#'   The default is `4L`.
#'
#' @param hermite.polar.cycles Positive integer giving the maximum number
#'   of alternating direction/radius cycles used by the polar Hermite
#'   optimizer.
#'
#' @param hermite.polar.v.iter Positive integer giving the maximum number
#'   of optimizer iterations for the direction step of each polar cycle.
#'
#' @param hermite.polar.r.iter Positive integer giving the maximum number
#'   of optimizer iterations for the radius step of each polar cycle.
#'
#' @param hermite.polar.tol Positive numeric convergence tolerance for the
#'   polar Hermite optimizer.
#'
#' @param hermite.polar.r.max Numeric value in `(0, 1)` giving the maximum
#'   admissible row norm for factor loadings in the polar
#'   parameterization.
#'
#' @param genz.conditional.Q Positive integer giving the number of
#'   conditional randomized quasi-Monte Carlo paths used for final
#'   Genz likelihood evaluation where applicable.
#'
#' @param genz.conditional.antithetic Logical indicating whether
#'   antithetic uniforms are used in the conditional RQMC evaluation.
#'
#' @param qmc.n_rep Positive integer giving the number of randomized QMC
#'   replicates used for final likelihood evaluation.
#'
#' @param qmc.seed Integer seed used for randomized QMC calculations and
#'   other reproducible stochastic computations within the fitting
#'   procedure.
#'
#' @param zero.dependence Character string specifying the factor family
#'   for the structural-zero layer when `L = 1`. Supported values are
#'   `"gaussian"` and `"biexp"`.
#'
#' @param zero.factor.skewness Numeric parameter controlling skewness of
#'   the standardized bilateral-exponential factor distribution when
#'   `zero.dependence = "biexp"`.
#'
#' @param biexp.root_tol Positive numeric tolerance used in numerical
#'   root finding required by the bilateral-exponential factor backend.
#'
#' @param psi0_start Optional numeric vector of raw starting parameters
#'   for the count-layer factor-loading matrix.
#'
#' @param rho0_start Optional numeric vector of raw starting parameters
#'   for the structural-zero factor-loading matrix.
#'
#' @param gh.progressive Logical. If `TRUE`, supported Gaussian
#'   count-factor fits may use progressively fitted lower-dimensional
#'   factor models as warm starts.
#'
#' @param use.gradient Logical indicating whether available analytic
#'   gradients are used during dependence optimization.
#'
#' @param progress Deprecated and ignored. Retained for backward API
#'   compatibility; public fitting is currently silent.
#'
#' @param trace Deprecated and ignored. Retained for backward API
#'   compatibility.
#'
#' @param eval.max Positive integer giving the maximum number of
#'   objective-function evaluations for the main dependence optimizer.
#'
#' @param iter.max Positive integer giving the maximum number of
#'   iterations for the main dependence optimizer.
#'
#' @param tiny Positive numeric lower bound used to protect probabilities
#'   before logarithms are evaluated.
#'
#' @param big_penalty Positive numeric objective value returned for
#'   invalid parameter states during optimization.
#'
#' @param rel.tol Positive numeric relative convergence tolerance for the
#'   main dependence optimizer.
#'
#' @return
#' An object of class `"multi_zipm_fit"` containing the fitted marginal
#' models and dependence structure. Important components include:
#'
#' \describe{
#'   \item{`method`}{
#'     Estimation method. Currently `"IFM"`.
#'   }
#'
#'   \item{`n`, `k`}{
#'     Number of observations and response margins.
#'   }
#'
#'   \item{`marginal`}{
#'     List of length `k` containing the fitted marginal information.
#'     Each element contains the margin index, mixing distribution,
#'     variance power, zero-inflation status, count coefficients (`be`),
#'     dispersion parameter (`phi`), zero-component coefficients (`gam`)
#'     when present, and marginal log-likelihood.
#'   }
#'
#'   \item{`marginal.dist`}{
#'     Mixing distribution used for each response margin.
#'   }
#'
#'   \item{`zero_inflation`}{
#'     Logical vector indicating which margins contain a zero-inflation
#'     component.
#'   }
#'
#'   \item{`count_layer`}{
#'     Fitted dependence structure for the mixed-Poisson count layer.
#'     Important components include `dependence`, `J`, `psi`, `A`,
#'     `implied_correlation`, and information on the likelihood method.
#'     The loading matrix `psi` has dimension `k` by `J`.
#'   }
#'
#'   \item{`zero_layer`}{
#'     Structural-zero dependence information. When zero-layer dependence
#'     is estimated, this includes the fitted loading matrix `rho` and its
#'     implied correlation matrix. The component also records cases in
#'     which no zero-layer dependence is estimated.
#'   }
#'
#'   \item{`copula_fit`}{
#'     Dependence-fitting results, including the fitted log-likelihood,
#'     optimization objective, convergence information, elapsed time,
#'     optimizer output, and available final likelihood-evaluation
#'     diagnostics.
#'   }
#'
#'   \item{`summary`}{
#'     Ordinary list containing the principal user-facing results,
#'     including dimensions, marginal-model summary, count- and
#'     zero-layer summaries, log-likelihood, convergence information,
#'     elapsed time, and final likelihood diagnostics.
#'   }
#'
#'   \item{`call`}{
#'     Matched function call.
#'   }
#' }
#'
#' @seealso
#' [select_multi_margin()] for selecting marginal models before
#' dependence fitting; [fit_zimp()] and [select_zimp()] for univariate
#' models; [fit_bizimp()] for the specialized bivariate model.
#'
#' @examples
#' set.seed(12345)
#'
#' n <- 15
#' k <- 3
#' x <- rnorm(n)
#' X <- cbind(Intercept = 1, x = x)
#'
#' beta <- c(0.5, 0.3)
#' mu <- exp(X %*% beta)
#'
#' Y <- cbind(
#'   y1 = rpois(n, mu),
#'   y2 = rpois(n, mu),
#'   y3 = rpois(n, mu)
#' )
#'
#' ## Directly specified mixed-Poisson margins
#' fit <- fit_multizimp(
#'   Y = Y,
#'   X = X,
#'   Z = NULL,
#'   mixing.distribution = "gamma",
#'   var.power = 2,
#'   J = 1,
#'   L = 0,
#'   d0 = 3,
#'   parallel = 1
#' )
#'
#' fit$summary
#' fit$count_layer$psi
#'
#' ## Alternatively, select the marginal models first
#' margins <- select_multi_margin(
#'   Y = Y,
#'   X = X,
#'   Z = NULL,
#'   mixing.distribution = c("gamma", "lnorm"),
#'   var.power = 2,
#'   compute.I = TRUE,
#'   I.d0 = 3,
#'   parallel = 1
#' )
#'
#' fit2 <- fit_multizimp(
#'   margin = margins,
#'   J = 1,
#'   L = 0,
#'   d0 = 3,
#'   parallel = 1
#' )
#' fit2$summary
#' @export
# ============================================================
# Backward-compatible convenience wrapper
# ============================================================
fit_multizimp <- function(
    # Code mới	
    # Response count matrix. Required only when `margin = NULL`.
    Y = NULL,
    # Count-model design matrix/list used during marginal fitting.
    X = NULL,
    # Zero-inflation design matrix/list used during marginal fitting.
    Z = NULL,
    # Pre-fitted object from select_multi_margin(); avoids refitting margins.
    margin = NULL,
    # Mixing distribution(s) for the marginal models. A scalar is shared
    # by all margins; otherwise supply one value per margin.
    mixing.distribution = "gamma",
    # Optional variance-power specification for marginal models.
    var.power = NULL,
    # Optional initial parameters for marginal fitting.
    margin.ini.par = NULL,
    # Parallel convention shared by marginal fitting and GH dependence fitting:
    # NULL or 1 = serial; integer > 1 = that many workers/threads.
    # With direct Y/X/Z input, the same value is used for marginal fitting.
    # With a pre-selected `margin`, marginal models are reused and the value
    # controls only the dependence GH backend.
    parallel = NULL,
    # Optional block size used by marginal fitting.
    margin.block_size = NULL,
    # Observation block size used by dependence-likelihood calculations.
    copula.block_size = 500L,
    # Checkerboard resolution for the latent mixed-Poisson layer.
    d0 = 10,
    # Number of Gaussian factors for the count/intensity layer.
    J = 1L,
    # Number of Gaussian factors for the structural-zero layer; currently 0 or 1.
    L = 0L,
    # Factor family for the count/intensity layer.
    latent.dependence = c("gaussian", "biexp"),
    # Skewness parameter used only for the BiExp count-factor family.
    latent.factor.skewness = 0.5,
    # Likelihood engine for Gaussian factors: NULL chooses GH for J<=2 and Hermite for J>=3.
    latent.loglike = NULL,
    # Single quadrature order used throughout the dependence fit:
    # GH count layer, GH zero layer (when present), Hermite projection,
    # and Gauss-Legendre quadrature for BiExp backends.
    gh_order = 16L,
    # Hermite polynomial order reserved for the joint count/zero branch; currently 2 or 4.
    approx.hermite.joint.order = 4L,
    # Hermite/Genz threads are derived from `parallel`; there is no separate
    # public thread-control argument in v1.
    # Number of alternating polar V|r -> r|V cycles in the Hermite optimizer.
    hermite.polar.cycles = 6L,
    # Maximum optimizer iterations for the direction V step of each polar cycle.
    hermite.polar.v.iter = 8L,
    # Maximum optimizer iterations for the radius r step of each polar cycle.
    hermite.polar.r.iter = 4L,
    # Relative convergence tolerance for the polar Hermite optimizer.
    hermite.polar.tol = 1e-5,
    # Maximum admissible row radius ||psi_j|| in the polar parameterization.
    hermite.polar.r.max = 0.995,
    # Number of conditional-RQMC paths used by the final Genz likelihood evaluation.
    genz.conditional.Q = 256L,
    # Whether antithetic uniforms are added to the final Genz conditional RQMC.
    genz.conditional.antithetic = TRUE,
    # Number of randomized QMC replicates used to estimate final likelihood uncertainty.
    qmc.n_rep = 4L,
    # Random seed used for QMC randomization and progressive sampling.
    qmc.seed = 1L,
    # Dependence family for the structural-zero factor when L=1.
    zero.dependence = c("gaussian", "biexp"),
    # Skewness parameter used only for a BiExp structural-zero factor.
    zero.factor.skewness = 0.5,
    # Root-finding tolerance used for BiExp quantiles.
    biexp.root_tol = 1e-10,
    # Optional raw starting parameters for the count-factor loading matrix Psi.
    psi0_start = NULL,
    # Optional raw starting parameters for the structural-zero loading matrix Rho.
    rho0_start = NULL,
    # Whether progressive factor-dimension warm starts are used.
    gh.progressive = TRUE,
    # v1 warm-starts every lower-J model on the FULL data.  No observation
    # subsampling/sample-size ladder is used.  These controls are retained
    # internally below only because the validated backend still accepts them.
    # Whether available analytic gradients are used during optimization.
    use.gradient = TRUE,
    # Whether progress messages are printed.
    progress = TRUE,
    # Optimizer trace level.
    trace = 0,
    # Maximum objective evaluations for the main optimizer.
    eval.max = 1000L,
    # Maximum iterations for the main optimizer.
    iter.max = 100L,
    # Numerical lower bound used to protect probabilities before taking logs.
    tiny = 1e-300,
    # Objective value returned when an invalid parameter state is encountered.
    big_penalty = 1e12,
    # Relative convergence tolerance for the main optimizer.
    rel.tol = 1e-4
 ) {
  # Public fitting is intentionally silent. These arguments are retained only
  # for backward API compatibility.
  progress <- FALSE
  trace <- 0

  # ---------------------------------------------------------------
  # v1 parallel convention: identical to select_multi_margin().
  # NULL / 1 -> serial; integer > 1 -> parallel.
  # The same requested parallelism is used for direct marginal fits
  # and for the exact GH dependence backend.
  # ---------------------------------------------------------------
  if (is.null(parallel)) {
    n_parallel <- 1L
  } else {
    if (!is.numeric(parallel) || length(parallel) != 1L ||
        !is.finite(parallel) || parallel < 1L ||
        abs(parallel - round(parallel)) > 1e-12) {
      stop("`parallel` must be NULL or a positive integer.")
    }
    n_parallel <- as.integer(round(parallel))
  }
  if (!is.numeric(gh_order) || length(gh_order) != 1L ||
      !is.finite(gh_order) || gh_order < 2L ||
      abs(gh_order - round(gh_order)) > 1e-12) {
    stop("`gh_order` must be a single integer >= 2.")
  }
  gh_order <- as.integer(round(gh_order))

  # One public quadrature-order control. Internal names are retained only
  # because the validated lower-level backends still use these arguments.
  gh_order_N <- gh_order
  gh_order_B <- gh_order
  approx.hermite.gh_order <- gh_order
  approx.hermite.joint.gh_order <- gh_order
  gl_order_N <- gh_order
  gl_order_B <- gh_order

  # A single public `parallel` argument controls every threaded dependence
  # backend in v1.  Keep the internal controls separate only for the lower-level
  # .fit_multizimp_copula() API.
  gh.threads <- n_parallel
  approx.hermite.threads <- n_parallel
  margin.n_workers <- n_parallel
  margin.parallel <- n_parallel > 1L

  prepared <- prepare_multi_zipm_I(
    Y = Y, X = X, Z = Z, margin = margin,
    mixing.distribution = mixing.distribution,
    var.power = var.power,
    margin.ini.par = margin.ini.par,
    margin.n_workers = margin.n_workers,
    margin.parallel = margin.parallel,
    margin.block_size = margin.block_size,
    d0 = d0,
    progress = progress
  )

  ## Internal constants retained only because validated legacy kernels still
  ## reference them. They are no longer public tuning options in this branch.
  qmc.Q <- 2^12
  qmc.chunk_size <- 4096L
  final.simulation <- "genz_conditional_qmc"
  simulation.rotation.energy <- 0.99
  simulation.eigen.tol <- 1e-4
  simulation.min.rank <- 1L
  simulation.max.rank <- NULL
  tilted.groups <- 8L
  tilted.cov.scale <- 1.25
  tilted.between.scale <- 1.0
  tilted.min.eigen <- 1e-5

  .fit_multizimp_copula(
    prepared = prepared,
    copula.block_size = copula.block_size,
    d0 = d0,
    J = J,
    L = L,
    latent.dependence = latent.dependence,
    latent.factor.skewness = latent.factor.skewness,
    latent.loglike = latent.loglike,
    gh_order_N = gh_order_N,
    gh.threads = gh.threads,
    approx.hermite.gh_order = approx.hermite.gh_order,
    approx.hermite.joint.gh_order = approx.hermite.joint.gh_order,
    approx.hermite.joint.order = approx.hermite.joint.order,
    approx.hermite.threads = approx.hermite.threads,
    hermite.optimizer = "polar",
    hermite.polar.cycles = hermite.polar.cycles,
    hermite.polar.v.iter = hermite.polar.v.iter,
    hermite.polar.r.iter = hermite.polar.r.iter,
    hermite.polar.tol = hermite.polar.tol,
    hermite.polar.r.max = hermite.polar.r.max,
    hermite.polar.joint.refine.iter = 0L,
    qmc.Q = qmc.Q,
    qmc.n_rep = qmc.n_rep,
    qmc.seed = qmc.seed,
    qmc.chunk_size = qmc.chunk_size,
    final.simulation = final.simulation,
    genz.conditional.Q = genz.conditional.Q,
    genz.conditional.antithetic = genz.conditional.antithetic,
    simulation.rotation.energy = simulation.rotation.energy,
    simulation.eigen.tol = simulation.eigen.tol,
    simulation.min.rank = simulation.min.rank,
    simulation.max.rank = simulation.max.rank,
    tilted.groups = tilted.groups,
    tilted.cov.scale = tilted.cov.scale,
    tilted.between.scale = tilted.between.scale,
    tilted.min.eigen = tilted.min.eigen,
    zero.dependence = zero.dependence,
    zero.factor.skewness = zero.factor.skewness,
    gh_order_B = gh_order_B,
    gl_order_N = gl_order_N,
    gl_order_B = gl_order_B,
    biexp.root_tol = biexp.root_tol,
    psi0_start = psi0_start,
    rho0_start = rho0_start,
    gh.progressive = gh.progressive,
    # v1: all progressive J stages use the full sample.
    progressive.sample.frac = 1,
    progressive.sample.n = NULL,
    progressive.sample.seed = qmc.seed,
    progressive.target.iter = 0L,
    progressive.target.ladder = FALSE,
    progressive.target.ladder.frac = numeric(0),
    progressive.target.ladder.iter = 0L,
    .gh_progressive_internal = FALSE,
    .gh_node_progressive_internal = FALSE,
    .skip_final_qmc = FALSE,
    use.gradient = use.gradient,
    progress = progress,
    trace = trace,
    eval.max = eval.max,
    iter.max = iter.max,
    tiny = tiny,
    big_penalty = big_penalty,
    rel.tol = rel.tol
  )
}

# ============================================================
# Source/version guard: progressive Hermite Q >= 2
# ============================================================
.fit_multizimp_progressive_Q_version <- "2026-09-26-Q2-v2"

.check_multizimp_progressive_Q_source <- function(verbose = TRUE) {
  b <- paste(deparse(body(.fit_multizimp_copula)), collapse = "\n")
  old5 <- grepl("approx\\.hermite\\.gh_order.*>= 5", b) ||
          grepl("approx\\.hermite\\.gh_order.*< 5L", b)
  q2 <- grepl("approx\\.hermite\\.gh_order.*< 2L", b)
  jq2 <- grepl("approx\\.hermite\\.joint\\.gh_order.*< 2L", b)
  ok <- !old5 && q2 && jq2
  if (!ok) {
    stop("The active `.fit_multizimp_copula()` is not the Q>=2 progressive-Hermite version. Re-source this file AFTER all older fit_multizimp/helper files.")
  }
  invisible(ok)
}

.check_multizimp_progressive_Q_source(verbose = TRUE)
