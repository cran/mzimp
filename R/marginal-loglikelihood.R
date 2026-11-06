neg_loglik_marginal <- function(
    y,
    X,
    Z = NULL,
    dist = c(
      "gamma",
      "lnorm",
      "ig",
      "weibull",
      "loglogistic",
      "betaprime"
    ),
    be,
    phi,
    al = NULL,
    zero_inflation = FALSE,
    use_custom_var_power = FALSE,
    var_power = NULL,
    n_gh = 30,
    gh = NULL,
    eps_prob = 1e-300
) {
  
  dist <- match.arg(dist)
  
  # ------------------------------------------------------------
  # Optional variance-power argument
  # ------------------------------------------------------------
  
  arg_vp <- if (use_custom_var_power) {
    
    if (is.null(var_power)) {
      stop(
        "var_power must be supplied when ",
        "use_custom_var_power = TRUE."
      )
    }
    
    list(
      var_power = var_power
    )
    
  } else {
    
    list()
  }
  
  # ------------------------------------------------------------
  # Gauss-Hermite quadrature for lognormal only
  # ------------------------------------------------------------
  
  arg_gh <- list()
  
  if (dist == "lnorm") {
    
    if (is.null(gh)) {
      
      gh <- statmod::gauss.quad(
        n = n_gh,
        kind = "hermite"
      )
    }
    
    arg_gh <- list(
      gh = gh
    )
  }
  
  # ------------------------------------------------------------
  # Zero-inflated model
  # ------------------------------------------------------------
  
  if (zero_inflation) {
    
    if (is.null(Z)) {
      stop(
        "Z must be supplied when zero_inflation = TRUE."
      )
    }
    
    if (is.null(al)) {
      stop(
        "al must be supplied when zero_inflation = TRUE."
      )
    }
    
    if (dist == "gamma") {
      
      re <- do.call(
        ZI_gamma,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp
        )
      )
      
    } else if (dist == "lnorm") {
      
      re <- do.call(
        ZI_lnorm,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp,
          arg_gh
        )
      )
      
    } else if (dist == "ig") {
      
      re <- do.call(
        ZI_ig,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp
        )
      )
      
    } else if (dist == "weibull") {
      
      re <- do.call(
        ZI_weibull,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp
        )
      )
      
    } else if (dist == "loglogistic") {
      
      re <- do.call(
        ZI_loglogistic,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp
        )
      )
      
    } else if (dist == "betaprime") {
      
      re <- do.call(
        ZI_betaprime,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp
        )
      )
    }
    
  } else {
    
    # ----------------------------------------------------------
    # Non-zero-inflated model
    # ----------------------------------------------------------
    
    if (dist == "gamma") {
      
      re <- do.call(
        I_gamma,
        c(
          list(
            z = y,
            X = X,
            be = be,
            phi = phi,
            d = 1L
          ),
          arg_vp
        )
      )
      
    } else if (dist == "lnorm") {
      
      re <- do.call(
        I_lnorm,
        c(
          list(
            z = y,
            X = X,
            be = be,
            phi = phi,
            d = 1L
          ),
          arg_vp,
          arg_gh
        )
      )
      
    } else if (dist == "ig") {
      
      re <- do.call(
        I_ig,
        c(
          list(
            z = y,
            X = X,
            be = be,
            phi = phi,
            d = 1L
          ),
          arg_vp
        )
      )
      
    } else if (dist == "weibull") {
      
      re <- do.call(
        I_weibull,
        c(
          list(
            z = y,
            X = X,
            be = be,
            phi = phi,
            d = 1L
          ),
          arg_vp
        )
      )
      
    } else if (dist == "loglogistic") {
      
      re <- do.call(
        I_loglogistic,
        c(
          list(
            z = y,
            X = X,
            be = be,
            phi = phi,
            d = 1L
          ),
          arg_vp
        )
      )
      
    } else if (dist == "betaprime") {
      
      re <- do.call(
        I_betaprime,
        c(
          list(
            z = y,
            X = X,
            be = be,
            phi = phi,
            d = 1L
          ),
          arg_vp
        )
      )
    }
  }
  
  # ------------------------------------------------------------
  # Negative log-likelihood
  # ------------------------------------------------------------
  
  re <- pmax(
    as.numeric(re),
    eps_prob
  )
  
  -sum(
    log(re)
  )
}


neg_loglik_marginal_grad <- function(
    y,
    X,
    Z = NULL,
    dist = c(
      "gamma",
      "lnorm",
      "ig",
      "weibull",
      "loglogistic",
      "betaprime"
    ),
    be,
    phi,
    al = NULL,
    zero_inflation = FALSE,
    use_custom_var_power = FALSE,
    var_power = NULL,
    n_gh = 20,
    gh = NULL,
    eps_prob = 1e-300
) {
  
  dist <- match.arg(dist)
  
  # ------------------------------------------------------------
  # Optional variance-power argument
  # ------------------------------------------------------------
  
  arg_vp <- if (use_custom_var_power) {
    
    if (is.null(var_power)) {
      stop(
        "var_power must be supplied when ",
        "use_custom_var_power = TRUE."
      )
    }
    
    list(
      var_power = var_power
    )
    
  } else {
    
    list()
  }
  
  # ------------------------------------------------------------
  # Gauss-Hermite quadrature for lognormal only
  # ------------------------------------------------------------
  
  arg_gh <- list()
  
  if (dist == "lnorm") {
    
    if (is.null(gh)) {
      
      gh <- statmod::gauss.quad(
        n = n_gh,
        kind = "hermite"
      )
    }
    
    arg_gh <- list(
      gh = gh
    )
  }
  
  # ============================================================
  # Zero-inflated case
  # ============================================================
  
  if (zero_inflation) {
    
    if (is.null(Z)) {
      stop(
        "Z must be supplied when zero_inflation = TRUE."
      )
    }
    
    if (is.null(al)) {
      stop(
        "al must be supplied when zero_inflation = TRUE."
      )
    }
    
    if (dist == "gamma") {

      # ZI_gamma_grad() returns both probability and gradient.
      # Therefore ZI_gamma() does not need to be evaluated separately.
      gr <- do.call(
        ZI_gamma_grad,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp
        )
      )

      re <- gr$prob

    } else if (dist == "lnorm") {

      # ZI_lnorm_grad() returns both probability and gradient.
      # Therefore ZI_lnorm() does not need to be evaluated separately.
      gr <- do.call(
        ZI_lnorm_grad,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp,
          arg_gh
        )
      )

      re <- gr$prob

    } else if (dist == "ig") {

      # ZI_ig_grad() returns both probability and gradient.
      # Therefore ZI_ig() does not need to be evaluated separately.
      gr <- do.call(
        ZI_ig_grad,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp
        )
      )

      re <- gr$prob

    } else if (dist == "weibull") {

      # ZI_weibull_grad() returns both probability and gradient.
      # Therefore ZI_weibull() does not need to be evaluated separately.
      gr <- do.call(
        ZI_weibull_grad,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp
        )
      )

      re <- gr$prob

    } else if (dist == "loglogistic") {

      # ZI_loglogistic_grad() returns both probability and gradient.
      # Therefore ZI_loglogistic() does not need to be evaluated separately.
      gr <- do.call(
        ZI_loglogistic_grad,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp
        )
      )

      re <- gr$prob

    } else if (dist == "betaprime") {

      # ZI_betaprime_grad() returns the probability in $value.
      # Add $prob so that beta-prime follows the same interface as
      # the other mixing distributions in this function.
      gr <- do.call(
        ZI_betaprime_grad,
        c(
          list(
            y = y,
            X = X,
            Z = Z,
            be = be,
            phi = phi,
            al = al
          ),
          arg_vp
        )
      )

      gr$prob <- gr$value
      re <- gr$prob
    }

    # ----------------------------------------------------------
    # Gradient of negative log-likelihood
    # ----------------------------------------------------------
    
    re <- pmax(
      as.numeric(re),
      eps_prob
    )
    
    l_phi <-
      -sum(
        as.numeric(gr$d_phi) /
        re
      )
    
    l_be <-
      colSums(
        -gr$d_be /
        re
      )
    
    l_al <-
      colSums(
        -gr$d_al /
        re
      )
    
    return(
      list(
        d_be = l_be,
        d_phi = l_phi,
        d_al = l_al
      )
    )
  }
  
  # ============================================================
  # Non-zero-inflated case
  # ============================================================
  
  if (dist == "gamma") {

    # I_gamma_grad() returns both probability and gradient.
    # Therefore I_gamma() does not need to be evaluated separately.
    gr <- do.call(
      I_gamma_grad,
      c(
        list(
          z = y,
          X = X,
          be = be,
          phi = phi,
          d = 1L
        ),
        arg_vp
      )
    )

    re <- gr$prob

  } else if (dist == "lnorm") {

    # I_lnorm_grad() returns both probability and gradient.
    # Therefore I_lnorm() does not need to be evaluated separately.
    gr <- do.call(
      I_lnorm_grad,
      c(
        list(
          z = y,
          X = X,
          be = be,
          phi = phi,
          d = 1L
        ),
        arg_vp,
        arg_gh
      )
    )

    re <- gr$prob

  } else if (dist == "ig") {

    # I_ig_grad() returns both probability and gradient.
    # Therefore I_ig() does not need to be evaluated separately.
    gr <- do.call(
      I_ig_grad,
      c(
        list(
          z = y,
          X = X,
          be = be,
          phi = phi,
          d = 1L
        ),
        arg_vp
      )
    )

    re <- gr$prob

  } else if (dist == "weibull") {

    # I_weibull_grad() returns both probability and gradient.
    # Therefore I_weibull() does not need to be evaluated separately.
    gr <- do.call(
      I_weibull_grad,
      c(
        list(
          z = y,
          X = X,
          be = be,
          phi = phi,
          d = 1L
        ),
        arg_vp
      )
    )

    re <- gr$prob

  } else if (dist == "loglogistic") {

    # I_loglogistic_grad() returns both probability and gradient.
    # Therefore I_loglogistic() does not need to be evaluated separately.
    gr <- do.call(
      I_loglogistic_grad,
      c(
        list(
          z = y,
          X = X,
          be = be,
          phi = phi,
          d = 1L
        ),
        arg_vp
      )
    )

    re <- gr$prob

  } else if (dist == "betaprime") {

    # I_betaprime_grad() returns the probability in $value.
    # Add $prob so that beta-prime follows the same interface as
    # the other mixing distributions in this function.
    gr <- do.call(
      I_betaprime_grad,
      c(
        list(
          z = y,
          X = X,
          be = be,
          phi = phi,
          d = 1L
        ),
        arg_vp
      )
    )

    gr$prob <- gr$value
    re <- gr$prob
  }

  # ------------------------------------------------------------
  # Gradient of negative log-likelihood
  # ------------------------------------------------------------
  
  re <- pmax(
    as.numeric(re),
    eps_prob
  )
  
  db <- gr$d_be
  
  if (length(dim(db)) == 3L) {
    db <- db[, , 1L]
  }
  
  dphi <- gr$d_phi
  
  if (is.matrix(dphi)) {
    dphi <- dphi[, 1L]
  }
  
  l_phi <-
    -sum(
      dphi /
      re
    )
  
  l_be <-
    colSums(
      -db /
      re
    )
  
  list(
    d_be = l_be,
    d_phi = l_phi
  )
}