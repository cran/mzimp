#============================================================================
# Test select_multi_margin
#============================================================================

test_that("select_multi_margin agrees with select_bizimp for two margins", {
  set.seed(12345)
  n <- 30
  
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  colnames(X.mat) <- c("x0", "x1", "x2")
  
  be1 <- c(0.5, 0.3, 0.7); be2 <- c(0.3, 0.5, 0.4)
  mu1 <- exp(X.mat %*% be1); mu2 <- exp(X.mat %*% be2)
  
  N <- simu_biv_pm(n = n, mu1 = mu1, phi1 = 0.2, dist1 = "gamma", mu2 = mu2, phi2 = 0.3, dist2 = "gamma", use_custom_var_power = TRUE, var_power = c(2, 2), family = 1, par = 0.5)
  Y <- cbind(N$N1, N$N2)
  
  multi <- select_multi_margin(Y = Y, X = X.mat, Z = NULL, mixing.distribution = c("gamma", "lnorm", "ig"), var.power = 2, parallel = 1)
  
  bi <- select_bizimp(Y = Y, X = X.mat, Z = NULL, mixing.distribution = c("gamma", "lnorm", "ig"), var.power = 2, margin.top = 1, mixing.copula = 1L)
  
  expect_s3_class(multi, "multi_margin_fit")
  expect_s3_class(bi, "bizimp_selection")
  
  expect_equal(multi$selected_margins[[1]]$distribution, bi$best_margin_1$distribution)
  expect_equal(multi$selected_margins[[2]]$distribution, bi$best_margin_2$distribution)
  
  expect_equal(multi$selected_margins[[1]]$var.power, bi$best_margin_1$var.power)
  expect_equal(multi$selected_margins[[2]]$var.power, bi$best_margin_2$var.power)
  
  expect_equal(multi$selected_margins[[1]]$zero_inflation, bi$best_margin_1$zero_inflation)
  expect_equal(multi$selected_margins[[2]]$zero_inflation, bi$best_margin_2$zero_inflation)
  
  expect_equal(multi$selected_margins[[1]]$be, bi$best_margin_1$be, tolerance = 1e-8)
  expect_equal(multi$selected_margins[[2]]$be, bi$best_margin_2$be, tolerance = 1e-8)
  
  expect_equal(multi$selected_margins[[1]]$phi, bi$best_margin_1$phi, tolerance = 1e-8)
  expect_equal(multi$selected_margins[[2]]$phi, bi$best_margin_2$phi, tolerance = 1e-8)
  
  expect_equal(multi$selected_margins[[1]]$log_lik_value, bi$best_margin_1$log_lik_value, tolerance = 1e-8)
  expect_equal(multi$selected_margins[[2]]$log_lik_value, bi$best_margin_2$log_lik_value, tolerance = 1e-8)
})


test_that("select_multi_margin supports all distributions and var.power values and is equivalent in serial and parallel", {
  set.seed(12345)
 
  n <- 30 ; k <- 3; J <- 1
  psi.N <- simu_psi_N(k,J)
  
  X <- cbind(1,matrix(rnorm(2*n),n,2))
  
  beta <- c(0.5,0.3,0.7)
  
  dist.vec <- c("gamma","lnorm", "ig", "loglogistic", "betaprime", "weibull")
  powers <- c(1L,2L,3L)
  
  dist.simu <- sample(dist.vec,3)
  
  sim <- simu_multizimp(X = X, beta = beta, Z = NULL, 
                        mixing.distribution = dist.simu,
                        psi_N = psi.N)
 
  serial <- select_multi_margin(Y = sim$Y, X = sim$X, Z = NULL, mixing.distribution = dist.vec , var.power = powers, compute.I = TRUE, parallel = 1)
 
  parallel <- select_multi_margin(Y = sim$Y, X = sim$X, Z = NULL, mixing.distribution = dist.vec , var.power = powers, parallel = 4)
  
  expect_s3_class(serial, "multi_margin_fit")
  expect_s3_class(parallel, "multi_margin_fit")
  
  expect_equal(nrow(serial$candidate_grid), 18L)
  expect_equal(nrow(parallel$candidate_grid), 18L)
  
  expect_setequal(serial$candidate_grid$distribution, dist.vec)
  expect_setequal(serial$candidate_grid$var.power, powers)
  
  expect_equal(length(serial$selected_margins), 3L)
  expect_equal(length(parallel$selected_margins), 3L)
  
  expect_equal(serial$selected_summary$distribution, parallel$selected_summary$distribution)
  expect_equal(serial$selected_summary$var.power, parallel$selected_summary$var.power)
  expect_equal(serial$selected_summary$zero_inflation, parallel$selected_summary$zero_inflation)
  expect_equal(serial$selected_summary$logLik, parallel$selected_summary$logLik, tolerance = 1e-8)
  expect_equal(serial$selected_summary$AIC, parallel$selected_summary$AIC, tolerance = 1e-8)
  expect_equal(serial$selected_summary$BIC, parallel$selected_summary$BIC, tolerance = 1e-8)
  
  ## TEST I matrices
  expect_true(serial$I_precomputed)
  expect_equal(serial$I_d0, 10L)
  expect_equal(length(serial$I_y_list), k)
  expect_equal(length(serial$I_0_list), k)
  expect_equal(length(serial$pb_list), k)
  expect_equal(length(serial$y0_list), k)
  
  for (j in seq_len(k)) {
    expect_true(is.matrix(serial$I_y_list[[j]]))
    expect_true(is.matrix(serial$I_0_list[[j]]))
    expect_equal(nrow(serial$I_y_list[[j]]), n)
    expect_equal(nrow(serial$I_0_list[[j]]), n)
    expect_equal(ncol(serial$I_y_list[[j]]), 10L)
    expect_equal(ncol(serial$I_0_list[[j]]), 10L)
    expect_true(all(is.finite(serial$I_y_list[[j]])))
    expect_true(all(is.finite(serial$I_0_list[[j]])))
    expect_equal(serial$pb_list[[j]], rep(1, n))
    expect_null(serial$y0_list[[j]])
  }
  
  for (j in 1:k) {
    expect_equal(serial$selected_margins[[j]]$distribution, parallel$selected_margins[[j]]$distribution)
    expect_equal(serial$selected_margins[[j]]$var.power, parallel$selected_margins[[j]]$var.power)
    expect_equal(serial$selected_margins[[j]]$be, parallel$selected_margins[[j]]$be, tolerance = 1e-8)
    expect_equal(serial$selected_margins[[j]]$phi, parallel$selected_margins[[j]]$phi, tolerance = 1e-8)
    expect_equal(serial$selected_margins[[j]]$log_lik_value, parallel$selected_margins[[j]]$log_lik_value, tolerance = 1e-8)
  }
})


test_that("select_multi_margin supports all distributions and var.power values with zero inflation", {
  set.seed(12345)
  
  n <- 30; k <- 3; J <- 1
  psi.N <- simu_psi_N(k, J)
  psi.B <- simu_psi_B(k)
  
  X.mat <- cbind(1, matrix(rnorm(2*n), n, 2))
  Z.mat <- matrix(rnorm(2*n), n, 2)
  beta <- c(0.5, 0.3, 0.7)
  gam <- c(0.4, 0.6)
  
  dist.vec <- c("gamma", "lnorm", "ig", "loglogistic", "betaprime", "weibull")
  powers <- c(1L, 2L, 3L)
  
  dist.simu <- sample(dist.vec, k)
  
  sim <- simu_multizimp(X = X.mat, beta = beta, zero_inflation = rep(TRUE, k), Z = Z.mat, gamma = gam, mixing.distribution = dist.simu, psi_N = psi.N, psi_B = psi.B)
  
  fit <- select_multi_margin(Y = sim$Y, X = X.mat, Z = Z.mat, mixing.distribution = dist.vec, var.power = powers, parallel = 1, compute.I=TRUE)
  
  expect_s3_class(fit, "multi_margin_fit")
  expect_equal(length(fit$selected_margins), k)
  expect_true(all(vapply(fit$selected_margins, function(x) !is.null(x), logical(1))))
  
  expected.grid <- expand.grid(distribution = dist.vec, var.power = powers, stringsAsFactors = FALSE)
  actual.grid <- unique(fit$candidate_grid[, c("distribution", "var.power")])
  
  expect_equal(nrow(actual.grid), 18L)
  expect_setequal(paste(actual.grid$distribution, actual.grid$var.power), paste(expected.grid$distribution, expected.grid$var.power))
  
  expect_true(all(fit$selected_summary$zero_inflation))
  expect_true(all(vapply(fit$selected_margins, function(x) isTRUE(x$zero_inflation), logical(1))))
  
  expect_true(all(fit$selected_summary$distribution %in% dist.vec))
  expect_true(all(fit$selected_summary$var.power %in% powers))
  expect_true(all(is.finite(fit$selected_summary$logLik)))
  expect_true(all(is.finite(fit$selected_summary$AIC)))
  expect_true(all(is.finite(fit$selected_summary$BIC)))
  
  for (j in seq_len(k)) {
    expect_equal(length(fit$selected_margins[[j]]$be), ncol(X.mat))
    expect_equal(length(fit$selected_margins[[j]]$gam), ncol(Z.mat))
    expect_true(is.finite(fit$selected_margins[[j]]$phi))
    expect_true(fit$selected_margins[[j]]$phi > 0)
    expect_true(all(is.finite(fit$selected_margins[[j]]$be)))
    expect_true(all(is.finite(fit$selected_margins[[j]]$gam)))
    expect_true(is.finite(fit$selected_margins[[j]]$log_lik_value))
  }
  
  ## TEST I matrices
  expect_true(fit$I_precomputed)
  expect_equal(fit$I_d0, 10L)
  expect_equal(length(fit$I_y_list), k)
  expect_equal(length(fit$I_0_list), k)
  expect_equal(length(fit$pb_list), k)
  expect_equal(length(fit$y0_list), k)
  
  for (j in seq_len(k)) {
    expect_true(is.matrix(fit$I_y_list[[j]]))
    expect_true(is.matrix(fit$I_0_list[[j]]))
    expect_equal(nrow(fit$I_y_list[[j]]), n)
    expect_equal(nrow(fit$I_0_list[[j]]), n)
    expect_equal(ncol(fit$I_y_list[[j]]), 10L)
    expect_equal(ncol(fit$I_0_list[[j]]), 10L)
    
    expect_true(all(is.finite(fit$I_y_list[[j]])))
    expect_true(all(is.finite(fit$I_0_list[[j]])))
    
    pb.expected <- plogis(as.vector(Z.mat %*% fit$selected_margins[[j]]$gam))
    expect_equal(fit$pb_list[[j]], pb.expected, tolerance=1e-10)
    expect_equal(fit$y0_list[[j]], sim$Y[, j] == 0L)
  }
  
})


test_that("select_multi_margin check.zero compares MP and intercept-only ZIMP independently for each margin", {
  set.seed(12345)
  
  n <- 50; k <- 3
  X <- cbind(1, matrix(rnorm(2*n), n, 2))
  beta <- c(0.5, 0.3, 0.7)
  
  sim <- simu_multizimp(X=X, beta=beta, Z=NULL, mixing.distribution=rep("gamma", k), psi_N=simu_psi_N(k, 1))
  
  fit <- select_multi_margin(Y=sim$Y, X=X, Z=NULL, check.zero=TRUE, mixing.distribution=c("gamma", "lnorm"), var.power=c(1L, 2L), parallel=1)
  
  expect_s3_class(fit, "multi_margin_fit")
  
  # 2 distributions x 2 var.power x 2 zero-inflation states
  expect_equal(nrow(fit$candidate_grid), 8L)
  expect_setequal(fit$candidate_grid$distribution, c("gamma", "lnorm"))
  expect_setequal(fit$candidate_grid$var.power, c(1L, 2L))
  expect_setequal(fit$candidate_grid$zero_inflation, c(FALSE, TRUE))
  
  expected <- expand.grid(distribution=c("gamma", "lnorm"), var.power=c(1L, 2L), zero_inflation=c(FALSE, TRUE), stringsAsFactors=FALSE)
  actual <- fit$candidate_grid[, c("distribution", "var.power", "zero_inflation")]
  
  expect_setequal(paste(actual$distribution, actual$var.power, actual$zero_inflation), paste(expected$distribution, expected$var.power, expected$zero_inflation))
  
  expect_equal(length(fit$selected_margins), k)
  expect_equal(nrow(fit$selected_summary), k)
  
  for (j in seq_len(k)) {
    mj <- fit$selected_margins[[j]]
    
    expect_true(mj$zero_inflation %in% c(FALSE, TRUE))
    
    if (isTRUE(mj$zero_inflation)) {
      expect_equal(length(mj$gam), 1L)
      expect_true(is.finite(mj$gam))
    } else {
      expect_null(mj$gam)
    }
  }
})


test_that("select_multi_margin supports margin-specific X and Z designs", {
  set.seed(12345)
  
  n <- 50; k <- 3
  
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  z1 <- rnorm(n); z2 <- rnorm(n); z3 <- rnorm(n)
  
  X <- list(cbind(1, x1), cbind(1, x1, x2), cbind(1, x1, x2, x3))
  Z <- list(matrix(1, n, 1), cbind(1, z1), cbind(1, z1, z2))
  
  beta <- list(c(0.5, 0.3), c(0.4, 0.2, -0.1), c(0.3, 0.2, -0.1, 0.15))
  gamma <- list(0.5, c(0.4, 0.3), c(0.3, 0.2, -0.2))
  
  Y <- matrix(0L, n, k)
  
  for (j in seq_len(k)) {
    mu <- exp(as.vector(X[[j]] %*% beta[[j]]))
    pb <- plogis(as.vector(Z[[j]] %*% gamma[[j]]))
    B <- rbinom(n, 1, pb)
    Y[, j] <- B * rpois(n, mu)
  }
  
  fit <- select_multi_margin(Y=Y, X=X, Z=Z, mixing.distribution="gamma", var.power=2, parallel=1)
  
  expect_s3_class(fit, "multi_margin_fit")
  expect_equal(length(fit$selected_margins), k)
  expect_true(fit$X_is_list)
  expect_true(fit$Z_is_list)
  
  for (j in seq_len(k)) {
    expect_true(fit$selected_margins[[j]]$zero_inflation)
    expect_equal(length(fit$selected_margins[[j]]$be), ncol(X[[j]]))
    expect_equal(length(fit$selected_margins[[j]]$gam), ncol(Z[[j]]))
    expect_true(all(is.finite(fit$selected_margins[[j]]$be)))
    expect_true(all(is.finite(fit$selected_margins[[j]]$gam)))
    expect_true(is.finite(fit$selected_margins[[j]]$phi))
    expect_true(is.finite(fit$selected_margins[[j]]$log_lik_value))
  }
})

test_that("select_multi_margin returns a stable public output contract and is silent", {
  set.seed(12345)
  
  n <- 30; k <- 3
  Y <- matrix(rpois(n*k, 2), n, k)
  X <- cbind(1, matrix(rnorm(2*n), n, 2))
  
  expect_silent(fit <- select_multi_margin(Y=Y, X=X, Z=NULL, mixing.distribution="gamma", var.power=2, parallel=1))
  
  expect_s3_class(fit, "multi_margin_fit")
  
  expected.names <- c(
    "Y", "X", "Z",
    "X_is_list", "Z_is_list", "X_columns", "Z_columns",
    "selected_margins",
    "benchmark.distribution", "benchmark_by_margin", "benchmark",
    "selected_summary",
    "I_precomputed", "I_d0", "I_y_list", "I_0_list", "pb_list", "y0_list",
    "mixing.distribution", "var.power", "candidate_grid",
    "check.zero", "selection.criterion",
    "n", "k", "response_names",
    "parallel", "n_workers", "performance",
    "summary"
  )
  
  expect_equal(names(fit), expected.names)
  
  expect_equal(fit$n, n)
  expect_equal(fit$k, k)
  expect_equal(dim(fit$Y), c(n, k))
  expect_equal(length(fit$selected_margins), k)
  expect_equal(nrow(fit$selected_summary), k)
  
  expect_false(fit$I_precomputed)
  expect_null(fit$I_d0)
  expect_null(fit$I_y_list)
  expect_null(fit$I_0_list)
  expect_null(fit$pb_list)
  expect_null(fit$y0_list)
  
  expect_type(fit$summary, "list")
  expect_equal(fit$summary$selected_summary, fit$selected_summary)
  expect_equal(fit$summary$benchmark, fit$benchmark)
  expect_equal(fit$summary$benchmark.distribution, fit$benchmark.distribution)
  expect_equal(fit$summary$selection.criterion, fit$selection.criterion)
  expect_equal(fit$summary$n, fit$n)
  expect_equal(fit$summary$k, fit$k)
  expect_equal(fit$summary$response_names, fit$response_names)
  expect_equal(fit$summary$performance, fit$performance)
})


test_that("select_multi_margin validates its inputs", {
  set.seed(12345)
  
  n <- 20; k <- 3
  Y <- matrix(rpois(n*k, 2), n, k)
  X <- cbind(1, matrix(rnorm(2*n), n, 2))
  Z <- cbind(1, rnorm(n))
  
  # Y must contain non-negative integer counts
  Y.neg <- Y; Y.neg[1,1] <- -1
  expect_error(select_multi_margin(Y=Y.neg, X=X, mixing.distribution="gamma", var.power=2), "nonnegative")
  
  Y.nonint <- Y; Y.nonint[1,1] <- 1.5
  expect_error(select_multi_margin(Y=Y.nonint, X=X, mixing.distribution="gamma", var.power=2), "integer")
  
  Y.na <- Y; Y.na[1,1] <- NA
  expect_error(select_multi_margin(Y=Y.na, X=X, mixing.distribution="gamma", var.power=2))
  
  # X and Z must have the same number of rows as Y
  expect_error(select_multi_margin(Y=Y, X=X[-1,,drop=FALSE], mixing.distribution="gamma", var.power=2), "nrow")
  expect_error(select_multi_margin(Y=Y, X=X, Z=Z[-1,,drop=FALSE], mixing.distribution="gamma", var.power=2), "nrow")
  
  # Margin-specific X/Z lists must have length k
  expect_error(select_multi_margin(Y=Y, X=list(X,X), mixing.distribution="gamma", var.power=2), "length")
  expect_error(select_multi_margin(Y=Y, X=X, Z=list(Z,Z), mixing.distribution="gamma", var.power=2), "length")
  
  
  # var.power only accepts 1, 2, or 3
  expect_error(select_multi_margin(Y=Y, X=X, mixing.distribution="gamma", var.power=4), "var.power")
  expect_error(select_multi_margin(Y=Y, X=X, mixing.distribution="gamma", var.power=c(1,2,4)), "var.power")
  
  # Unsupported mixing distribution
  expect_error(select_multi_margin(Y=Y, X=X, mixing.distribution="invalid_distribution", var.power=2))
  
  # Invalid zero-check specification
  expect_error(select_multi_margin(Y=Y, X=X, check.zero=NA, mixing.distribution="gamma", var.power=2), "check.zero")
  expect_error(select_multi_margin(Y=Y, X=X, check.zero=1, mixing.distribution="gamma", var.power=2), "check.zero")
  
  # Invalid checkerboard resolution
  expect_error(select_multi_margin(Y=Y, X=X, mixing.distribution="gamma", var.power=2, compute.I=TRUE, I.d0=1), "I.d0")
  
  # Invalid parallel specification
  expect_error(select_multi_margin(Y=Y, X=X, mixing.distribution="gamma", var.power=2, parallel=0))
  expect_error(select_multi_margin(Y=Y, X=X, mixing.distribution="gamma", var.power=2, parallel=-1))
  
  # check.zero is ignored when Z is explicitly supplied
  expect_warning(
    fit <- select_multi_margin(Y=Y, X=X, Z=Z, check.zero=TRUE, mixing.distribution="gamma", var.power=2, parallel=1),
    "check.zero is ignored"
  )
  expect_true(all(fit$candidate_grid$zero_inflation))
})

test_that("select_multi_margin supports mixed ZIMP and MP margins through Z list", {
  set.seed(12345)
  n <- 30; k <- 3
  X <- cbind(1, rnorm(n)); Z1 <- cbind(1, rnorm(n)); Z2 <- cbind(1, rnorm(n))
  beta <- c(0.5, 0.3); gamma1 <- c(0.5, 0.3); gamma2 <- c(0.4, -0.2)
  mu <- exp(as.vector(X %*% beta)); pb1 <- plogis(as.vector(Z1 %*% gamma1)); pb2 <- plogis(as.vector(Z2 %*% gamma2))
  Y <- cbind(rbinom(n, 1, pb1)*rpois(n, mu), rbinom(n, 1, pb2)*rpois(n, mu), rpois(n, mu))
  
  fit <- select_multi_margin(Y=Y, X=X, Z=list(Z1,Z2,NULL), mixing.distribution="gamma", var.power=2, compute.I=TRUE, I.d0=3, parallel=1)
  
  expect_s3_class(fit, "multi_margin_fit"); expect_equal(fit$n, n); expect_equal(fit$k, k)
  expect_true(fit$selected_margins[[1]]$zero_inflation); expect_true(fit$selected_margins[[2]]$zero_inflation); expect_false(fit$selected_margins[[3]]$zero_inflation)
  expect_equal(length(fit$selected_margins[[1]]$gam), ncol(Z1)); expect_equal(length(fit$selected_margins[[2]]$gam), ncol(Z2)); expect_null(fit$selected_margins[[3]]$gam)
  expect_true(all(is.finite(fit$selected_margins[[1]]$gam))); expect_true(all(is.finite(fit$selected_margins[[2]]$gam)))
  expect_true(all(fit$pb_list[[1]] > 0 & fit$pb_list[[1]] < 1)); expect_true(all(fit$pb_list[[2]] > 0 & fit$pb_list[[2]] < 1)); expect_equal(fit$pb_list[[3]], rep(1, n))
  expect_equal(length(fit$y0_list[[1]]), n); expect_equal(length(fit$y0_list[[2]]), n); expect_null(fit$y0_list[[3]])
})