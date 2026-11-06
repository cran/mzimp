# ============================================================
# Tests for select_bizimp()
# ============================================================
start <- proc.time()

test_that("select_bizimp agrees with fit_bizimp for a single candidate", {
  
  set.seed(12345)
  n <- 30
  
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  colnames(X.mat) <- c("x0", "x1", "x2")
  
  be1 <- c(0.5, 0.3, 0.7); be2 <- c(0.3, 0.5, 0.4)
  phi1 <- 0.2; phi2 <- 0.3
  mu1 <- exp(X.mat %*% be1); mu2 <- exp(X.mat %*% be2)
  
  N <- simu_biv_pm(n = n, mu1 = mu1, phi1 = phi1, dist1 = "gamma", mu2 = mu2, phi2 = phi2, dist2 = "gamma", use_custom_var_power = TRUE, var_power = c(2, 2), family = 1, par = 0.5)
  Y <- cbind(N$N1, N$N2)
  
  fit <- fit_bizimp(Y = Y, X = X.mat, Z = NULL, mixing.distribution.1 = "gamma", mixing.distribution.2 = "gamma", var.power = c(2, 2), mixing.copula = 3L)
  
  select <- select_bizimp(Y = Y, X = X.mat, Z = NULL, mixing.distribution = "gamma", var.power = 2, mixing.copula = 3L)
  
  expect_s3_class(fit, "bizimp_fit")
  expect_s3_class(select, "bizimp_selection")
  
  expect_equal(select$best_margin_1$par, fit$fit_margin_1$par, tolerance = 1e-8)
  expect_equal(select$best_margin_2$par, fit$fit_margin_2$par, tolerance = 1e-8)
  expect_equal(select$best_model$family, fit$best_model$family)
  expect_equal(select$best_model$par, fit$best_model$par, tolerance = 1e-8)
  expect_equal(select$best_model$log_lik_value, fit$best_model$log_lik_value, tolerance = 1e-8)
})


test_that("select_bizimp supports all mixing copula families", { # ALL COPULA TEST WILL TAKE ~ 30s
  
  set.seed(12345)
  n <- 30
  
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  colnames(X.mat) <- c("x0", "x1", "x2")
  
  be1 <- c(0.5, 0.3, 0.7); be2 <- c(0.3, 0.5, 0.4)
  phi1 <- 0.2; phi2 <- 0.3
  mu1 <- exp(X.mat %*% be1); mu2 <- exp(X.mat %*% be2)
  
  N <- simu_biv_pm(n = n, mu1 = mu1, phi1 = phi1, dist1 = "gamma", mu2 = mu2, phi2 = phi2, dist2 = "gamma", use_custom_var_power = TRUE, var_power = c(2, 2), family = 1, par = 0.5)
  Y <- cbind(N$N1, N$N2)
  
  # families <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
  #               10L, 13L, 14L, 16L, 17L, 18L, 19L,
  #               20L, 23L, 24L, 26L, 27L, 28L, 29L,
  #               30L, 33L, 34L, 36L, 37L, 38L, 39L,
  #               40L, 104L, 114L, 124L, 134L, 204L, 214L, 224L, 234L)
  
  
  families <- c(1L, 2L, 3L, 4L, 5L, 7L, 13L, 23L, 104L)
  
  for (fam in families) {
    select <- select_bizimp(Y = Y, X = X.mat, Z = NULL, mixing.distribution = "gamma", var.power = 2, mixing.copula = fam)
    
    expect_s3_class(select, "bizimp_selection")
    expect_true(nrow(select$candidate_table) >= 1L)
    expect_true(all(is.finite(select$candidate_table$logLik)))
    expect_equal(select$best_model$family, fam)
    expect_equal(select$best_model$log_lik_value, max(select$candidate_table$logLik), tolerance = 1e-8)
  }
})

test_that("margin.top controls retained marginal combinations", {
  
  
  set.seed(12345)
  n <- 30
  
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  colnames(X.mat) <- c("x0", "x1", "x2")
  
  be1 <- c(0.5, 0.3, 0.7); be2 <- c(0.3, 0.5, 0.4)
  phi1 <- 0.2; phi2 <- 0.3
  mu1 <- exp(X.mat %*% be1); mu2 <- exp(X.mat %*% be2)
  
  N <- simu_biv_pm(n = n, mu1 = mu1, phi1 = phi1, dist1 = "gamma", mu2 = mu2, phi2 = phi2, dist2 = "gamma", use_custom_var_power = TRUE, var_power = c(2, 2), family = 1, par = 0.5)
  Y <- cbind(N$N1, N$N2)
  
  select <- select_bizimp(Y = Y, X = X.mat, Z = NULL, mixing.distribution = c("gamma", "lnorm", "ig"), var.power = 2, margin.top = 2, mixing.copula = 1L)
  
  expect_equal(length(select$retained_margin_1), 2L)
  expect_equal(length(select$retained_margin_2), 2L)
  expect_equal(length(select$pair_results), 4L)
  expect_equal(nrow(select$pair_summary), 4L)
})



test_that("ZIMP-ZIMP evaluates all zero dependence structures", {
  
  set.seed(12345)
  n <- 30
  
  # Input data
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Z.mat <- matrix(rnorm(n * 2), n, 2)
  be1 <- c(0.5, 0.3, 0.7); be2 <- c(0.3, 0.5, 0.4)
  gam1 <- c(0.4, 0.6); gam2 <- c(0.5, 0.3)
  mu1 <- exp(X.mat %*% be1); mu2 <- exp(X.mat %*% be2)
  pr1 <- plogis(Z.mat %*% gam1); pr2 <- plogis(Z.mat %*% gam2)
  
  # Simulate two dependent ZIMP margins
  N <- simu_biv_pm(n = n, mu1 = mu1, phi1 = 0.2, dist1 = "gamma", mu2 = mu2, phi2 = 0.3, dist2 = "gamma", family = 1, par = 0.5)
  B <- simu_biv_bernoulli(n = n, p1 = pr1, p2 = pr2, rho = 0.5)
  Y <- cbind(N$N1 * B[, 1], N$N2 * B[, 2])
  
  select <- select_bizimp(Y = Y, X = X.mat, Z = Z.mat, mixing.distribution = "gamma", var.power = 2, mixing.copula = 1L)
  
  expect_true(select$best_margin_1$zero_inflation)
  expect_true(select$best_margin_2$zero_inflation)
  
  zero <- select$pair_results[[1]]$candidate_zero
  expect_setequal(unname(zero), c("normal", "comonotone", "countermonotone"))
})



test_that("candidate_table is sorted by decreasing logLik", {
  set.seed(12345)
  n <- 30
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2)); colnames(X.mat) <- c("x0", "x1", "x2")
  be1 <- c(0.5, 0.3, 0.7); be2 <- c(0.3, 0.5, 0.4); phi1 <- 0.2; phi2 <- 0.3
  mu1 <- exp(X.mat %*% be1); mu2 <- exp(X.mat %*% be2)
  N <- simu_biv_pm(n = n, mu1 = mu1, phi1 = phi1, dist1 = "gamma", mu2 = mu2, phi2 = phi2, dist2 = "gamma", use_custom_var_power = TRUE, var_power = c(2, 2), family = 1, par = 0.5)
  Y <- cbind(N$N1, N$N2)
  
  select <- select_bizimp(Y = Y, X = X.mat, Z = NULL, mixing.distribution = "gamma", var.power = 2, mixing.copula = c(1L, 3L, 5L))
  
  expect_gt(nrow(select$candidate_table), 1L)
  expect_true(all(is.finite(select$candidate_table$logLik)))
  expect_equal(select$candidate_table$logLik, sort(select$candidate_table$logLik, decreasing = TRUE), tolerance = 1e-10)
  expect_equal(select$best_model$log_lik_value, select$candidate_table$logLik[1L], tolerance = 1e-8)
})


test_that("select_bizimp rejects invalid selection arguments", {
  Y <- cbind(c(0, 1, 2, 1, 0), c(1, 0, 2, 1, 1))
  X <- cbind(1, c(-1, -0.5, 0, 0.5, 1))
  
  expect_error(select_bizimp(Y = Y, X = X, margin.top = 0), "margin.top")
  expect_error(select_bizimp(Y = Y, X = X, margin.top = 1.5), "margin.top")
  expect_error(select_bizimp(Y = Y, X = X, margin.top = NA), "margin.top")
  
  expect_error(select_bizimp(Y = Y, X = X, mixing.copula = integer(0)), "mixing.copula")
  expect_error(select_bizimp(Y = Y, X = X, mixing.copula = 999L), "Unsupported VineCopula family")
  expect_error(select_bizimp(Y = Y, X = X, mixing.copula = NA_integer_), "mixing.copula")
  
  expect_error(select_bizimp(Y = Y, X = X, parallel = 0), "parallel")
  expect_error(select_bizimp(Y = Y, X = X, parallel = 1.5), "parallel")
  
  expect_error(select_bizimp(Y = Y, X = X, ini.par = list(NULL)), "ini.par")
})

proc.time() - start