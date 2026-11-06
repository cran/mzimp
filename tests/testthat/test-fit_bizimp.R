# ============================================================
# Tests for fit_bizimp()
# ============================================================

test_that("fit_bizimp accepts margin-specific X and Z inputs", {
  
  set.seed(12345)
  n <- 30
  
  # Margin-specific design matrices
  X1 <- cbind(1, x1 = rnorm(n), x2 = rnorm(n))
  X2 <- cbind(1, x1 = rnorm(n))
  Z1 <- cbind(z1 = rnorm(n), z2 = rnorm(n))
  Z2 <- cbind(z1 = rnorm(n))
  
  # Parameters
  be1 <- c(0.5, 0.3, 0.7); be2 <- c(0.3, 0.5)
  gam1 <- c(0.4, 0.6); gam2 <- 0.5
  mu1 <- exp(X1 %*% be1); mu2 <- exp(X2 %*% be2)
  pr1 <- plogis(Z1 %*% gam1); pr2 <- plogis(Z2 %*% gam2)
  
  # Simulate two ZIMP margins
  N <- simu_biv_pm(n = n, mu1 = mu1, phi1 = 0.2, dist1 = "gamma", mu2 = mu2, phi2 = 0.3, dist2 = "gamma", family = 1, par = 0.5)
  B <- simu_biv_bernoulli(n = n, p1 = pr1, p2 = pr2, rho = 0.5)
  Y1 <- N$N1 * B[, 1]; Y2 <- N$N2 * B[, 2]
  
  # Margin-specific list input
  fit <- fit_bizimp(Y = list(Y1, Y2), X = list(X1, X2), Z = list(Z1, Z2), mixing.distribution.1 = "gamma", mixing.distribution.2 = "gamma", mixing.copula = 1, progress = FALSE)
  
  # Check that margin-specific designs are correctly used
  expect_s3_class(fit, "bizimp_fit")
  expect_equal(length(fit$fit_margin_1$be), ncol(X1))
  expect_equal(length(fit$fit_margin_2$be), ncol(X2))
  expect_equal(length(fit$fit_margin_1$gam), ncol(Z1))
  expect_equal(length(fit$fit_margin_2$gam), ncol(Z2))
  expect_true(is.finite(fit$best_model$log_lik_value))
})

test_that("fit_bizimp rejects invalid inputs", {
  
  set.seed(12345)
  n <- 30
  
  Y <- cbind(rpois(n, 2), rpois(n, 3))
  X <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Z <- matrix(rnorm(n * 2), n, 2)
  
  # Y must contain exactly two margins
  expect_error(fit_bizimp(Y = Y[, 1], X = X, Z = Z, progress = FALSE))
  expect_error(fit_bizimp(Y = cbind(Y, rpois(n, 2)), X = X, Z = Z, progress = FALSE))
  
  # Y must contain non-negative integer counts
  Y.neg <- Y; Y.neg[1, 1] <- -1
  expect_error(fit_bizimp(Y = Y.neg, X = X, Z = Z, progress = FALSE))
  
  Y.nonint <- Y; Y.nonint[1, 1] <- 1.5
  expect_error(fit_bizimp(Y = Y.nonint, X = X, Z = Z, progress = FALSE))
  
  # X and Z must have compatible numbers of observations
  expect_error(fit_bizimp(Y = Y, X = X[-1, ], Z = Z, progress = FALSE))
  expect_error(fit_bizimp(Y = Y, X = X, Z = Z[-1, ], progress = FALSE))
  
  # List inputs must contain exactly two margins
  expect_error(fit_bizimp(Y = list(Y[, 1]), X = list(X, X), Z = list(Z, Z), progress = FALSE))
  expect_error(fit_bizimp(Y = list(Y[, 1], Y[, 2]), X = list(X), Z = list(Z, Z), progress = FALSE))
  expect_error(fit_bizimp(Y = list(Y[, 1], Y[, 2]), X = list(X, X), Z = list(Z), progress = FALSE))
  
  # Invalid marginal distributions
  expect_error(fit_bizimp(Y = Y, X = X, Z = Z, mixing.distribution.1 = "invalid", progress = FALSE))
  expect_error(fit_bizimp(Y = Y, X = X, Z = Z, mixing.distribution.2 = "invalid", progress = FALSE))
  
  # Invalid variance-power specification
  expect_error(fit_bizimp(Y = Y, X = X, Z = Z, var.power = 4, progress = FALSE))
  expect_error(fit_bizimp(Y = Y, X = X, Z = Z, var.power = c(1, 2, 3), progress = FALSE))
  
  # Invalid mixing copula
  expect_error(fit_bizimp(Y = Y, X = X, Z = Z, mixing.copula = "gaussian", progress = FALSE))
  
  # Invalid zero-dependence specification
  expect_error(fit_bizimp(Y = Y, X = X, Z = Z, zero.dependence = "invalid", progress = FALSE))
  
  # Formula and direct interfaces cannot be used simultaneously
  dat <- data.frame(y1 = Y[, 1], y2 = Y[, 2], x1 = X[, 2], x2 = X[, 3], z1 = Z[, 1], z2 = Z[, 2])
  expect_error(fit_bizimp(y1 ~ 1 + x1 + x2 | z1 + z2, y2 ~ 1 + x1 + x2 | z1 + z2, data = dat, Y = Y, X = X, Z = Z, progress = FALSE))
})

test_that("fit_bizimp preserves independent fit_zimp marginal estimates under IFM", {
  
  set.seed(12345)
  n <- 30
  dist.vec <- c("lnorm", "gamma", "ig","weibull", "loglogistic", "betaprime")
  
  # Input data
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Z.mat <- matrix(rnorm(n * 2), n, 2)
  colnames(X.mat) <- c("x0", "x1", "x2")
  colnames(Z.mat) <- c("z1", "z2")
  
  # Input coefficients
  be1 <- c(0.5, 0.3, 0.7); be2 <- c(0.3, 0.5, 0.4)
  phi1 <- 0.2; phi2 <- 0.3
  gam1 <- c(0.4, 0.6); gam2 <- c(0.5, 0.3)
  mu1 <- exp(X.mat %*% be1); mu2 <- exp(X.mat %*% be2)
  pr1 <- exp(Z.mat %*% gam1) / (1 + exp(Z.mat %*% gam1))
  pr2 <- exp(Z.mat %*% gam2) / (1 + exp(Z.mat %*% gam2))
  dist.test <- sample(dist.vec, 2, replace = FALSE)
  
  
  # Simulate two ZIMP margins
  N <- simu_biv_pm(n = n, mu1 = mu1, phi1 = phi1, dist1 = dist.test[1],
                   mu2 = mu2, phi2 = phi2, , dist2 = dist.test[2], 
                   use_custom_var_power = TRUE,
                   var_power = c(1,3),
                   family = 1, par = 0.75)
  
  B <- simu_biv_bernoulli(n = n, p1 = pr1, p2 = pr2, rho = 0.5)
  Y1 <- N$N1* B[,1]
  Y2 <- N$N2* B[,2]
  Y <- cbind(Y1,Y2)
  
  # Fit margins independently
  fit1 <- fit_zimp(Y = Y1, X = X.mat, Z = Z.mat, mixing.distribution = dist.test[1], var.power = 1)
  fit2 <- fit_zimp(Y = Y2, X = X.mat, Z = Z.mat, mixing.distribution = dist.test[2], var.power = 3)
  
  # Fit bivariate model
  bifit <- fit_bizimp(Y = Y, X = X.mat, Z = Z.mat, mixing.distribution.1 = dist.test[1], mixing.distribution.2 = dist.test[2], var.power = c(1, 3), mixing.copula = 1, progress = FALSE)
  
  # Margin 1
  expect_equal(bifit$fit_margin_1$be, fit1$be, tolerance = 1e-10)
  expect_equal(bifit$fit_margin_1$phi, fit1$phi, tolerance = 1e-10)
  expect_equal(bifit$fit_margin_1$gam, fit1$gam, tolerance = 1e-10)
  expect_equal(bifit$fit_margin_1$par, fit1$par, tolerance = 1e-10)
  expect_equal(bifit$fit_margin_1$log_lik_value, fit1$log_lik_value, tolerance = 1e-10)
  
  # Margin 2
  expect_equal(bifit$fit_margin_2$be, fit2$be, tolerance = 1e-10)
  expect_equal(bifit$fit_margin_2$phi, fit2$phi, tolerance = 1e-10)
  expect_equal(bifit$fit_margin_2$gam, fit2$gam, tolerance = 1e-10)
  expect_equal(bifit$fit_margin_2$par, fit2$par, tolerance = 1e-10)
  expect_equal(bifit$fit_margin_2$log_lik_value, fit2$log_lik_value, tolerance = 1e-10)
})

test_that("fit_bizimp supports all four MP/ZIMP configurations", {
  
  set.seed(12345)
  n <- 30
  dist.vec <- c("lnorm", "gamma", "ig", "weibull", "loglogistic", "betaprime")
  
  # Input data
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Z.mat <- matrix(rnorm(n * 2), n, 2)
  colnames(X.mat) <- c("x0", "x1", "x2")
  colnames(Z.mat) <- c("z1", "z2")
  
  # Input coefficients
  be1 <- c(0.5, 0.3, 0.7); be2 <- c(0.3, 0.5, 0.4)
  phi1 <- 0.2; phi2 <- 0.3
  gam1 <- c(0.4, 0.6); gam2 <- c(0.5, 0.3)
  mu1 <- exp(X.mat %*% be1); mu2 <- exp(X.mat %*% be2)
  pr1 <- plogis(Z.mat %*% gam1); pr2 <- plogis(Z.mat %*% gam2)
  dist.test <- sample(dist.vec, 2, replace = FALSE)
  
  # Simulate dependent mixed-Poisson counts
  N <- simu_biv_pm(n = n, mu1 = mu1, phi1 = phi1, dist1 = dist.test[1], mu2 = mu2, phi2 = phi2, dist2 = dist.test[2], use_custom_var_power = TRUE, var_power = c(1, 3), family = 1, par = 0.75)
  
  # Simulate dependent Bernoulli activation indicators
  B <- simu_biv_bernoulli(n = n, p1 = pr1, p2 = pr2, rho = 0.5)
  
  # Four possible observed margins
  Y1.MP <- N$N1; Y2.MP <- N$N2
  Y1.ZIMP <- N$N1 * B[, 1]; Y2.ZIMP <- N$N2 * B[, 2]
  
  # 1. MP - MP
  fit.mm <- fit_bizimp(Y = cbind(Y1.MP, Y2.MP), X = X.mat, Z = NULL, 
                       mixing.distribution.1 = dist.test[1], mixing.distribution.2 = dist.test[2], var.power = c(1, 3), mixing.copula = 1, progress = FALSE)
  
  # 2. ZIMP - MP
  fit.zm <- fit_bizimp(Y = cbind(Y1.ZIMP, Y2.MP), X = X.mat, Z = list(Z.mat, NULL), 
                       mixing.distribution.1 = dist.test[1], mixing.distribution.2 = dist.test[2], var.power = c(1, 3), mixing.copula = 1, progress = FALSE)
  
  # 3. MP - ZIMP
  fit.mz <- fit_bizimp(Y = cbind(Y1.MP, Y2.ZIMP), X = X.mat, Z = list(NULL, Z.mat), 
                       mixing.distribution.1 = dist.test[1], mixing.distribution.2 = dist.test[2], var.power = c(1, 3), mixing.copula = 1, progress = FALSE)
  
  # 4. ZIMP - ZIMP
  fit.zz <- fit_bizimp(Y = cbind(Y1.ZIMP, Y2.ZIMP), X = X.mat, Z = Z.mat, 
                       mixing.distribution.1 = dist.test[1], mixing.distribution.2 = dist.test[2], var.power = c(1, 3), mixing.copula = 1, progress = FALSE)
  
  # Check model classes
  expect_s3_class(fit.mm, "bizimp_fit")
  expect_s3_class(fit.zm, "bizimp_fit")
  expect_s3_class(fit.mz, "bizimp_fit")
  expect_s3_class(fit.zz, "bizimp_fit")
  
  # Check marginal structures
  expect_equal(fit.mm$zero_inflation, c(FALSE, FALSE))
  expect_equal(fit.zm$zero_inflation, c(TRUE, FALSE))
  expect_equal(fit.mz$zero_inflation, c(FALSE, TRUE))
  expect_equal(fit.zz$zero_inflation, c(TRUE, TRUE))
  
  # Check finite likelihoods
  expect_true(is.finite(fit.mm$best_model$log_lik_value))
  expect_true(is.finite(fit.zm$best_model$log_lik_value))
  expect_true(is.finite(fit.mz$best_model$log_lik_value))
  expect_true(is.finite(fit.zz$best_model$log_lik_value))
  
  # Zero dependence is estimated only for ZIMP-ZIMP
  expect_false(fit.mm$zero_dependence_estimated)
  expect_false(fit.zm$zero_dependence_estimated)
  expect_false(fit.mz$zero_dependence_estimated)
  expect_true(fit.zz$zero_dependence_estimated)
})

test_that("fit_bizimp supports all mixing copula families", {
  
  set.seed(12345)
  n <- 30
  
  # Input data
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  colnames(X.mat) <- c("x0", "x1", "x2")
  
  # Input coefficients
  be1 <- c(0.5, 0.3, 0.7); be2 <- c(0.3, 0.5, 0.4)
  phi1 <- 0.2; phi2 <- 0.3
  mu1 <- exp(X.mat %*% be1); mu2 <- exp(X.mat %*% be2)
  
  # Simulate two dependent MP margins
  N <- simu_biv_pm(n = n, mu1 = mu1, phi1 = phi1, dist1 = "gamma", mu2 = mu2, phi2 = phi2, dist2 = "gamma", use_custom_var_power = TRUE, var_power = c(2, 2), family = 1, par = 0.5)
  Y <- cbind(N$N1, N$N2)
  
  # All supported VineCopula families
  families <- c(1L, 2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L,
                10L, 13L, 14L, 16L, 17L, 18L, 19L,
                20L, 23L, 24L, 26L, 27L, 28L, 29L,
                30L, 33L, 34L, 36L, 37L, 38L, 39L,
                40L, 104L, 114L, 124L, 134L, 204L, 214L, 224L, 234L)
  
  
  for (fam in families) {
    
    fit <- fit_bizimp(Y = Y, X = X.mat, Z = NULL, mixing.distribution.1 = "gamma", mixing.distribution.2 = "gamma", var.power = c(2, 2), mixing.copula = fam, progress = FALSE)
    
    # Model successfully returned
    expect_s3_class(fit, "bizimp_fit")
    
    # Requested family was actually fitted
    expect_equal(as.integer(fit$best_model$family), fam)
    
    # Valid likelihood
    expect_true(is.finite(fit$best_model$log_lik_value))
    
    # Information criteria are available
    expect_true(is.finite(fit$AIC))
    expect_true(is.finite(fit$BIC))
  }
})

test_that("fit_bizimp supports mixing copulas 1-5 with all zero-dependence structures", {
  
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
  
  mixing.families <- 1:5
  zero.families <- c("gaussian", "comonotone", "countermonotone")
  
  for (fam in mixing.families) {
    for (zero.fam in zero.families) {
      
      fit <- fit_bizimp(Y = Y, X = X.mat, Z = Z.mat, mixing.distribution.1 = "gamma", mixing.distribution.2 = "gamma", var.power = c(2, 2), mixing.copula = fam, zero.dependence = zero.fam, progress = FALSE)
      
      expect_s3_class(fit, "bizimp_fit")
      expect_equal(as.integer(fit$best_model$family), fam)
      expect_true(is.finite(fit$best_model$log_lik_value))
      expect_true(is.finite(fit$AIC))
      expect_true(is.finite(fit$BIC))
    }
  }
})