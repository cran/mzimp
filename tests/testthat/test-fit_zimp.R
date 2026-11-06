# ============================================================
# Tests for fit_zimp()
# ============================================================

n <- 100
test_that("fit_zimp supports all mixing distributions and variance powers for Mixed Poisson models", {
  
  set.seed(12345)
  
  # Input data
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  colnames(X.mat) <- c("x0", "x1", "x2")
  
  # Input coefficients
  be <- c(0.5, 0.3, 0.7)
  phi <- 0.2
  mu.vec <- exp(X.mat %*% be)
  
  # Mixing distributions and variance powers
  dist.vec <- c("lnorm", "gamma", "ig", "weibull", "loglogistic", "betaprime")
  var.power.vec <- c(1, 2, 3)
  
  # Simulation by lognormal
  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = "lnorm")$Y
  
  # Test all mixing distribution × variance power combinations
  for (dist.check in dist.vec) {
    for (p.check in var.power.vec) {
      
      fit <- fit_zimp(Y = Y, X = X.mat, mixing.distribution = dist.check, var.power = p.check)
      
      # Test log-likelihood and parameter dimensions
      expect_true(is.finite(fit$log_lik_value))
      expect_length(fit$be, ncol(X.mat))
      expect_length(fit$phi, 1)
      
      # Test estimated parameters
      expect_true(all(is.finite(fit$be)))
      expect_true(all(is.finite(fit$phi)))
      expect_true(all(is.finite(fit$hessian)))
      expect_true(fit$phi > 0)
      
      # Test variance power
      expect_equal(fit$var.power, p.check)
    }
  }
})

test_that("fit_zimp supports all mixing distributions and variance powers for Zero-Inflated Mixed Poisson models", {
  
  set.seed(12345)
  
  # Input data
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  Z.mat <- matrix(rnorm(n * 2), n, 2)
  colnames(X.mat) <- c("x0", "x1", "x2")
  colnames(Z.mat) <- c("z1", "z2")
  
  # Input coefficients
  be <- c(0.5, 0.3, 0.7)
  phi <- 0.2
  gam <- c(0.4, 0.6)
  mu.vec <- exp(X.mat %*% be)
  pr.vec <- exp(Z.mat %*% gam) / (1 + exp(Z.mat %*% gam))
  
  # Mixing distributions and variance powers
  dist.vec <- c("lnorm", "gamma", "ig", "weibull", "loglogistic", "betaprime")
  var.power.vec <- c(1, 2, 3)
  
  # Simulation by lognormal
  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = "lnorm",
                 zero_inflation = TRUE, prob = pr.vec)$Y
  
  # Test all mixing distribution × variance power combinations
  for (dist.check in dist.vec) {
    for (p.check in var.power.vec) {
      
      fit <- fit_zimp(Y = Y, X = X.mat, Z = Z.mat,
                      mixing.distribution = dist.check, var.power = p.check)
      
      # Test log-likelihood and parameter dimensions
      expect_true(is.finite(fit$log_lik_value))
      expect_length(fit$be, ncol(X.mat))
      expect_length(fit$gam, ncol(Z.mat))
      expect_length(fit$phi, 1)
      
      # Test estimated parameters
      expect_true(all(is.finite(fit$be)))
      expect_true(all(is.finite(fit$gam)))
      expect_true(all(is.finite(fit$phi)))
      expect_true(all(is.finite(fit$hessian)))
      expect_true(fit$phi > 0)
      
      # Test variance power
      expect_equal(fit$var.power, p.check)
    }
  }
})

test_that("fit_zimp formula and matrix interfaces are equivalent", {

  set.seed(12345)
  X.mat <- matrix(rnorm(n*2),n,2)
  Z.mat <- matrix(rnorm(n*2),n,2)
  colnames(X.mat) <- c("x1","x2")
  colnames(Z.mat) <- c("z1","z2")
  
  # Test.1 Mixed Poisson #
  
  ## Test 1.a. WITHOUT Intercept -------------------------------------------------#
  
  ### Simulation MP by lnorm
  be <- c(0.3, 0.7)
  phi <- 0.2
  mu.vec <- exp(X.mat %*% be)
  dist.simu <- "gamma"
  
  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = dist.simu)$Y
  ### Fit using matrix interface
  fit.matrix <- fit_zimp(
    Y = Y, X = X.mat, mixing.distribution = "lnorm"
  )
  ### Fit using formula interface
  dat <- data.frame(y = Y, x1 = X.mat[,1], x2 = X.mat[,2])
  fit.formula <- fit_zimp(
    y ~ x1 + x2,
    data = dat,
    mixing.distribution = "lnorm"
  )
  
  ### Compare results
  expect_equal(unname(fit.formula$be), unname(fit.matrix$be),tolerance = 1e-8)
  expect_equal(fit.formula$phi,fit.matrix$phi,tolerance = 1e-8)
  expect_equal(fit.formula$log_lik_value,fit.matrix$log_lik_value,tolerance = 1e-8)
  expect_equal(fit.formula$AIC,fit.matrix$AIC,tolerance = 1e-8)
  expect_equal(fit.formula$BIC,fit.matrix$BIC,tolerance = 1e-8)
  
  #-------------------------------------------------------------------------------#
  
  
  ## Test 1.b. WITH Intercept ----------------------------------------------------#
  
  ### Simulation ZIMP by lnorm
  be <- c(0.5, 0.3, 0.7)
  phi <- 0.2
  mu.vec <- exp(cbind(1,X.mat) %*% be) # add col 1 to X.mat
  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = dist.simu)$Y
  ### Fit using matrix interface
  fit.matrix <- fit_zimp(
    Y = Y, X = cbind(1,X.mat), mixing.distribution = "lnorm"
  )
  ### Fit using formula interface
  dat <- data.frame(y = Y, x1 = X.mat[,1], x2 = X.mat[,2])
  fit.formula <- fit_zimp(
    y ~ 1 + x1 + x2,
    data = dat,
    mixing.distribution = "lnorm"
  )
  
  ### Compare results
  expect_equal(unname(fit.formula$be), unname(fit.matrix$be),tolerance = 1e-8)
  expect_equal(fit.formula$phi,fit.matrix$phi,tolerance = 1e-8)
  expect_equal(fit.formula$log_lik_value,fit.matrix$log_lik_value,tolerance = 1e-8)
  expect_equal(fit.formula$AIC,fit.matrix$AIC,tolerance = 1e-8)
  expect_equal(fit.formula$BIC,fit.matrix$BIC,tolerance = 1e-8)
  
  
  # Test.2 Zero Inflated Mixed Poisson #
  
  ## Test 2.a. WITHOUT Intercept + WITHOUT Intercept ---------------------------#
  
  ### Simulation ZIMP by lnorm
  be <- c(0.3, 0.7)
  phi <- 0.2
  mu.vec <- exp(X.mat %*% be)
  gam <- c(0.4,0.6)
  pr.vec <- exp(Z.mat %*% gam)/(1+exp(Z.mat %*% gam))

  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = dist.simu,
                 zero_inflation = TRUE, prob = pr.vec)$Y
  
  ### Fit using matrix interface
  fit.matrix <- fit_zimp(
    Y = Y, X = X.mat, Z = Z.mat, mixing.distribution = "lnorm"
  )
  ### Fit using formula interface
  dat <- data.frame(y = Y, x1 = X.mat[,1], x2 = X.mat[,2], z1 = Z.mat[,1], z2 = Z.mat[,2])
  fit.formula <- fit_zimp(
    y ~ x1 + x2|z1 + z2,
    data = dat,
    mixing.distribution = "lnorm"
  )
  
  ### Compare results
  expect_equal(unname(fit.formula$be), unname(fit.matrix$be),tolerance = 1e-8)
  expect_equal(fit.formula$phi,fit.matrix$phi,tolerance = 1e-8)
  expect_equal(unname(fit.formula$gam), unname(fit.matrix$gam),tolerance = 1e-8)
  expect_equal(fit.formula$log_lik_value,fit.matrix$log_lik_value,tolerance = 1e-8)
  expect_equal(fit.formula$AIC,fit.matrix$AIC,tolerance = 1e-8)
  expect_equal(fit.formula$BIC,fit.matrix$BIC,tolerance = 1e-8)
  
  #-------------------------------------------------------------------------------#
  
  ## Test 2.b. WITH Intercept + WITHOUT Intercept ---------------------------#
  
  ### Simulation ZIMP by lnorm
  be <- c(0.5, 0.3, 0.7)
  phi <- 0.2
  mu.vec <- exp(cbind(1,X.mat) %*% be)
  gam <- c(0.4,0.6)
  pr.vec <- exp(Z.mat %*% gam)/(1+exp(Z.mat %*% gam))
  
  
  
  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = dist.simu,
                 zero_inflation = TRUE, prob = pr.vec)$Y
  
  ### Fit using matrix interface
  fit.matrix <- fit_zimp(
    Y = Y, X = cbind(1,X.mat), Z = Z.mat, mixing.distribution = "lnorm"
  )
  ### Fit using formula interface
  dat <- data.frame(y = Y, x1 = X.mat[,1], x2 = X.mat[,2], z1 = Z.mat[,1], z2 = Z.mat[,2])
  fit.formula <- fit_zimp(
    y ~ 1 + x1 + x2|z1 + z2,
    data = dat,
    mixing.distribution = "lnorm"
  )
  
  ### Compare results
  expect_equal(unname(fit.formula$be), unname(fit.matrix$be),tolerance = 1e-8)
  expect_equal(fit.formula$phi,fit.matrix$phi,tolerance = 1e-8)
  expect_equal(unname(fit.formula$gam), unname(fit.matrix$gam),tolerance = 1e-8)
  expect_equal(fit.formula$log_lik_value,fit.matrix$log_lik_value,tolerance = 1e-8)
  expect_equal(fit.formula$AIC,fit.matrix$AIC,tolerance = 1e-8)
  expect_equal(fit.formula$BIC,fit.matrix$BIC,tolerance = 1e-8)
  
  #-------------------------------------------------------------------------------#
  
  
  ## Test 2.c. WITHOUT Intercept + WITH Intercept ---------------------------#
  
  ### Simulation ZIMP by lnorm
  be <- c(0.3, 0.7)
  phi <- 0.2
  mu.vec <- exp(X.mat %*% be)
  gam <- c(0.5, 0.4, 0.6)
  pr.vec <- exp(cbind(1,Z.mat) %*% gam)/(1+exp(cbind(1,Z.mat) %*% gam))
  
  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = dist.simu,
                 zero_inflation = TRUE, prob = pr.vec)$Y
  
  ### Fit using matrix interface
  fit.matrix <- fit_zimp(
    Y = Y, X = X.mat, Z = cbind(1,Z.mat), mixing.distribution = "lnorm"
  )
  ### Fit using formula interface
  dat <- data.frame(y = Y, x1 = X.mat[,1], x2 = X.mat[,2], z1 = Z.mat[,1], z2 = Z.mat[,2])
  fit.formula <- fit_zimp(
    y ~ x1 + x2|1 + z1 + z2,
    data = dat,
    mixing.distribution = "lnorm"
  )
  
  ### Compare results
  expect_equal(unname(fit.formula$be), unname(fit.matrix$be),tolerance = 1e-8)
  expect_equal(fit.formula$phi,fit.matrix$phi,tolerance = 1e-8)
  expect_equal(unname(fit.formula$gam), unname(fit.matrix$gam),tolerance = 1e-8)
  expect_equal(fit.formula$log_lik_value,fit.matrix$log_lik_value,tolerance = 1e-8)
  expect_equal(fit.formula$AIC,fit.matrix$AIC,tolerance = 1e-8)
  expect_equal(fit.formula$BIC,fit.matrix$BIC,tolerance = 1e-8)
  
  #-------------------------------------------------------------------------------#
  
  ## Test 2.d. WITH Intercept + WITH Intercept ---------------------------#
  
  ### Simulation ZIMP by lnorm
  be <- c(0.5, 0.3, 0.7)
  phi <- 0.2
  mu.vec <- exp(cbind(1,X.mat) %*% be)
  gam <- c(0.5, 0.4, 0.6)
  pr.vec <- exp(cbind(1,Z.mat) %*% gam)/(1+exp(cbind(1,Z.mat) %*% gam))
  
  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = dist.simu,
                 zero_inflation = TRUE, prob = pr.vec)$Y
  
  ### Fit using matrix interface
  fit.matrix <- fit_zimp(
    Y = Y, X = cbind(1,X.mat), Z = cbind(1,Z.mat), mixing.distribution = "lnorm"
  )
  ### Fit using formula interface
  dat <- data.frame(y = Y, x1 = X.mat[,1], x2 = X.mat[,2], z1 = Z.mat[,1], z2 = Z.mat[,2])
  fit.formula <- fit_zimp(
    y ~ 1 + x1 + x2|1 + z1 + z2,
    data = dat,
    mixing.distribution = "lnorm"
  )
  
  ### Compare results
  expect_equal(unname(fit.formula$be), unname(fit.matrix$be),tolerance = 1e-8)
  expect_equal(fit.formula$phi,fit.matrix$phi,tolerance = 1e-8)
  expect_equal(unname(fit.formula$gam), unname(fit.matrix$gam),tolerance = 1e-8)
  expect_equal(fit.formula$log_lik_value,fit.matrix$log_lik_value,tolerance = 1e-8)
  expect_equal(fit.formula$AIC,fit.matrix$AIC,tolerance = 1e-8)
  expect_equal(fit.formula$BIC,fit.matrix$BIC,tolerance = 1e-8)
  
  #-------------------------------------------------------------------------------#

})

test_that("fit_zimp handles count offset vectors correctly", {
  
  set.seed(12345)
  
  # Input data
  n <- 100
  offset <- runif(n, -0.5, 0.5)
  X.mat <- matrix(rnorm(2*n),n,2)
  X.mat <- cbind(1, X.mat)
  colnames(X.mat) <- c("x0", "x1", "x2")
  
  # Input coefficients
  be <- c(0.5, 0.3, 0.7)
  phi <- 0.2
  mu.vec <- exp(X.mat %*% be + offset)
  
  # Simulation
  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = "lnorm")$Y
  
  # Fit with offset vector
  fit.offset <- fit_zimp(Y = Y, X = X.mat, offset.count = offset,
                         mixing.distribution = "lnorm")
  
  # Basic checks
  expect_true(is.finite(fit.offset$log_lik_value))
  expect_length(fit.offset$be, ncol(X.mat))
  expect_true(all(is.finite(fit.offset$be)))
  expect_true(is.finite(fit.offset$phi))
  expect_true(fit.offset$phi > 0)
  expect_true(all(is.finite(fit.offset$hessian)))
  expect_equal(fit.offset$offset.count, offset, tolerance = 1e-12)
})

test_that("fit_zimp rejects invalid inputs", {
  
  set.seed(12345)
  Y <- rpois(n, 2)
  X <- cbind(1, matrix(rnorm(2*n), n, 2))
  Z <- cbind(1, matrix(rnorm(2*n), n, 2))
  
  # Invalid response
  expect_error(fit_zimp(Y = c(-1, Y[-1]), X = X))
  expect_error(fit_zimp(Y = c(1.5, Y[-1]), X = X))
  expect_error(fit_zimp(Y = c(NA, Y[-1]), X = X))
  expect_error(fit_zimp(Y = c(Inf, Y[-1]), X = X))
  
  # Invalid design matrices
  expect_error(fit_zimp(Y = Y, X = X[-1, ]))
  expect_error(fit_zimp(Y = Y, X = X, Z = Z[-1, ]))
  expect_error(fit_zimp(Y = Y, X = X, Z = Z, offset.count = rep(0, n - 1)))
  
  # Invalid model specification
  expect_error(fit_zimp(Y = Y, X = X, mixing.distribution = "invalid"))
  expect_error(fit_zimp(Y = Y, X = X, mixing.distribution = "lnorm", var.power = 4))
  
  # Invalid initial parameters
  expect_error(fit_zimp(Y = Y, X = X, ini.par = c(0, 0)))
  
})