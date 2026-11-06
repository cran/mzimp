# ============================================================
# Tests for select_zimp()
# ============================================================

start <- proc.time()
n <- 100
test_that("select_zimp and fit_zimp provide the same result when there is only one option", {
  
  set.seed(12345)
  
  # Input data
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  colnames(X.mat) <- c("x0", "x1", "x2")
  
  # Input coefficients
  be <- c(0.5, 0.3, 0.7)
  phi <- 0.2
  mu.vec <- exp(X.mat %*% be)
  
  # Mixing distributions and variance powers
  dist.check <- c("lnorm")
  
  # Simulation by lognormal
  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = "lnorm")$Y
  
  # Test all mixing distribution × variance power combinations
  fit <- fit_zimp(Y = Y, X = X.mat, mixing.distribution = dist.check)
  select <- select_zimp(Y = Y, X = X.mat, mixing.distribution = dist.check)    
  
  best <- select$best_model
  
  expect_equal(best$be, fit$be, tolerance = 1e-10)
  expect_equal(best$phi, fit$phi, tolerance = 1e-10)
  
  expect_equal(best$par, fit$par, tolerance = 1e-10)
  
  expect_equal(
    best$log_lik_value,
    fit$log_lik_value,
    tolerance = 1e-10
  )
})
    
test_that("select_zimp supports all mixing distributions and variance powers and selects the best model", {
  
  set.seed(12345)
  
  # Input data
  X.mat <- cbind(1, matrix(rnorm(n * 2), n, 2))
  colnames(X.mat) <- c("x0", "x1", "x2")
  
  # Input coefficients
  be <- c(0.5, 0.3, 0.7)
  phi <- 0.2
  mu.vec <- exp(X.mat %*% be)
  
  # Mixing distributions and variance powers
  dist.vec <- c("lnorm", "gamma", "ig","weibull", "loglogistic", "betaprime")
  var.power.vec <- c(1, 2, 3)
  
  # Simulation by lognormal
  Y <- simu_zipm(n, mu = mu.vec, phi = phi, latent.distribution = "lnorm")$Y
  
  # Fit all candidate models
  select <- select_zimp(Y = Y, X = X.mat, mixing.distribution = dist.vec, var.power = var.power.vec)
  
  tab <- select$summary_table
  best <- select$best_model
  
  # 6 distributions × 3 variance powers = 18 candidate models
  expect_equal(nrow(tab), length(dist.vec) * length(var.power.vec))
  
  # All requested distributions were fitted
  expect_setequal(
    unique(tab$distribution),
    dist.vec
  )
  
  # All requested variance powers were fitted
  expect_setequal(
    unique(tab$var.power),
    var.power.vec
  )
  
  # All candidate models returned finite log-likelihoods
  expect_true(all(is.finite(tab$log_lik)))
  expect_true(all(is.finite(tab$neg_log_lik)))
  expect_true(all(is.finite(tab$phi)))
  expect_true(all(tab$phi > 0))
  
  # Best model has the largest log-likelihood
  expect_equal(
    best$log_lik_value,
    max(tab$log_lik),
    tolerance = 1e-10
  )
  
  # Best model agrees with the first row of the ranked summary table
  expect_equal(
    best$mixing.distribution,
    tab$distribution[1]
  )
  
  expect_equal(
    best$var.power,
    tab$var.power[1]
  )
  
  expect_equal(
    best$log_lik_value,
    tab$log_lik[1],
    tolerance = 1e-10
  )
  
  expect_equal(
    best$phi,
    tab$phi[1],
    tolerance = 1e-10
  )
  
  # Best model contains valid estimates and inference
  expect_length(best$be, ncol(X.mat))
  expect_length(best$phi, 1)
  expect_true(all(is.finite(best$be)))
  expect_true(is.finite(best$phi))
  expect_true(best$phi > 0)
  expect_true(all(is.finite(best$hessian)))
})

test_that("select_zimp supports all mixing distributions and variance powers and selects the best ZIMP model", {
  
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
  dist.vec <- c("lnorm", "gamma", "ig",
                "weibull", "loglogistic", "betaprime")
  var.power.vec <- c(1, 2, 3)
  
  # Simulation by lognormal
  Y <- simu_zipm(
    n,
    mu = mu.vec,
    phi = phi,
    latent.distribution = "lnorm",
    zero_inflation = TRUE,
    prob = pr.vec
  )$Y
  
  # Fit all candidate ZIMP models
  select <- select_zimp(
    Y = Y,
    X = X.mat,
    Z = Z.mat,
    mixing.distribution = dist.vec,
    var.power = var.power.vec
  )
  
  tab <- select$summary_table
  best <- select$best_model
  
  # 6 distributions × 3 variance powers = 18 candidate models
  expect_equal(
    nrow(tab),
    length(dist.vec) * length(var.power.vec)
  )
  
  # All requested distributions were fitted
  expect_setequal(
    unique(tab$distribution),
    dist.vec
  )
  
  # All requested variance powers were fitted
  expect_setequal(
    unique(tab$var.power),
    var.power.vec
  )
  
  # All candidate models returned valid results
  expect_true(all(is.finite(tab$log_lik)))
  expect_true(all(is.finite(tab$neg_log_lik)))
  expect_true(all(is.finite(tab$phi)))
  expect_true(all(tab$phi > 0))
  
  # Best model has the largest log-likelihood
  expect_equal(
    best$log_lik_value,
    max(tab$log_lik),
    tolerance = 1e-10
  )
  
  # Best model agrees with the first row of the ranked summary table
  expect_equal(
    best$mixing.distribution,
    tab$distribution[1]
  )
  
  expect_equal(
    best$var.power,
    tab$var.power[1]
  )
  
  expect_equal(
    best$log_lik_value,
    tab$log_lik[1],
    tolerance = 1e-10
  )
  
  expect_equal(
    best$phi,
    tab$phi[1],
    tolerance = 1e-10
  )
  
  # Best model has the correct parameter dimensions
  expect_length(best$be, ncol(X.mat))
  expect_length(best$gam, ncol(Z.mat))
  expect_length(best$phi, 1)
  
  # Best model contains finite parameter estimates
  expect_true(all(is.finite(best$be)))
  expect_true(all(is.finite(best$gam)))
  expect_true(is.finite(best$phi))
  expect_true(best$phi > 0)
  
  # Inference was computed for the selected model
  expect_true(all(is.finite(best$hessian)))
  
})

test_that("select_zimp supports warm-start distributions without gamma as candidate", {
  
  
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
  dist.vec <- c("lnorm", "gamma", "ig",
                "weibull", "loglogistic", "betaprime")
  var.power.vec <- c(1, 2, 3)
  
  # Simulation by lognormal
  Y <- simu_zipm(
    n,
    mu = mu.vec,
    phi = phi,
    latent.distribution = "lnorm",
    zero_inflation = TRUE,
    prob = pr.vec
  )$Y
  
  for (dist.check in c("weibull", "loglogistic", "betaprime")) {
    
    select <- select_zimp(
      Y = Y,
      X = X.mat,
      mixing.distribution = dist.check,
      var.power = c(1, 2, 3)
    )
    
    expect_equal(nrow(select$summary_table), 3)
    
    expect_true(
      all(select$summary_table$distribution == dist.check)
    )
    
    expect_false(
      "gamma" %in% select$summary_table$distribution
    )
    
    expect_true(
      all(is.finite(select$summary_table$log_lik))
    )
    
    expect_equal(
      select$best_model$log_lik_value,
      max(select$summary_table$log_lik),
      tolerance = 1e-10
    )
  }
})

test_that("select_zimp returns the expected compact output", {
  
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
  dist.vec <- c("lnorm", "gamma", "ig",
                "weibull", "loglogistic", "betaprime")
  var.power.vec <- c(1, 2, 3)
  
  # Simulation by lognormal
  Y <- simu_zipm(
    n,
    mu = mu.vec,
    phi = phi,
    latent.distribution = "lnorm",
    zero_inflation = TRUE,
    prob = pr.vec
  )$Y
  
  select <- select_zimp(
    Y = Y,
    X = X.mat,
    mixing.distribution = c("gamma", "lnorm"),
    var.power = 2
  )
  
  expect_named(
    select,
    c("summary_table", "best_model", "summary", "call"),
    ignore.order = FALSE
  )
  
  expect_s3_class(select, "zimp_selection")
  expect_s3_class(select$best_model, "zimp_fit")
  expect_s3_class(select$summary, "zimp_selection_summary")
})

test_that("select_zimp formula and matrix interfaces give the same selection", {
  
  
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
  dist.vec <- c("lnorm", "gamma", "ig",
                "weibull", "loglogistic", "betaprime")
  var.power.vec <- c(1, 2, 3)
  
  # Simulation by lognormal
  Y <- simu_zipm(
    n,
    mu = mu.vec,
    phi = phi,
    latent.distribution = "lnorm",
    zero_inflation = TRUE,
    prob = pr.vec
  )$Y
  
  dat <- data.frame(
    y = Y,
    x1 = X.mat[, 2],
    x2 = X.mat[, 3]
  )
  
  select.mat <- select_zimp(
    Y = Y,
    X = X.mat,
    mixing.distribution = c("gamma", "lnorm", "ig"),
    var.power = c(1, 2, 3)
  )
  
  select.form <- select_zimp(
    y ~ 1 + x1 + x2,
    data = dat,
    mixing.distribution = c("gamma", "lnorm", "ig"),
    var.power = c(1, 2, 3)
  )
  
  expect_equal(
    select.form$summary_table$log_lik,
    select.mat$summary_table$log_lik,
    tolerance = 1e-10
  )
  
  expect_equal(
    unname(select.form$best_model$par),
    unname(select.mat$best_model$par),
    tolerance = 1e-10
  )
})

proc.time() - start