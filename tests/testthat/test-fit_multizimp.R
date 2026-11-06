# ================================================================
# TEST fit_multizimp
# ================================================================

test_that("fit_multizimp supports direct MP and ZIMP inputs with matrix and list designs", {
  set.seed(12345)
  n <- 15; k <- 3; J <- 1
  psi.N <- simu_psi_N(k, J)
  
  ## 1. MP: common X matrix
  X <- cbind(1, matrix(rnorm(2*n), n, 2)); beta <- c(0.5, 0.3, 0.7)
  sim <- simu_multizimp(X=X, beta=beta, Z=NULL, mixing.distribution=rep("gamma", k), psi_N=psi.N)
  fit <- fit_multizimp(Y=sim$Y, X=X, Z=NULL, mixing.distribution="gamma", var.power=2, J=J)
  
  expect_equal(fit$n, n); expect_equal(fit$k, k); expect_equal(length(fit$marginal), k)
  for (j in seq_len(k)) {
    expect_false(fit$marginal[[j]]$zero_inflation); expect_equal(length(fit$marginal[[j]]$be), ncol(X)); expect_null(fit$marginal[[j]]$gam)
    expect_true(all(is.finite(fit$marginal[[j]]$be))); expect_true(is.finite(fit$marginal[[j]]$phi))
  }
  
  ## 2. ZIMP: common X and Z matrices
  Z <- cbind(1, rnorm(n)); gamma <- c(0.5, 0.3); psi.B <- simu_psi_B(k)
  sim <- simu_multizimp(X=X, beta=beta, Z=Z, gamma=gamma, zero_inflation=rep(TRUE, k), mixing.distribution=rep("gamma", k), psi_N=psi.N, psi_B=psi.B)
  fit <- fit_multizimp(Y=sim$Y, X=X, Z=Z, mixing.distribution="gamma", var.power=2, J=J)
  
  expect_equal(fit$n, n); expect_equal(fit$k, k); expect_equal(length(fit$marginal), k)
  for (j in seq_len(k)) {
    expect_true(fit$marginal[[j]]$zero_inflation); expect_equal(length(fit$marginal[[j]]$be), ncol(X)); expect_equal(length(fit$marginal[[j]]$gam), ncol(Z))
    expect_true(all(is.finite(fit$marginal[[j]]$be))); expect_true(all(is.finite(fit$marginal[[j]]$gam))); expect_true(is.finite(fit$marginal[[j]]$phi))
  }
  
  ## 3. MP: margin-specific X lists
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  X.list <- list(cbind(1, x1), cbind(1, x1, x2), cbind(1, x1, x2, x3))
  beta.list <- list(c(0.5, 0.3), c(0.4, 0.2, -0.1), c(0.3, 0.2, -0.1, 0.15))
  Y <- matrix(0L, n, k)
  for (j in seq_len(k)) { mu <- exp(as.vector(X.list[[j]] %*% beta.list[[j]])); Y[,j] <- rpois(n, mu) }
  
  fit <- fit_multizimp(Y=Y, X=X.list, Z=NULL, mixing.distribution="gamma", var.power=2, J=J)
  
  expect_equal(fit$n, n); expect_equal(fit$k, k); expect_equal(length(fit$marginal), k)
  for (j in seq_len(k)) {
    expect_false(fit$marginal[[j]]$zero_inflation); expect_equal(length(fit$marginal[[j]]$be), ncol(X.list[[j]])); expect_null(fit$marginal[[j]]$gam)
    expect_true(all(is.finite(fit$marginal[[j]]$be))); expect_true(is.finite(fit$marginal[[j]]$phi))
  }
  
  ## 4. ZIMP: margin-specific X and Z lists
  z1 <- rnorm(n); z2 <- rnorm(n)
  Z.list <- list(matrix(1, n, 1), cbind(1, z1), cbind(1, z1, z2))
  gamma.list <- list(0.5, c(0.4, 0.3), c(0.3, 0.2, -0.2))
  Y <- matrix(0L, n, k)
  for (j in seq_len(k)) { mu <- exp(as.vector(X.list[[j]] %*% beta.list[[j]])); pb <- plogis(as.vector(Z.list[[j]] %*% gamma.list[[j]])); Y[,j] <- rbinom(n, 1, pb) * rpois(n, mu) }
  
  fit <- fit_multizimp(Y=Y, X=X.list, Z=Z.list, mixing.distribution="gamma", var.power=2, J=J)
  
  expect_equal(fit$n, n); expect_equal(fit$k, k); expect_equal(length(fit$marginal), k)
  for (j in seq_len(k)) {
    expect_true(fit$marginal[[j]]$zero_inflation); expect_equal(length(fit$marginal[[j]]$be), ncol(X.list[[j]])); expect_equal(length(fit$marginal[[j]]$gam), ncol(Z.list[[j]]))
    expect_true(all(is.finite(fit$marginal[[j]]$be))); expect_true(all(is.finite(fit$marginal[[j]]$gam))); expect_true(is.finite(fit$marginal[[j]]$phi))
  }
})



test_that("fit_multizimp accepts marginal models from select_multi_margin", {
  set.seed(12345)
  n <- 15; k <- 3; J <- 1
  X <- cbind(1, rnorm(n)); beta <- c(0.5, 0.3); psi.N <- simu_psi_N(k, J)
  
  ## 1. MP
  sim <- simu_multizimp(X=X, beta=beta, Z=NULL, mixing.distribution=rep("gamma", k), psi_N=psi.N)
  sel <- select_multi_margin(Y=sim$Y, X=X, Z=NULL, mixing.distribution=c("gamma","lnorm"), var.power=2, compute.I=TRUE, I.d0=3, parallel=1)
  fit <- fit_multizimp(margin=sel, J=J, d0=3)
  
  expect_s3_class(sel, "multi_margin_fit"); expect_s3_class(fit, "multi_zipm_fit")
  expect_equal(fit$n, n); expect_equal(fit$k, k)
  expect_equal(length(fit$marginal), k)
  for (j in seq_len(k)) { expect_false(fit$marginal[[j]]$zero_inflation); expect_equal(fit$marginal[[j]]$distribution, sel$selected_margins[[j]]$distribution); expect_equal(fit$marginal[[j]]$var.power, sel$selected_margins[[j]]$var.power); expect_equal(fit$marginal[[j]]$be, sel$selected_margins[[j]]$be); expect_equal(fit$marginal[[j]]$phi, sel$selected_margins[[j]]$phi); expect_null(fit$marginal[[j]]$gam) }
  
  ## 2. ZIMP
  Z <- cbind(1, rnorm(n)); gamma <- c(0.5, 0.4); psi.B <- simu_psi_B(k)
  sim <- simu_multizimp(X=X, beta=beta, Z=Z, gamma=gamma, zero_inflation=rep(TRUE, k), mixing.distribution=rep("gamma", k), psi_N=psi.N, psi_B=psi.B)
  sel <- select_multi_margin(Y=sim$Y, X=X, Z=Z, mixing.distribution=c("gamma","lnorm"), var.power=2, compute.I=TRUE, I.d0=3, parallel=1)
  fit <- fit_multizimp(margin=sel, J=J, d0=3)
  
  expect_s3_class(sel, "multi_margin_fit"); expect_s3_class(fit, "multi_zipm_fit")
  expect_equal(fit$n, n); expect_equal(fit$k, k)
  expect_equal(length(fit$marginal), k)
  for (j in seq_len(k)) { expect_true(fit$marginal[[j]]$zero_inflation); expect_equal(fit$marginal[[j]]$distribution, sel$selected_margins[[j]]$distribution); expect_equal(fit$marginal[[j]]$var.power, sel$selected_margins[[j]]$var.power); expect_equal(fit$marginal[[j]]$be, sel$selected_margins[[j]]$be); expect_equal(fit$marginal[[j]]$phi, sel$selected_margins[[j]]$phi); expect_equal(fit$marginal[[j]]$gam, sel$selected_margins[[j]]$gam) }
})

test_that("fit_multizimp estimates count dependence after marginal selection", {
  set.seed(12345)
  n <- 15; k <- 3; J <- 1
  X <- cbind(1, rnorm(n)); beta <- c(0.5, 0.3); psi.N <- simu_psi_N(k, J)
  
  sim <- simu_multizimp(X=X, beta=beta, Z=NULL, mixing.distribution=rep("gamma", k), psi_N=psi.N)
  sel <- select_multi_margin(Y=sim$Y, X=X, Z=NULL, mixing.distribution="gamma", var.power=2, compute.I=TRUE, I.d0=3, parallel=1)
  fit <- fit_multizimp(margin=sel, J=J, d0=3)
  
  expect_true(!is.null(fit$count_layer))
  expect_true(!is.null(fit$count_layer$psi))
  expect_equal(dim(fit$count_layer$psi), c(k, J))
  expect_true(all(is.finite(fit$count_layer$psi)))
  expect_true(all(rowSums(fit$count_layer$psi^2) < 1))
  expect_true(is.finite(fit$summary$fit$logLik))
})


test_that("fit_multizimp estimates count dependence with different factor dimensions", {
  set.seed(12345)
  n <- 15; k <- 4
  X <- cbind(1, rnorm(n)); 
  beta <- c(0.5, 0.3); psi.N <- simu_psi_N(k, 1)
  
  sim <- simu_multizimp(X=X, beta=beta, Z=NULL, mixing.distribution=rep("gamma", k), psi_N=psi.N)
  sel <- select_multi_margin(Y=sim$Y, X=X, Z=NULL, mixing.distribution="gamma", var.power=2, compute.I=TRUE, I.d0=3, parallel=1)
  
  fit1 <- fit_multizimp(margin=sel, J=1, d0=3)
  expect_equal(dim(fit1$count_layer$psi), c(k, 1)); expect_true(all(is.finite(fit1$count_layer$psi))); 
  expect_true(all(rowSums(fit1$count_layer$psi^2) < 1))
  
  # fit2 <- fit_multizimp(margin=sel, J=2, d0=3)
  # expect_equal(dim(fit2$count_layer$psi), c(k, 2)); expect_true(all(is.finite(fit2$count_layer$psi))); 
  # expect_true(all(rowSums(fit2$count_layer$psi^2) < 1))
  
  fit3 <- fit_multizimp(margin=sel, J=3, d0=3)
  expect_equal(dim(fit3$count_layer$psi), c(k, 3)); expect_true(all(is.finite(fit3$count_layer$psi))); 
  expect_true(all(rowSums(fit3$count_layer$psi^2) < 1))
  
})


test_that("fit_multizimp estimates count dependence for ZIMP margins with L = 0", {
  set.seed(12345)
  n <- 15; k <- 4
  X <- cbind(1, rnorm(n)); Z <- cbind(1, rnorm(n))
  beta <- c(0.5, 0.3); gamma <- c(0.5, 0.3); psi.N <- simu_psi_N(k, 1); psi.B <- simu_psi_B(k)
  
  sim <- simu_multizimp(X=X, beta=beta, Z=Z, gamma=gamma, zero_inflation=rep(TRUE, k), mixing.distribution=rep("gamma", k), 
                        psi_N=psi.N, psi_B = psi.B)
  
  sel <- select_multi_margin(Y=sim$Y, X=X, Z=Z, mixing.distribution="gamma", var.power=2, compute.I=TRUE, I.d0=3, parallel=1)
  
  fit1 <- fit_multizimp(margin=sel, J=1, L=0, d0=3)
  expect_equal(dim(fit1$count_layer$psi), c(k, 1)); 
  expect_true(all(is.finite(fit1$count_layer$psi))); expect_true(all(rowSums(fit1$count_layer$psi^2) < 1))
  
  # fit2 <- fit_multizimp(margin=sel, J=2, L=0, d0=3)
  # expect_equal(dim(fit2$count_layer$psi), c(k, 2)); 
  # expect_true(all(is.finite(fit2$count_layer$psi))); expect_true(all(rowSums(fit2$count_layer$psi^2) < 1))
  
  fit3 <- fit_multizimp(margin=sel, J=3, L=0, d0=3)
  expect_equal(dim(fit3$count_layer$psi), c(k, 3)); 
  expect_true(all(is.finite(fit3$count_layer$psi))); expect_true(all(rowSums(fit3$count_layer$psi^2) < 1))
})


test_that("fit_multizimp estimates count dependence for ZIMP margins with L = 1", {
  set.seed(12345)
  n <- 15; k <- 4
  X <- cbind(1, rnorm(n)); Z <- cbind(1, rnorm(n))
  beta <- c(0.5, 0.3); gamma <- c(0.5, 0.3); psi.N <- simu_psi_N(k, 1); psi.B <- simu_psi_B(k)
  
  sim <- simu_multizimp(X=X, beta=beta, Z=Z, gamma=gamma, zero_inflation=rep(TRUE, k), mixing.distribution=rep("gamma", k), 
                        psi_N=psi.N, psi_B = psi.B)
  
  sel <- select_multi_margin(Y=sim$Y, X=X, Z=Z, mixing.distribution="gamma", var.power=2, compute.I=TRUE, I.d0=3, parallel=1)
  
  fit1 <- fit_multizimp(margin=sel, J=1, L=1, d0=3)
  expect_equal(dim(fit1$count_layer$psi), c(k, 1)); 
  expect_true(all(is.finite(fit1$count_layer$psi))); expect_true(all(rowSums(fit1$count_layer$psi^2) < 1))
  
  # fit2 <- fit_multizimp(margin=sel, J=2, L=1, d0=3)
  # expect_equal(dim(fit2$count_layer$psi), c(k, 2)); 
  # expect_true(all(is.finite(fit2$count_layer$psi))); expect_true(all(rowSums(fit2$count_layer$psi^2) < 1))
  
  fit3 <- fit_multizimp(margin=sel, J=3, L=1, d0=3)
  expect_equal(dim(fit3$count_layer$psi), c(k, 3)); 
  expect_true(all(is.finite(fit3$count_layer$psi))); expect_true(all(rowSums(fit3$count_layer$psi^2) < 1))
})

test_that("fit_multizimp handles two ZIMP margins and one MP margin with L = 0 and L = 1", {
  set.seed(12345)
  n <- 20; k <- 3; J <- 1
  X <- cbind(1, rnorm(n)); Z1 <- cbind(1, rnorm(n)); Z2 <- cbind(1, rnorm(n))
  beta <- c(0.5, 0.3); gamma1 <- c(0.5, 0.3); gamma2 <- c(0.4, -0.2)
  psi.N <- simu_psi_N(k, J); psi.B <- simu_psi_B(2)
  
  mu <- exp(as.vector(X %*% beta)); pb1 <- plogis(as.vector(Z1 %*% gamma1)); pb2 <- plogis(as.vector(Z2 %*% gamma2))
  Y <- matrix(0L, n, k); Y[,1] <- rbinom(n, 1, pb1) * rpois(n, mu); Y[,2] <- rbinom(n, 1, pb2) * rpois(n, mu); Y[,3] <- rpois(n, mu)
  
  sel <- select_multi_margin(Y=Y, X=list(X,X,X), Z=list(Z1,Z2,NULL),  mixing.distribution="gamma", var.power=2, compute.I=TRUE, I.d0=3, parallel=1)
  
  expect_true(sel$selected_margins[[1]]$zero_inflation); expect_true(sel$selected_margins[[2]]$zero_inflation); expect_false(sel$selected_margins[[3]]$zero_inflation)
  
  ## L = 0: no dependence in zero layer
  fit0 <- fit_multizimp(margin=sel, J=J, L=0, d0=3)
  expect_equal(dim(fit0$count_layer$psi), c(k, J)); expect_true(all(is.finite(fit0$count_layer$psi))); expect_true(all(rowSums(fit0$count_layer$psi^2) < 1))
  expect_true(fit0$marginal[[1]]$zero_inflation); expect_true(fit0$marginal[[2]]$zero_inflation); expect_false(fit0$marginal[[3]]$zero_inflation)
  
  ## L = 1: dependence between the two ZIMP zero indicators
  fit1 <- fit_multizimp(margin=sel, J=J, L=1, d0=3)
  expect_equal(dim(fit1$count_layer$psi), c(k, J)); expect_true(all(is.finite(fit1$count_layer$psi))); expect_true(all(rowSums(fit1$count_layer$psi^2) < 1))
  expect_true(fit1$marginal[[1]]$zero_inflation); expect_true(fit1$marginal[[2]]$zero_inflation); expect_false(fit1$marginal[[3]]$zero_inflation)
  expect_true(!is.null(fit1$zero_layer))
})

