#' mzimp: Zero-Inflated Mixed Poisson Models
#'
#' @useDynLib mzimp, .registration = TRUE
#' @importFrom Rcpp evalCpp
#'
#' @importFrom stats binomial coef dnorm glm model.matrix nlminb
#' @importFrom stats pgamma plogis pnorm poisson qbeta qgamma
#' @importFrom stats qlnorm qlogis qnorm qt qweibull
#' @importFrom stats rbinom rexp rgamma rnorm rpois runif sd uniroot var
#'
#' @importFrom utils flush.console getFromNamespace head tail
#'
#' @keywords internal
"_PACKAGE"