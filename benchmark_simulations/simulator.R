#' Simulation script for Leaky IV project
#' 
#' This function simulates data from a linear structural equation model with 
#' some leaky instrument(s).
#' 
#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param z_cnt Logical. If \code{TRUE}, Z is multivariate normal; else, 
#'   multivariate binomial.
#' @param z_rho Autocorrelation of the Toeplitz matrix for Z.
#' @param rho Correlation between residuals for X and Y.
#' @param theta Average treatment effect of X on Y.
#' @param r2_x Proportion of variance explained for X.
#' @param r2_y Proportion of variance explained for Y.
#' @param pr_valid Proportion of candidate instruments that are valid, 
#'   i.e. have no direct effect on Y.

# Simulate data (inspired by Hartwig et al., 2017)
sim_dat <- function(n, d_z, z_cnt, z_rho, rho, theta, r2_x, r2_y, pr_valid) {
  # What proportion are valid?
  valid_cnt <- round(pr_valid * d_z)
  pr_valid <- valid_cnt / d_z
  # Draw Z's
  if (z_rho != 0) {
    Sigma_z <- toeplitz(z_rho^(0:(d_z - 1)))
  }
  if (isTRUE(z_cnt)) {
    z <- matrix(rnorm(n * d_z), ncol = d_z)
    if (z_rho == 0) {
      Sigma_z <- diag(rep(1, d_z))
    } else {
      z <- z %*% chol(Sigma_z)
    }
    z <- z / d_z
  } else {
    pr_z <- runif(d_z, min = 0.1, max = 0.9)
    var_z <- pr_z * (1 - pr_z)
    if (z_rho == 0) {
      z <- sapply(pr_z, function(p) rbinom(n, 1, p))
      Sigma_z <- diag(var_z)
    } else {
      z <- draw.correlated.binary(no.row = n, d = d_z, prop.vec = pr_z, 
                                  corr.mat = Sigma_z)
      Sigma_z <- diag(sqrt(var_z)) %*% Sigma_z %*% diag(sqrt(var_z))
    }
  }
  colnames(z) <- paste0('z', seq_len(d_z))
  # Simulate standardized residual vectors
  Sigma_eps <- matrix(c(1, rho, rho, 1), ncol = 2)
  eps <- matrix(rnorm(n * 2), ncol = 2)
  eps <- eps %*% chol(Sigma_eps)
  # Draw random weights for Z (adapted from Hartford et al., 2020)
  # Then calculate eta_x from r2_x
  nu1 <- runif(d_z, min = 0.01, max = 0.2) 
  sigma_zx <- sd(sqrt(0.1) * as.numeric(z %*% nu1))
  beta <- (sqrt(0.1) / sigma_zx) * nu1
  var_x <- as.numeric(t(beta) %*% Sigma_z %*% beta) / r2_x 
  eta_x <- sqrt(var_x * (1 - r2_x))
  x <- as.numeric(z %*% beta) + eps[, 1] * eta_x
  # And again for Y, although we need to solve a quadratic equation for eta_y
  nu2 <- runif(d_z, min = 0.01, max = 0.2)
  sigma_zy <- sd(sqrt(0.1) * as.numeric(z %*% nu2))
  gamma <- (sqrt(0.1) / sigma_zy) * nu2 * (1 - pr_valid)
  if (pr_valid > 0) {
    gamma[sample(d_z, valid_cnt)] <- 0
  }
  var_mu <- as.numeric(t(gamma) %*% Sigma_z %*% gamma) + theta^2 * var_x + 
    2 * theta * as.numeric(beta %*% Sigma_z %*% gamma)
  var_y <- var_mu / r2_y
  b <- 2 * theta * rho * eta_x
  eta_y <- (-b + sqrt(b^2 - 4 * (var_mu - var_y))) / 2
  y <- as.numeric(z %*% gamma) + x * theta + eps[, 2] * eta_y
  # Export results
  dat <- data.table(z, 'x' = x, 'y' = y)
  #fwrite(dat, paste0('./simulations/sim', s_idx, '.csv'))
  params <- list(
    'theta' = theta, 'beta' = beta, 'gamma' = gamma, 'eta_x' = eta_x, 
    'eta_y' = eta_y, 'rho' = rho, 'r2_x' = r2_x, 'r2_y' = r2_y
  )
  #saveRDS(params, paste0('./simulations/params', s_idx, '.rds'))
  out <- list('dat' = dat, 'params' = params)
  return(out)
}




