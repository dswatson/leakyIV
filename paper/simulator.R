#' Simulation scripts for Leaky IV project
#' 
#' This function simulates data from a linear SEM with some leaky instrument(s).
#' 
#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param z_rho Autocorrelation of the Toeplitz matrix for Z.
#' @param rho Correlation between residuals for X and Y.
#' @param theta Average treatment effect of X on Y.
#' @param snr_x Signal to noise ratio for X.
#' @param snr_y Signal to noise ratio for Y.
#' @param pr_valid Proportion of candidate instruments that are valid, 
#'   i.e. have no direct effect on Y.
#'   

sim_dat <- function(n, d_z, z_rho, rho, theta, snr_x, snr_y, pr_valid) {
  
  # What proportion are valid?
  valid_cnt <- round(pr_valid * d_z)
  pr_valid <- valid_cnt / d_z
  
  # Draw Z's
  var_z <- 1 / d_z
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  Sigma_z <- toeplitz(z_rho^(0:(d_z - 1))) * var_z
  z <- z %*% chol(Sigma_z)
  colnames(z) <- paste0('z', seq_len(d_z))
  
  # Simulate standardized residual vectors
  Sigma_eps <- matrix(c(1, rho, rho, 1), ncol = 2)
  eps <- matrix(rnorm(n * 2), ncol = 2)
  eps <- eps %*% chol(Sigma_eps)
  
  # Draw random weights for beta, calculate var_x and eta_x
  pve_x <- snr_x / (1 + snr_x)
  beta <- rnorm(d_z)
  var_x <- as.numeric(t(beta) %*% Sigma_z %*% beta) / pve_x 
  eta_x <- sqrt(var_x * (1 - pve_x))
  x <- as.numeric(z %*% beta) + eps[, 1] * eta_x
  
  # Draw random weights for gamma, tune magnitude to ensure target SNR
  var_y <- 10
  pve_y <- snr_y / (1 + snr_y)
  gamma <- rnorm(d_z)
  if (pr_valid > 0) {
    gamma[sample(d_z, valid_cnt)] <- 0
  }
  fn <- function(lambda) {
    var_mu <- as.numeric(t(gamma * lambda) %*% Sigma_z %*% (gamma * lambda)) +
      theta^2 * var_x + 2 * theta * as.numeric(beta %*% Sigma_z %*% (gamma * lambda))
    (var_mu / var_y - pve_y)^2
  }
  lambda <- optim(1, fn, method = 'Brent', lower = -10, upper = 10)$par
  gamma <- gamma * lambda
  var_mu <- as.numeric(t(gamma) %*% Sigma_z %*% (gamma)) + theta^2 * var_x + 
    2 * theta * as.numeric(beta %*% Sigma_z %*% (gamma))
  
  # Solve a quadratic equation for eta_y
  b <- 2 * theta * rho * eta_x
  eta_y <- (-b + sqrt(b^2 - 4 * (var_mu - var_y))) / 2
  y <- as.numeric(z %*% gamma) + x * theta + eps[, 2] * eta_y
  
  # Export results
  dat <- data.table('x' = x, 'y' = y, z)
  params <- list(
    'theta' = theta, 'beta' = beta, 'gamma' = gamma, 'eta_x' = eta_x, 
    'eta_y' = eta_y, 'rho' = rho, 'snr_x' = snr_x, 'snr_y' = snr_y
  )
  out <- list('dat' = dat, 'params' = params)
  return(out)
}

################################################################################

#' Sample data from the leaky IV model with var_x = var_y = 1
#' 
#' This function simulates the covariance matrix of a linear SEM with some 
#' leaky instrument(s).
#' 
#' @param d_z Dimensionality of Z.
#' @param z_rho Autocorrelation of the Toeplitz matrix for Z.
#' @param rho Correlation between residuals for X and Y.
#' @param theta Average treatment effect of X on Y.
#' @param snr_x Signal to noise ratio for X.
#' @param snr_y Signal to noise ratio for Y.
#' @param pr_valid Proportion of candidate instruments that are valid, 
#'   i.e. have no direct effect on Y.
#' 

sim_cov <- function(d_z, z_rho, rho, theta, snr_x, snr_y, pr_valid) {
  
  # What proportion are valid?
  valid_cnt <- round(pr_valid * d_z)
  pr_valid <- valid_cnt / d_z
  
  # Draw Z's
  var_z <- 1 / d_z
  Sigma_z <- toeplitz(z_rho^(0:(d_z - 1))) * var_z
  
  # Draw random weights for beta, calculate var_x, eta_x, Sigma_zx
  pve_x <- snr_x / (1 + snr_x)
  beta <- rnorm(d_z)
  var_x <- as.numeric(t(beta) %*% Sigma_z %*% beta) / pve_x 
  eta_x <- sqrt(var_x * (1 - pve_x))
  Sigma_xz <- beta %*% Sigma_z
  Sigma_zx <- t(Sigma_xz)
  
  # Draw random weights for gamma, tune magnitude to ensure target SNR
  var_y <- 10
  pve_y <- snr_y / (1 + snr_y)
  gamma <- rnorm(d_z)
  if (pr_valid > 0) {
    gamma[sample(d_z, valid_cnt)] <- 0
  }
  fn <- function(lambda) {
    var_mu <- as.numeric(t(gamma * lambda) %*% Sigma_z %*% (gamma * lambda)) +
      theta^2 * var_x + 2 * theta * as.numeric(beta %*% Sigma_z %*% (gamma * lambda))
    (var_mu / var_y - pve_y)^2
  }
  lambda <- optim(1, fn, method = 'Brent', lower = -10, upper = 10)$par
  gamma <- gamma * lambda
  var_mu <- as.numeric(t(gamma) %*% Sigma_z %*% (gamma)) + theta^2 * var_x + 
    2 * theta * as.numeric(beta %*% Sigma_z %*% (gamma))
  
  # Solve a quadratic equation for eta_y, calculate Y's covariance with X, Z
  b <- 2 * theta * rho * eta_x
  eta_y <- (-b + sqrt(b^2 - 4 * (var_mu - var_y))) / 2
  Sigma_yz <- gamma %*% Sigma_z + theta * Sigma_xz
  sigma_yx <- theta * var_x + tcrossprod(gamma, Sigma_xz) + rho * eta_x * eta_y
  
  # Format, export
  block1 <- matrix(c(var_x, sigma_yx, Sigma_zx, 
                     sigma_yx, var_y, Sigma_yz), ncol = 2)
  z_plus <- rbind(Sigma_xz, Sigma_yz, Sigma_z)
  Sigma <- cbind(block1, z_plus)
  colnames(Sigma) <- rownames(Sigma) <- c('x', 'y', paste0('z', seq_len(d_z)))
  params <- list(
    'theta' = theta, 'beta' = beta, 'gamma' = gamma, 'eta_x' = eta_x, 
    'eta_y' = eta_y, 'rho' = rho, 'snr_x' = snr_x, 'snr_y' = snr_y
  )
  out <- list('Sigma' = Sigma, 'params' = params)
  return(out)
  
}

################################################################################

#' Sample data from the leaky IV model with var_x = var_y = 1
#' 
#' This function simulates data from a linear SEM with some leaky instrument(s).
#' 
#' @param n Sample size.
#' @param d_z Dimensionality of Z.
#' @param z_rho Autocorrelation of the Toeplitz matrix for Z.
#' @param rho Correlation between residuals for X and Y.
#' @param theta Average treatment effect of X on Y.
#' @param r2_x Proportion of variance explained for X.
#' @param r2_y Proportion of variance explained for Y.
#' @param pr_valid Proportion of candidate instruments that are valid, 
#'   i.e. have no direct effect on Y.
#' 

sim_dat2 <- function(n, d_z, z_cnt, z_rho, rho, theta, r2_x, r2_y, pr_valid) {
  # What proportion are valid?
  valid_cnt <- round(pr_valid * d_z)
  pr_valid <- valid_cnt / d_z
  # Draw Z's
  var_z <- 1 / d_z
  z <- matrix(rnorm(n * d_z), ncol = d_z)
  Sigma_z <- toeplitz(z_rho^(0:(d_z - 1))) * var_z
  z <- z %*% chol(Sigma_z)
  colnames(z) <- paste0('z', seq_len(d_z))
  # Simulate standardized residual vectors
  #rho <- rho * sample(c(-1, 1), 1)
  Sigma_eps <- matrix(c(1, rho, rho, 1), ncol = 2)
  eps <- matrix(rnorm(n * 2), ncol = 2)
  eps <- eps %*% chol(Sigma_eps)
  # Draw random beta, set var_mu = r2_x
  beta <- runif(d_z, min = -1, 1)
  fn <- function(lambda) {
    var_mu <- as.numeric(t(beta * lambda) %*% Sigma_z %*% (beta * lambda))
    (var_mu - r2_x)^2
  }
  lambda <- optim(1, fn, method = 'Brent', lower = 0, upper = 10)$par
  beta <- beta * lambda
  var_x <- 1
  eta_x <- sqrt(var_x * (1 - r2_x))
  x <- as.numeric(z %*% beta) + eps[, 1] * eta_x
  # Set theta to account for a third, a half or two thirds of signal variance
  if (theta == 'low') {
    theta <- sqrt(r2_y / 3)
  } else if (theta == 'med') {
    theta <- sqrt(r2_y / 2) 
  } else if (theta == 'high') {
    theta <- sqrt(2 * r2_y / 3) 
  }
  #theta <- theta * sample(c(-1, 1), 1)
  # Draw random gamma, set var_mu = r2_y
  gamma <- runif(d_z, min = -1, 1) / (d_z * (1 - pr_valid))
  if (pr_valid > 0) {
    gamma[sample(d_z, valid_cnt)] <- 0
  }
  fn <- function(lambda) {
    var_mu <- as.numeric(t(gamma * lambda) %*% Sigma_z %*% (gamma * lambda)) +
      theta^2 * var_x + 2 * theta * as.numeric(beta %*% Sigma_z %*% (gamma * lambda))
    (var_mu - r2_y)^2
  }
  lambda <- optim(1, fn, method = 'Brent', lower = 0, upper = 10)$par
  gamma <- gamma * lambda
  # Solve a quadratic equation for eta_y
  var_mu <- r2_y
  var_y <- 1
  b <- 2 * theta * rho * eta_x
  eta_y <- (-b + sqrt(b^2 - 4 * (var_mu - var_y))) / 2
  y <- as.numeric(z %*% gamma) + x * theta + eps[, 2] * eta_y
  # Export results
  dat <- data.table('x' = x, 'y' = y, z)
  params <- list(
    'theta' = theta, 'beta' = beta, 'gamma' = gamma, 'eta_x' = eta_x, 
    'eta_y' = eta_y, 'rho' = rho, 'r2_x' = r2_x, 'r2_y' = r2_y
  )
  out <- list('dat' = dat, 'params' = params)
  return(out)
}



