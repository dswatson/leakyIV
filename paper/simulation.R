#' Simulation script for Leaky IV project
# Set working directory
setwd('~/Documents/Kings/leakyIV')

# Load libraries, register cores
library(data.table)
library(doMC)
registerDoMC(8)

# Set seed
set.seed(123)

#' Sample data from the leaky IV model
#' 
#' This function simulates data from a linear structural equation model with 
#' some leaky instrument(s).
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
#' @param s_idx Simulation index.
#' 

sim_dat <- function(n, d_z, z_rho, rho, theta, r2_x, r2_y, pr_valid, sim_idx) {
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
  rho <- rho * sample(c(-1, 1), 1)
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
  theta <- theta * sample(c(-1, 1), 1)
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
  dat <- data.table(z, x, y)
  fwrite(dat, paste0('./simulations/', sim_idx, '.csv'))
}

# Create simulation grid
sim_grd <- as.data.table(expand.grid(
  'd_z' = c(4, 20, 100),
  'z_rho' = c(0, 1/2),
  'rho' = c(1/4, 3/4),
  'theta' = c('low', 'med', 'high'),
  'pr_valid' = c(0, 1/2)
))
sim_grd[, sim_idx := .I]
setcolorder(sim_grd, 'sim_idx')
fwrite(sim_grd, './simulations/simulation_grid.csv')

# Sample a dataset of n=10k for each setting, to later be partitioned into 
# 10 datasets of n=1k
sapply(seq_len(nrow(sim_grd)), function(i) {
  tmp <- sim_grd[i, ]
  sim_dat(n = 1e4L, d_z = tmp$d_z, z_rho = tmp$z_rho, rho = tmp$rho,
          theta = tmp$theta, r2_x = 0.8, r2_y = 0.8, pr_valid = tmp$pr_valid,
          sim_idx = i)
})



