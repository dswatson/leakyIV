# Set working directory
setwd('~/Documents/UCL/soft_instruments')

# Load libraries, register cores
library(data.table)
library(matrixStats)
library(tidyverse)
library(doMC)
registerDoMC(8)

# Simulation script
source('simulation.R')

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

soft_iv <- function(params, n_rho) {
  # Assume diagonal covariance for Z
  d_z <- length(params$beta)
  var_z <- 1 / d_z
  Sigma_z <- diag(var_z, d_z)
  # Extract all other parameters
  alpha <- params$alpha
  beta <- params$beta
  gamma <- params$gamma
  rho <- params$rho
  eta_x <- params$eta_x
  eta_y <- params$eta_y
  var_x <- var_y <- 1
  # Data covariance
  Theta_z <- solve(Sigma_z)
  Theta_z2 <- Theta_z %*% Theta_z
  Sigma_zx <- Sigma_z %*% beta
  Sigma_xz <- t(Sigma_zx)
  Sigma_zy <- Sigma_z %*% gamma + Sigma_zx %*% alpha
  Sigma_yz <- t(Sigma_zy)
  sigma_xy <- Sigma_xz %*% gamma + alpha * var_x + eta_x * eta_y * rho
  # Define
  aa <- as.numeric(Sigma_xz %*% Theta_z2 %*% Sigma_zx)
  bb <- as.numeric(Sigma_xz %*% Theta_z2 %*% Sigma_zy)
  cc <- as.numeric(Sigma_yz %*% Theta_z2 %*% Sigma_zx)
  dd <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
  ee <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
  ff <- as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
  # Find tau-feasible region
  l2 <- sum(gamma^2)
  qd_bnd <- cc - (bb^2 / aa)
  tau <- ifelse(l2 > qd_bnd, l2, qd_bnd)
  lo <- (bb - sqrt(aa * (tau - cc) + bb^2)) / aa
  hi <- (bb + sqrt(aa * (tau - cc) + bb^2)) / aa
  # Compute alpha as function of rho
  alpha_fn <- function(rho) {
    gg <- (ee - var_x) * (1 + ((ee - var_x) / (eta_x^2 * rho^2)))
    hh <- -(dd - sigma_xy) * (1 + ((ee - var_x) / (eta_x^2 * rho^2)))
    ii <- ((dd - sigma_xy)^2 / (eta_x^2 * rho^2)) + ff - var_y
    if (rho < 0) {
      alpha <- (-hh + sqrt(hh^2 - gg * ii)) / gg 
    } else {
      alpha <- (-hh - sqrt(hh^2 - gg * ii)) / gg 
    }
    return(alpha)
  }
  rhos <- seq(-0.99, 0.99, length.out = n_rho)
  alphas <- sapply(rhos, function(r) alpha_fn(r))
  alphas <- ifelse(alphas >= lo & alphas <= hi, alphas, NA_real_)
  # Export
  out <- data.frame('rho_in' = rhos, 'alpha' = alphas)
  out <- na.omit(out)
  return(out)
}
loop_fn <- function(idx_b, alpha_b, rho_b, r2_xb, r2_yb, prop_b) {
  tmp <- sim_dat(
    n = 1000, d_z = 4, z_cnt = TRUE, rho = rho_b, alpha = alpha_b, 
    r2_x = r2_xb, r2_y = r2_yb, pr_valid = prop_b, idx_b
  )
  soft_iv(tmp$params, n_rho = 200) %>%
    mutate('ACE' = alpha_b, 'rho' = rho_b, 'r2_x' = r2_xb, 'r2_y' = r2_yb, 
           'pr_valid' = prop_b) %>%
    return(.)
}
# Execute in parallel
df <- foreach(alphas = c(-1, -0.5, 0.5, 1), .combine = rbind) %:%
  foreach(rhos = c(-0.75, -0.25, 0.25, 0.75), .combine = rbind) %dopar%
  loop_fn(1, alphas, rhos, 0.5, 0.5, 0)
# Plot results
df %>%
  ggplot(aes(rho_in, alpha)) + 
  geom_line() + 
  geom_hline(aes(yintercept = ACE), color = 'red', linetype = 'dashed') + 
  geom_vline(aes(xintercept = rho), color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 0, color = 'blue') + 
  theme_bw() + 
  facet_grid(ACE ~ rho, scales = 'free_y') + 
  labs(x = 'Unobserved Confounding', y = 'Average Causal Effect')


fn <- function(tau) {
  lo <- (bb - sqrt(aa * (tau - cc) + bb^2)) / aa
  hi <- (bb + sqrt(aa * (tau - cc) + bb^2)) / aa
  out <- data.frame('tau' = tau, 'lo' = lo, 'hi' = hi)
  return(out)
}
min_tau <- sum(gamma^2) - 0.05
df <- foreach(taus = seq(min_tau, 4, length.out = 500), .combine = rbind) %dopar%
  fn(taus)
df <- na.omit(df)
ggplot(df) +
  geom_ribbon(aes(tau, ymin = lo, ymax =  hi), alpha = 0.5) + 
  theme_bw() + 
  labs(y = 'alpha')













