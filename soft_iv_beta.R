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

# Define soft_iv function
#' @param dat Input dataset with instruments z, treatment x, and outcome y.
#' @param p Power of the gamma norm.
#' @param tau Upper bound on the gamma norm.
#' @param n_rho Number of rhos to evaluate.
#' @param n_boot Number of bootstrap replicates.
#' @param bayes Use Bayesian bootstrap?
#' @param parallel Compute bootstrap estimates in parallel?

soft_iv <- function(dat, p, tau, n_rho, n_boot, bayes, parallel) {
  # Prelimz
  n <- nrow(dat)
  p <- ncol(dat)
  d_z <- sum(grepl('z', colnames(dat)))
  rhos <- seq(-0.99, 0.99, length.out = n_rho)
  # Define bootstrap loop
  boot_loop <- function(b) {
    if (b == 0) {
      Sigma <- cov(dat)
    } else {
      if (bayes) {
        # Draw Dirichlet weights
        wts <- rexp(n)
        wts <- (wts / sum(wts)) * n
        # Estimate eta_x
        f1 <- lm(x ~ ., data = select(dat, -y), weights = wts)
        eta_x <- sqrt(weighted.mean(x = residuals(f1)^2, w = wts))
        # Estimate data covariance
        Sigma <- cov.wt(dat, wt = wts)$cov
      } else {
        # Draw bootstrap sample
        tmp <- dat[sample.int(n, replace = TRUE)]
        # Estimate eta_x
        f1 <- lm(x ~ ., data = select(tmp, -y))
        eta_x <- sqrt(mean((residuals(f1)^2)))
        # Estimate data covariance
        Sigma <- cov(tmp)
      }
    }
    # Extract elements of covariance matrix
    Sigma_z <- Sigma[seq_len(d_z), seq_len(d_z)]
    Theta_z <- solve(Sigma_z)
    Theta_z2 <- Theta_z %*% Theta_z
    Sigma_zy <- matrix(Sigma[seq_len(d_z), p], ncol = 1)
    Sigma_yz <- t(Sigma_zy)
    Sigma_zx <- matrix(Sigma[seq_len(d_z), p - 1], ncol = 1)
    Sigma_xz <- t(Sigma_zx)
    sigma_xy <- Sigma[p, p - 1]
    var_x <- Sigma[p - 1, p - 1]
    var_y <- Sigma[p, p]
    # Define
    dd <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
    ee <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
    ff <- as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
    # Compute alpha and gamma as functions of rho
    from_rho <- function(rho_in) {
      gg <- (ee - var_x) * (1 + ((ee - var_x) / (eta_x^2 * rho_in^2)))
      hh <- -(dd - sigma_xy) * (1 + ((ee - var_x) / (eta_x^2 * rho_in^2)))
      ii <- ((dd - sigma_xy)^2 / (eta_x^2 * rho_in^2)) + ff - var_y
      if (rho_in < 0) {
        alpha_hat <- as.numeric((-hh + sqrt(hh^2 - gg * ii)) / gg)
      } else {
        alpha_hat <- as.numeric((-hh - sqrt(hh^2 - gg * ii)) / gg)
      }
      gamma_hat <- as.numeric(Theta_z %*% (Sigma_zy - alpha_hat * Sigma_zx))
      # Take the norm
      norm <- (sum(abs(gamma_hat)^p))^(1 / p)
      out <- data.frame(
        'rho' = rho_in, 'alpha' = alpha_hat, 'norm' = norm
      )
      return(out)
    }
    rho_loop <- foreach(r = rhos, .combine = rbind) %do% from_rho(r) %>%
      filter(norm <= tau) %>%
      na.omit(.)
    out <- data.table(
      'b' = b,
      'bound' = c('lo', 'hi'), 
      'alpha' = c(min(rho_loop$alpha), max(rho_loop$alpha))
    )
    return(out)
  }
  # Run bootstrap
  if (parallel) {
    boots <- foreach(i = seq_len(n_boot), .combine = rbind) %dopar% 
      boot_loop(i)
  } else {
    boots <- foreach(i = seq_len(n_boot), .combine = rbind) %do% 
      boot_loop(i)
  }
  out <- data.frame(
    'bound' = c('lo', 'hi'),
    'alpha' = c(boots[bound == 'lo', mean(alpha)], 
                boots[bound == 'hi', mean(alpha)]),
    'se' = c(boots[bound == 'lo', sd(alpha)],
             boots[bound == 'hi', sd(alpha)])
  )
  return(out)
}


# Loop over different data configurations
loop_fn <- function(idx_b, alpha_b, rho_b, r2_xb, r2_yb, prop_b) {
  tmp <- sim_dat(
    n = 1000, d_z = 4, z_cnt = TRUE, rho = rho_b, alpha = alpha_b, 
    r2_x = r2_xb, r2_y = r2_yb, pr_valid = prop_b, idx_b
  )
  # Assume oracle access to tau
  soft_iv(tmp$dat, sum(tmp$params$gamma^2), n_rho = 100, 
          n_boot = 200, bayes = TRUE, parallel = FALSE) %>%
    mutate('ACE' = alpha_b, 'rho' = rho_b,  'r2_x' = r2_xb, 'r2_y' = r2_yb, 
           'pr_valid' = prop_b) %>%
    return(.)
}
# Execute in parallel
df <- foreach(alphas = c(-0.5, -0.25, 0.25, 0.5), .combine = rbind) %:%
  foreach(rhos = c(-2/3, -1/3, 1/3, 2/3), .combine = rbind) %dopar%
  loop_fn(1, alphas, rhos, 0.5, 0.5, 0)

# Plot results
df %>%
  ggplot(aes(rho_in, alpha)) + 
  geom_ribbon(aes(ymin = alpha - qnorm(0.975) * se, 
                  ymax = alpha + qnorm(0.975) * se), alpha = 0.2) +
  geom_line(size = 0.25) + 
  geom_hline(aes(yintercept = ACE), color = 'red', linetype = 'dashed') + 
  geom_vline(aes(xintercept = rho), color = 'red', linetype = 'dashed') +
  geom_hline(yintercept = 0, color = 'blue') + 
  theme_bw() + 
  facet_grid(ACE ~ rho, scales = 'free_y') + 
  labs(x = 'Unobserved Confounding', y = 'Average Causal Effect')


