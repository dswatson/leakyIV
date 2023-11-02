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
#' @param tau Upper bound on the L2 norm of gamma.
#' @param n_rho Number of rhos to evaluate.
#' @param n_boot Number of bootstrap replicates.
#' @param bayes Use Bayesian bootstrap?
#' @param parallel Compute bootstrap estimates in parallel?

soft_iv <- function(dat, tau, n_rho, n_boot, bayes, parallel) {
  n <- nrow(dat)
  p <- ncol(dat)
  d_z <- sum(grepl('z', colnames(dat)))
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
    aa <- as.numeric(Sigma_xz %*% Theta_z2 %*% Sigma_zx)
    bb <- as.numeric(Sigma_xz %*% Theta_z2 %*% Sigma_zy)
    cc <- as.numeric(Sigma_yz %*% Theta_z2 %*% Sigma_zx)
    dd <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
    ee <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
    ff <- as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
    # Find tau-feasible region
    lo_tau <- (bb - sqrt(aa * (tau - cc) + bb^2)) / aa
    hi_tau <- (bb + sqrt(aa * (tau - cc) + bb^2)) / aa
    # Find rho-feasible region
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
    lo_rho <- alpha_fn(0.9999)
    hi_rho <- alpha_fn(-0.9999)
    # Find the tightest bound
    alpha_min <- max(c(lo_tau, lo_rho))
    alpha_max <- min(c(hi_tau, hi_rho))
    # Compute corresponding values of rho
    rho_fn <- function(alpha) {
      RETURN_RHO
    }
    rho_min <- rho_fn(alpha_max)
    rho_max <- rho_fn(alpha_min)
    # Optionally return rho-alpha grid
    if (is.null(n_rho)) {
      out <- data.frame('alpha' = c(alpha_min, alpha_max), 
                        'rho' = c(rho_max, rho_min))
    } else {
      rhos <- seq(rho_max, rho_min, length.out = n_rho)
      out <- data.frame('alpha' = sapply(rhos, function(r) alpha_fn(r)), 
                        'rho' = rhos)
    }
    return(out)
  }
  # Run bootstrap
  if (is.null(n_boot)) {
    boot_out <- boot_loop(0)
    if(is.null(n_rho)) {
      out <- data.frame('bound' = c('lower', 'upper'),
                        'alpha' = boot_out$alpha, 
                        'se' = NA_real_,
                        'rho' = boot_out$rho)
    } else {
      
      
      
    }
    out <- data.frame('alpha' = alpha, 'se' = NA_real_, 
                      'tau_in' = tau, 'rho_in' = rhos)
    out <- out[!is.na(out$alpha), ]
  } else {
    if (parallel) {
      alpha <- foreach(i = seq_len(n_boot), .combine = cbind) %dopar% 
        boot_loop(i)
    } else {
      alpha <- sapply(seq_len(n_boot), function(i) boot_loop(i))
    }
    out <- data.frame('alpha' = rowMeans(alpha, na.rm = TRUE), 
                      'se' = rowSds(alpha, na.rm = TRUE), 
                      'tau_in' = tau, 'rho_in' = rhos)
    out <- na.omit(out)
  }
  return(out)
}

# Loop over different data configurations
loop_fn <- function(idx_b, alpha_b, z_rho_b, rho_b, r2_xb, r2_yb, prop_b) {
  tmp <- sim_dat(
    n = 1000, d_z = 4, z_cnt = TRUE, z_rho = z_rho_b,
    rho = rho_b, alpha = alpha_b, r2_x = r2_xb, r2_y = r2_yb, 
    pr_valid = prop_b, idx_b
  )
  # Assume oracle access to tau
  soft_iv(tmp$dat, sum(tmp$params$gamma^2), n_rho = 100, 
          n_boot = 200, bayes = FALSE, parallel = FALSE) %>%
    mutate('ACE' = alpha_b, 'z_rho' = z_rho_b, 'rho' = rho_b,  
           'r2_x' = r2_xb, 'r2_y' = r2_yb, 'pr_valid' = prop_b) %>%
    return(.)
}
# Execute in parallel
df <- foreach(alphas = c(-2, -1, 1, 2), .combine = rbind) %:%
  foreach(rhos = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75), .combine = rbind) %dopar%
  loop_fn(1, alphas, 0, rhos, 0.75, 0.75, 0)

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

# Can we write a rho_fn that takes alpha as input?



