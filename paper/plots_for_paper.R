# Set working directory
setwd('~/Documents/Kings/leakyIV/paper')

# Load libraries
library(data.table)
library(ggplot2)
library(cowplot)
library(ggh4x)
library(doMC)
registerDoMC(8)
source('simulator.R')

# Set seed
set.seed(123)

# Data generating function

fn <- function(n, d_z, z_cnt, z_rho, rho, theta, r2_x, r2_y, pr_valid, n_boot) {
  
  # Generate data
  sim <- sim_dat(n, d_z, z_cnt, z_rho, rho, theta, r2_x, r2_y, pr_valid)
  dat <- sim$dat
  n <- nrow(dat)
  d <- ncol(dat)
  d_z <- ncol(dat) - 2
  if (n_boot > 0) {
    b_mat <- matrix(sample(n, n * n_boot, replace = TRUE), ncol = n_boot)
  }
  
  # Boot loop
  loop <- function(boot_idx) {
    
    # Draw bootstrap sample
    if (boot_idx > 0) {
      dat <- dat[b_mat[, boot_idx], ]
    }
    
    # Covariance parameters
    Sigma <- cov(dat)
    Theta_z <- solve(Sigma[seq_len(d_z), seq_len(d_z)])
    Sigma_z <- Sigma[seq_len(d_z), seq_len(d_z)]
    Sigma_zy <- matrix(Sigma[seq_len(d_z), d], ncol = 1L)
    Sigma_yz <- t(Sigma_zy)
    Sigma_zx <- matrix(Sigma[seq_len(d_z), d - 1L], ncol = 1L)
    Sigma_xz <- t(Sigma_zx)
    sigma_xy <- Sigma[d, d - 1L]
    var_x <- Sigma[d - 1L, d - 1L]
    var_y <- Sigma[d, d]
    
    # Our magic variables
    eta_x2 <- var_x - as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
    phi2 <- var_y - as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
    psi <- sigma_xy - as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
    
    # Compute theta as a function of rho, norm as a function of theta
    theta_fn <- function(rho) {
      theta <- (psi / eta_x2) - sign(rho) * 
        ((sqrt((1 - 1/rho^2) * (psi^2 - phi2 * eta_x2)))) / 
        (-eta_x2 * (1 - 1/rho^2))
      return(theta)
    }
    norm_fn <- function(theta, p) {
      gamma <- as.numeric(Theta_z %*% (Sigma_zy - theta * Sigma_zx))
      norm <- (sum(abs(gamma)^p))^(1 / p)
      return(norm)
    }
    
    # Loop through a bunch of rho-values
    df <- data.table(rho = seq(-0.99, 0.99, length.out = 1000))
    df[, theta := sapply(rho, theta_fn)]
    df[, l1 := sapply(theta, norm_fn, 1)]
    df[, l2 := sapply(theta, norm_fn, 2)]
    
    # Export
    df[, b := boot_idx]
    setcolorder(df, 'b')
    return(df)
  }
  
  # Compute in parallel, export
  if (n_boot > 0) {
    df <- foreach(bb = seq_len(n_boot), .combine = rbind) %dopar% loop(bb)
  } else {
    df <- loop(0)
  }
  return(df)
  
}

################################################################################

# Plot the rho-theta curve with bootstraps
df <- fn(n = 1000, d_z = 4, z_cnt = TRUE, z_rho = 1/2, 
         rho = 1/4, theta = 1, r2_x = 2/3, r2_y = 2/3, 
         pr_valid = 0, n_boot = 1000)
df[, theta_mu := mean(theta), by = rho]
df[, theta_se := sd(theta), by = rho]
df[, l2_mu := mean(l2), by = rho]
df[, l2_se := sd(l2), by = rho]
df <- unique(df[, .(rho, theta_mu, theta_se, l2_mu, l2_se)])
p1 <- ggplot(df, aes(rho, theta_mu)) + 
  geom_ribbon(aes(ymin = theta_mu - theta_se * qnorm(0.95),
                  ymax = theta_mu + theta_se * qnorm(0.95)),
              color = 'grey90', alpha = 0.25) +
  geom_line(linewidth = 0.75) + 
  labs(x = expression(paste('Confounding Coefficient ', rho)), 
       y = expression(paste('Average Treatment Effect ', theta))) +
  theme_bw() +
  theme(axis.title = element_text(size = 22))

# Plot the theta-L2 curve with bootstraps
p2 <- ggplot(df, aes(theta_mu, l2_mu)) + 
  geom_ribbon(aes(ymin = l2_mu - l2_se * qnorm(0.95),
                  ymax = l2_mu + l2_se * qnorm(0.95)),
              color = 'grey90', alpha = 0.25) +
  geom_line(linewidth = 0.75) + 
  labs(x = expression(paste('Average Treatment Effect ', theta)),
       y = expression(paste(L[2], ' Norm, Leakage Weights ', gamma))) +
  theme_bw() +
  theme(axis.title = element_text(size = 22))

plot_grid(p1, p2, labels = c('A', 'B'), label_size = 22)
ggsave2('./plots/lemmas.pdf', height = 7, width = 14)
















