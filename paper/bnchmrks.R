# Set working directory
setwd('~/Documents/Kings/leakyIV/paper')

# Load libraries, register cores
library(data.table)
library(broom)
library(sisVIVE)
library(leakyIV)
library(ggplot2)
library(ggsci)
library(doMC)
registerDoMC(8)

# Load scripts
source('bayesian_baseline.R')
source('MBE.R')
source('simulator.R')


# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Benchmark against 2SLS, sisVIVE, MBE, and Bayes:
bnchmrk <- function(z_rho, rho, pr_valid, n = 1000, n_sim = 100) {
  
  # Generate data, extract "population" data
  sim <- sim_dat(n = 1e5, d_z = 20, z_cnt = TRUE, z_rho, rho, 
                 theta = 1, r2_x = 3/4, r2_y = 3/4, pr_valid)
  d_z <- 20
  d <- d_z + 2
  l2 <- sqrt(sum(sim$params$gamma^2))
  tau <- 1.5 * l2
  
  # Inner loop
  inner_loop <- function(b) {
    
    # Draw data
    tmp <- sim$dat[sample(1e5, n), ]
    
    # Run 2SLS
    f0 <- lm(x ~ ., data = tmp[, -c('y')])
    tmp[, x_hat := fitted(f0)]
    f1 <- lm(y ~ x_hat, data = tmp)
    ate_2sls <- as.numeric(tidy(f1)$estimate[2])
    
    # Run sisVIVE
    sisvive <- cv.sisVIVE(tmp$y, tmp$x, as.matrix(tmp[, -c('x', 'y')]))
    ate_sisvive <- sisvive$beta
    
    # Run MBE
    beta_hat <- as.numeric(tidy(f0)$estimate[2:(d_z + 1)])
    se_beta <- as.numeric(tidy(f0)$std.error[2:(d_z + 1)])
    f2 <- lm(y ~ ., data = tmp[, -c('x_hat', 'x')])
    gamma_hat <- tidy(f2)$estimate[2:(d_z + 1)]
    se_gamma <- tidy(f2)$std.error[2:(d_z + 1)]
    mbe <- MBE(beta_hat, gamma_hat, se_beta, se_gamma, phi = 1, n_boot = 1)
    ate_mbe <- mbe$Estimate[2]
    
    # # Bayesian baseline
    # bb_res <- baseline_gaussian_sample(
    #   m = 2000, dat = as.matrix(tmp), pos_z = 1:d_z, pos_x = d_z + 1, 
    #   pos_y = d_z + 1, prop_var = 0.01, lvar_z_mu = -1, lvar_z_var = 3,
    #   beta_var = 5, pre_gamma_var = 1, theta_var = 5, lvar_error_mu = -1,
    #   lvar_error_var = 1, tau = tau, alpha = 0.05
    # )
    # ate_bb <- mean(bb_res$theta)
    # ate_lo_bb <- bb_res$lb_interval[2]
    # ate_hi_bb <- bb_res$ub_interval[1]
    
    # LeakyIV
    suppressWarnings(
      ate_bnds <- leakyIV(tmp$x, tmp$y, tmp[, -c('x', 'y')], tau = tau,
                          method = 'glasso', rho = 0.01)
    )
    ate_lo <- ate_bnds$ATE_lo
    ate_hi <- ate_bnds$ATE_hi
    
    # Export
    out <- data.table(b,
      method = c('TSLS', 'sisVIVE', 'MBE', 'Lower', 'Upper'),
      theta = c(ate_2sls, ate_sisvive, ate_mbe, ate_lo, ate_hi)
    )
    return(out)
  }
  out <- foreach(bb = seq_len(n_sim), .combine = rbind) %do% inner_loop(bb)
  
  # Export
  out[, z_rho := z_rho][, rho := rho][, pr_valid := pr_valid]
  return(out)
  
}

# Execute in parallel
df <- foreach(zz = c(0, 0.5), .combine = rbind) %:%
  foreach(rr = c(1/4, 3/4), .combine = rbind) %:%
  foreach(pp = seq(0, 0.95, by = 0.05), .combine = rbind) %dopar%
  bnchmrk(zz, rr, pp)

################################################################################

# Plot it
df[, mu := mean(theta, na.rm = TRUE), by = .(method, z_rho, rho, pr_valid)]
df[, se := sd(theta, na.rm = TRUE), by = .(method, z_rho, rho, pr_valid)]
tmp <- unique(df[, .(method, mu, se, z_rho, rho, pr_valid)])
tmp[, z_rho := fifelse(z_rho == 0, 'Diagonal', 'Toeplitz')]
tmp[, rho := fifelse(rho == 0.75, 'Strong Confounding', 'Weak Confounding')]
#tmp[method == 'TSLS', method := '2SLS']
#tmp[, method := factor(method, levels = c(''))]

p1 <- ggplot(tmp, aes(pr_valid, mu, color = method, fill = method)) + 
  geom_ribbon(aes(ymin = mu + se, ymax = mu - se), alpha = 0.4) + 
  geom_line() + 
  geom_hline(yintercept = 1, linewidth = 1, color = 'grey60') +
  scale_color_d3() +
  scale_fill_d3() +
  labs(x = 'Proportion of Valid Instruments',
       y = expression(paste('Average Treatment Effect ', theta))) +
  facet_grid(z_rho ~ rho) +
  theme_bw() + 
  theme(legend.position = 'bottom')



theme(axis.title = element_text(size = 24),
      legend.title = element_text(size = 24),
      legend.text = element_text(size = 24),
      legend.position = 'bottom')
  























