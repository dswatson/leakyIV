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

# Benchmark against backdoor adjustment, 2SLS, sisVIVE, and MBE:
bnchmrk <- function(z_rho, rho, pr_valid, n, n_sim) {
  
  
  # Generate data, extract "population" data
  sim <- sim_dat(n = 1e5, d_z = 20, z_cnt = TRUE, z_rho, rho,
                 theta = 1, r2_x = 1/4, r2_y = 3/4, pr_valid)
  # sim <- sim_dat2(n = 1e5, d_z = 20, z_cnt = TRUE, z_rho, rho,
  #                 theta = 'high', r2_x = 4/5, r2_y = 4/5, pr_valid)
  d_z <- 20
  d <- d_z + 2
  l2 <- sqrt(sum(sim$params$gamma^2))
  tau <- 1.2 * l2
  
  # Inner loop
  inner_loop <- function(b) {
    
    # Draw data
    tmp <- sim$dat[sample(1e5, n), ]
    
    # Run backdoor adjustment
    f0 <- lm(y ~ ., data = tmp)
    ate_bda <- as.numeric(tidy(f0)$estimate[d_z + 2])
    
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
                          method = 'glasso', rho = 0.005)
      # ate_bnds <- leakyIV(tmp$x, tmp$y, tmp[, -c('x', 'y')], tau = tau,
      #                     method = 'mle')
    )
    ate_lo <- ate_bnds$ATE_lo
    ate_hi <- ate_bnds$ATE_hi
    
    # Export
    out <- data.table(b,
      method = c('Backdoor', 'TSLS', 'sisVIVE', 'MBE', 'Lower', 'Upper'),
      theta = c(ate_bda, ate_2sls, ate_sisvive, ate_mbe, ate_lo, ate_hi)
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
  foreach(rr = c(0.25, 0.75), .combine = rbind) %:%
  foreach(pp = seq(0, 0.95, by = 0.05), .combine = rbind) %dopar%
  #foreach(pp = seq(0, 0.9, by = 0.1), .combine = rbind) %dopar%
  bnchmrk(zz, rr, pp, n = 1000, n_sim = 20)

################################################################################

# Spot check
df[method %in% c('Lower', 'Upper'), sum(is.na(theta)), by = .(z_rho, rho, pr_valid)]

# Plot it
df[, mu := mean(theta, na.rm = TRUE), by = .(method, z_rho, rho, pr_valid)]
df[, se := sd(theta, na.rm = TRUE), by = .(method, z_rho, rho, pr_valid)]
tmp <- unique(df[, .(method, mu, se, z_rho, rho, pr_valid)])
tmp[, z_rho := fifelse(z_rho == 0, 'Diagonal', 'Toeplitz')]
tmp[, rho := fifelse(rho > 0.5, 'Strong Confounding', 'Weak Confounding')]
tmp[, rho := factor(rho, levels = c('Weak Confounding', 'Strong Confounding'))]
tmp[method == 'TSLS', method := '2SLS']
tmp[method == 'Lower', method := 'LeakyIV_lo']
tmp[method == 'Upper', method := 'LeakyIV_hi']
tmp[, method := factor(
  method, levels = c('Backdoor', '2SLS', 'sisVIVE', 'MBE', 'LeakyIV_lo', 'LeakyIV_hi')
)]
setnames(tmp, 'method', 'Method')

p1 <- ggplot(tmp, aes(pr_valid, mu, fill = Method)) + 
  geom_ribbon(aes(ymin = mu + se, ymax = mu - se), alpha = 0.4) + 
  geom_line(aes(color = Method)) + 
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
  







################################################################################

# Benchmark against Bayesian baseline

# Initialize
out <- data.table(
  b = NA_real_, method = NA_character_, theta = NA_real_
)
saveRDS(out, './bnchmrks/bayesian.rds')

bayes_bnchmrk <- function(z_rho, rho, pr_valid, n = 1000, n_sim = 100) {
  
  # Generate data, extract "population" data
  sim <- sim_dat(n = 1e5, d_z = 20, z_cnt = TRUE, z_rho, rho,
                 theta = 1, r2_x = 3/4, r2_y = 3/4, pr_valid)
  d_z <- 20
  d <- d_z + 2
  l2 <- sqrt(sum(sim$params$gamma^2))
  tau <- 1.1 * l2
  
  # Inner loop
  inner_loop <- function(b) {
    
    # Draw data
    tmp <- sim$dat[sample(1e5, n), ]
    
    # Bayesian baseline
    bb_res <- baseline_gaussian_sample(
      m = 2000, dat = as.matrix(tmp), pos_z = 1:d_z, pos_x = d_z + 1,
      pos_y = d_z + 1, prop_var = 0.01, lvar_z_mu = -1, lvar_z_var = 3,
      beta_var = 5, pre_gamma_var = sqrt(tau) / d_z, theta_var = 5, 
      lvar_error_mu = -1, lvar_error_var = 1, tau = tau, alpha = 0.05
    )
    ate_bb <- mean(bb_res$theta)
    ate_lo_bb <- bb_res$lb_interval[2]
    ate_hi_bb <- bb_res$ub_interval[1]
    
    # LeakyIV
    suppressWarnings(
      ate_bnds <- leakyIV(tmp$x, tmp$y, tmp[, -c('x', 'y')], tau = tau,
                          method = 'glasso', rho = 0.005)
    )
    ate_lo <- ate_bnds$ATE_lo
    ate_hi <- ate_bnds$ATE_hi
    
    # Export
    out <- data.table(b, 
                      method = c('Bayes_mean', 'Bayes_lo', 'Bayes_hi',
                                 'leakyIV_lo', 'leakyIV_hi'),
                      theta = c(ate_bb, ate_lo_bb, ate_hi_bb, ate_lo, ate_hi)
    )
    return(out)
  }
  out <- foreach(bb = seq_len(n_sim), .combine = rbind) %do% inner_loop(bb)
  
  # Import, export
  old <- readRDS('./bnchmrks/bayesian.rds')
  new <- rbind(old, out)
  saveRDS(new, './bnchmrks/bayesian.rds')
}

# Execute in parallel over a grid
foreach(zz = c(0, 0.5), .combine = rbind) %:%
  foreach(rr = c(1/4, 3/4), .combine = rbind) %:%
  foreach(pp = seq(0, 0.95, by = 0.05), .combine = rbind) %dopar%
  bnchmrk(zz, rr, pp, n_sim = 10)




# Weak confounding is hard because it means that tau is near its theoretical 
# minimum. We actually do better witih more confounding, weirdly









