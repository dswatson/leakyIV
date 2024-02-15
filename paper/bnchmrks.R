# Set working directory
setwd('~/Documents/Kings/leakyIV/paper')

# Load libraries, register cores
library(data.table)
library(mvnfast)
library(broom)
library(sisVIVE)
library(leakyIV)
library(ggplot2)
library(ggsci)
library(cowplot)
library(doParallel)
registerDoParallel(8)

# Load scripts
source('bayesian_baseline.R')
source('MBE.R')
source('simulator.R')

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Benchmark against backdoor adjustment, 2SLS, sisVIVE, and MBE:
bnchmrk <- function(d_z, z_rho, rho, snr_x, snr_y, pr_valid, n, n_sim) {
  
  # Update
  cat(paste('Running:', d_z, z_rho, rho, snr_x, snr_y, '\n'))
  
  # Simulate population covariance matrix
  sim <- sim_cov(d_z, z_rho, rho, theta = 1, snr_x, snr_y, pr_valid)
  Sigma <- sim$Sigma
  l2 <- sqrt(sum(sim$params$gamma^2))
  d <- d_z + 2
  tau <- 1.1 * l2
  
  # Treat this as ground truth
  dat <- as.data.table(rmvn(n * n_sim, mu = rep(0, d), Sigma))
  colnames(dat) <- c('x', 'y', paste0('z', seq_len(d_z)))
  
  # Inner loop
  inner_loop <- function(b) {
    
    # Draw data
    tmp <- dat[((b - 1) * n + 1):(b * n), ]
    
    # Backdoor adjustment
    f0 <- lm(y ~ ., data = tmp)
    tidy_f <- as.data.table(tidy(f0))
    ate_bda <- tidy_f[term == 'x', estimate]
    
    # 2SLS
    f0 <- lm(x ~ ., data = tmp[, -c('y')])
    x_hat <- fitted(f0)
    f1 <- lm(tmp$y ~ x_hat)
    tidy_f1 <- as.data.table(tidy(f1))
    ate_2sls <- tidy_f1[term == 'x_hat', estimate]
    
    # sisVIVE
    sisvive <- cv.sisVIVE(tmp$y, tmp$x, as.matrix(tmp[, -c('x', 'y')]))
    ate_sisvive <- sisvive$beta
    
    # MBE
    tidy_f0 <- as.data.table(tidy(f0))
    beta_hat <- tidy_f0[grepl('z', term), estimate]
    se_beta <- tidy_f0[grepl('z', term), std.error]
    f1 <- lm(y ~ ., data = tmp[, -c('x')])
    tidy_f1 <- as.data.table(tidy(f1))
    gamma_hat <- tidy_f1[grepl('z', term), estimate]
    se_gamma <- tidy_f1[grepl('z', term), std.error]
    mbe <- MBE(beta_hat, gamma_hat, se_beta, se_gamma, phi = 1, n_boot = 1)
    ate_mbe <- mbe$Estimate[2]
    
    # Bayesian baseline
    bb_res <- baseline_gaussian_sample(
      m = 2000, dat = as.matrix(tmp), pos_z = 3:d, pos_x = 1, pos_y = 2, 
      prop_var = 0.01, lvar_z_mu = -1, lvar_z_var = 3,
      beta_var = 5, pre_gamma_var = 1, theta_var = 5, 
      lvar_error_mu = -1, lvar_error_var = 1, tau = tau, alpha = 0.025
    )
    ate_bb_lo <- bb_res$lb_interval[2]
    ate_bb_hi <- bb_res$ub_interval[1]
    
    # LeakyIV
    suppressWarnings(
      ate_leaky <- leakyIV(tmp, tau = tau, method = 'mle', parallel = FALSE)
    )
    ate_leaky_lo <- ate_leaky$ATE_lo
    ate_leaky_hi <- ate_leaky$ATE_hi
    
    # Export
    out <- data.table(b,
      method = c('Backdoor', 'TSLS', 'sisVIVE', 'MBE', 'Bayes_lo', 'Bayes_hi', 
                 'leaky_lo', 'leaky_hi'),
      theta = c(ate_bda, ate_2sls, ate_sisvive, ate_mbe, ate_bb_lo, ate_bb_hi,
                ate_leaky_lo, ate_leaky_hi)
    )
    return(out)
  }
  out <- foreach(bb = seq_len(n_sim), .combine = rbind) %do% inner_loop(bb)
  
  # Export
  out[, d_z := d_z][, z_rho := z_rho][, rho := rho][, 
    snr_x := snr_x][, snr_y := snr_y][, pr_valid := pr_valid]
  return(out)
  
}

# Execute in parallel
df <- foreach(dd = c(5, 10), .combine = rbind) %:%
  foreach(zz = c(0, 0.5), .combine = rbind) %:%
  foreach(rr = seq(-0.9, 0.9, 0.1), .combine = rbind) %:%
  foreach(s_x = c(1/2, 1, 2), .combine = rbind) %:%
  foreach(s_y = c(1/2, 1, 2), .combine = rbind) %dopar%
  bnchmrk(dd, zz, rr, s_x, s_y, pr_valid = 1/5, n = 1000, n_sim = 50)
saveRDS(df, './bnchmrks/big_sim.rds')

################################################################################

# Spot check
df[method %in% c('Lower', 'Upper'), sum(is.na(theta)), 
   by = .(d_z, z_rho, rho, snr_x, snr_y)]

# Prepare for plotting
df[, mu := mean(theta, na.rm = TRUE), by = .(method, d_z, z_rho, rho, snr_x, snr_y)]
df[, se := sd(theta, na.rm = TRUE), by = .(method, d_z, z_rho, rho, snr_x, snr_y)]
tmp <- unique(df[, .(method, mu, se, d_z, z_rho, rho, snr_x, snr_y)])
tmp[, z_rho := fifelse(z_rho == 0, 'Diagonal', 'Toeplitz')]
tmp[snr_y == 1/2, SNR := 'SNR = 1/2']
tmp[snr_y == 1, SNR := 'SNR = 1']
tmp[snr_y == 2, SNR := 'SNR = 2']
tmp[, SNR := factor(SNR, levels = c('SNR = 1/2', 'SNR = 1', 'SNR = 2'))]
tmp[method == 'TSLS', method := '2SLS']
setnames(tmp, 'method', 'Method')
tmp[, lo := .SD[Method == 'Lower', mu], by = .(d_z, z_rho, rho, snr_x, snr_y)]
tmp[, hi := .SD[Method == 'Upper', mu], by = .(d_z, z_rho, rho, snr_x, snr_y)]
tmp2 <- tmp[!Method %in% c('Lower', 'Upper')]
tmp2[, Method := factor(Method, levels = c('Backdoor', '2SLS', 'sisVIVE', 'MBE'))]

# We facet d_z rows and SNR_y columns
# Two values of Sigma_z x three values of SNR_x = 6 plots
plot_fn <- function(z, s) {
  p1 <- ggplot(tmp2[z_rho == z & snr_x == s], aes(rho, mu, fill = Method)) + 
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = 'LeakyIV'), alpha = 0.25) +
    geom_ribbon(aes(ymin = mu - se, ymax = mu + se), alpha = 0.5) + 
    geom_line(aes(color = Method)) + 
    geom_hline(yintercept = 1, linewidth = 0.5, color = 'black') +
    scale_color_d3(guide = 'none') +
    scale_fill_d3() +
    labs(x = expression(paste('Confounding Coefficient ', rho)),
         y = expression(paste('Average Treatment Effect ', theta))) +
    facet_grid(d_z ~ SNR, scales = 'free', 
               labeller = label_bquote(italic(d[Z])==.(d_z))) +
    theme_bw() + 
    theme(axis.title = element_text(size = 16),
          strip.text.x = element_text(size = 16),
          strip.text.y = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.position = 'bottom')
  ggsave(paste0('./plots/benchmark_', z, '_SNRx=', s, '.pdf'), width = 10)
}
foreach(zz = c('Diagonal', 'Toeplitz')) %:%
  foreach(ss = c(1/2, 1, 2)) %do% plot_fn(zz, ss)

################################################################################

# Proportion of valid instruments
df <- foreach(zz = c(0, 0.5), .combine = rbind) %:%
  foreach(s_x = c(1/2, 1, 2), .combine = rbind) %:%
  foreach(s_y = c(1/2, 1, 2), .combine = rbind) %:%
  foreach(pp = seq(0, 0.95, 0.05), .combine = rbind) %dopar%
  bnchmrk(d_z = 20, zz, rho = 0.5, s_x, s_y, pp, n = 1000, n_sim = 50)
saveRDS(df, './bnchmrks/pr_valid_sim.rds')

# Prepare for plotting
df[, mu := mean(theta, na.rm = TRUE), by = .(method, z_rho, snr_x, snr_y, pr_valid)]
df[, se := sd(theta, na.rm = TRUE), by = .(method, z_rho, snr_x, snr_y, pr_valid)]
tmp <- unique(df[, .(method, mu, se, z_rho, snr_x, snr_y, pr_valid)])
tmp[, z_rho := fifelse(z_rho == 0, 'Diagonal', 'Toeplitz')]
tmp[snr_y == 1/2, SNR := 'SNR = 1/2']
tmp[snr_y == 1, SNR := 'SNR = 1']
tmp[snr_y == 2, SNR := 'SNR = 2']
tmp[, SNR := factor(SNR, levels = c('SNR = 1/2', 'SNR = 1', 'SNR = 2'))]
tmp[method == 'TSLS', method := '2SLS']
setnames(tmp, 'method', 'Method')
tmp[, lo := .SD[Method == 'Lower', mu], by = .(z_rho, snr_x, snr_y, pr_valid)]
tmp[, hi := .SD[Method == 'Upper', mu], by = .(z_rho, snr_x, snr_y, pr_valid)]
tmp2 <- tmp[!Method %in% c('Lower', 'Upper')]
tmp2[, Method := factor(Method, levels = c('Backdoor', '2SLS', 'sisVIVE', 'MBE'))]

plot_fn <- function(s) {
  p1 <- ggplot(tmp2[snr_x == s], aes(pr_valid, mu, fill = Method)) + 
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = 'LeakyIV'), alpha = 0.25) +
    geom_ribbon(aes(ymin = mu - se, ymax = mu + se), alpha = 0.5) + 
    geom_line(aes(color = Method)) + 
    geom_hline(yintercept = 1, linewidth = 0.5, color = 'black') +
    scale_color_d3(guide = 'none') +
    scale_fill_d3() +
    labs(x = 'Proportion of Valid Instruments',
         y = expression(paste('Average Treatment Effect ', theta))) +
    facet_grid(z_rho ~ SNR, scales = 'free') +
    theme_bw() + 
    theme(axis.title = element_text(size = 16),
          strip.text.x = element_text(size = 16),
          strip.text.y = element_text(size = 16),
          legend.title = element_text(size = 16),
          legend.text = element_text(size = 16),
          legend.position = 'bottom')
  ggsave(paste0('./plots/prop_valid_SNRx=', s, '.pdf'), width = 10)
}
sapply(c(1/2, 1, 2), plot_fn)

################################################################################

# Benchmark against Bayesian baseline

# Initialize
out <- data.table(
  b = numeric(), method = character(), theta = numeric(),
  d_z = numeric(), z_rho = numeric(), rho = numeric(), 
  snr_x = numeric(), snr_y = numeric(), pr_valid = numeric()
)
saveRDS(out, './bnchmrks/bayesian.rds')

bayes_bnchmrk <- function(d_z, z_rho, rho, snr_x, snr_y, pr_valid, n, n_sim) {
  
  # Update
  cat(paste('Running:', d_z, z_rho, rho, snr_x, snr_y, '\n'))
  
  # Simulate population covariance matrix
  sim <- sim_cov(d_z, z_rho, rho, theta = 1, snr_x, snr_y, pr_valid)
  Sigma <- sim$Sigma
  l2 <- sqrt(sum(sim$params$gamma^2))
  d <- d_z + 2
  tau <- 1.1 * l2
  
  # Treat this as ground truth
  dat <- as.data.table(rmvn(n * n_sim, mu = rep(0, d), Sigma))
  colnames(dat) <- c('x', 'y', paste0('z', seq_len(d_z)))
  
  # Inner loop
  inner_loop <- function(b) {
    
    # Draw data
    tmp <- dat[((b - 1) * n + 1):(b * n), ]
    
    # Bayesian baseline
    bb_res <- baseline_gaussian_sample(
      m = 2000, dat = as.matrix(tmp), pos_z = 3:d, pos_x = 1,
      pos_y = 2, prop_var = 0.01, lvar_z_mu = -1, lvar_z_var = 3,
      beta_var = 5, pre_gamma_var = 1, theta_var = 5, 
      lvar_error_mu = -1, lvar_error_var = 1, tau = tau, alpha = 0.05
    )
    ate_bb <- mean(bb_res$theta)
    ate_lo_bb <- bb_res$lb_interval[2]
    ate_hi_bb <- bb_res$ub_interval[1]
    
    # LeakyIV
    suppressWarnings(
      ate_bnds <- leakyIV(tmp, tau = tau, method = 'mle')
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
  out[, d_z := d_z][, z_rho := z_rho][, rho := rho][, 
    snr_x := snr_x][, snr_y := snr_y][, pr_valid := pr_valid]
  old <- readRDS('./bnchmrks/bayesian.rds')
  new <- rbind(old, out)
  saveRDS(new, './bnchmrks/bayesian.rds')
}

# Execute in parallel over a grid
foreach(dd = c(5, 10), .combine = rbind) %:%
  foreach(zz = c(0, 1/2), .combine = rbind) %:%
  foreach(rr = seq(-0.9, 0.9, 0.1), .combine = rbind) %:%
  foreach(s_y = c(1/2, 1, 2), .combine = rbind) %dopar%
  bayes_bnchmrk(dd, zz, rr, snr_x = 2, s_y, pr_valid = 1/5, 
                n = 1000, n_sim = 50)




df <- foreach(dd = c(5, 10), .combine = rbind) %:%
  foreach(zz = c(0, 0.5), .combine = rbind) %:%
  foreach(rr = seq(-0.9, 0.9, 0.1), .combine = rbind) %:%
  foreach(s_x = c(1/2, 1, 2), .combine = rbind) %:%
  foreach(s_y = c(1/2, 1, 2), .combine = rbind) %dopar%
  bnchmrk(dd, zz, rr, s_x, s_y, pr_valid = 1/5, n = 1000, n_sim = 50)

################################################################################

# Spot check
df[method %in% c('leakyIV_lo', 'leakyIV_hi'), sum(is.na(theta)), 
   by = .(d_z, z_rho, rho, snr_x, snr_y)]

# Plot it
df[, mu := mean(theta, na.rm = TRUE), by = .(method, d_z, z_rho, rho, snr_x, snr_y)]
df[, se := sd(theta, na.rm = TRUE), by = .(method, d_z, z_rho, rho, snr_x, snr_y)]
tmp <- unique(df[, .(method, mu, se, d_z, z_rho, rho, snr_x, snr_y)])
tmp[, z_rho := fifelse(z_rho == 0, 'Diagonal', 'Toeplitz')]
tmp[snr_y == 1/2, SNR := 'SNR = 1/2']
tmp[snr_y == 1, SNR := 'SNR = 1']
tmp[snr_y == 2, SNR := 'SNR = 2']
tmp[, SNR := factor(SNR, levels = c('SNR = 1/2', 'SNR = 1', 'SNR = 2'))]
setnames(tmp, 'method', 'Method')
tmp[, leaky_lo := .SD[Method == 'leakyIV_lo', mu], by = .(d_z, z_rho, rho, snr_x, snr_y)]
tmp[, leaky_hi := .SD[Method == 'leakyIV_hi', mu], by = .(d_z, z_rho, rho, snr_x, snr_y)]
tmp[, bayes_lo := .SD[Method == 'Bayes_lo', mu], by = .(d_z, z_rho, rho, snr_x, snr_y)]
tmp[, bayes_hi := .SD[Method == 'Bayes_hi', mu], by = .(d_z, z_rho, rho, snr_x, snr_y)]
tmp2 <- tmp[!(grepl('lo', Method) | grepl('hi', Method))]


p2 <- ggplot(tmp2[z_rho == 'Diagonal' & snr_x == 2], aes(rho, mu)) + 
  geom_ribbon(aes(ymin = leaky_lo, ymax = leaky_hi, fill = 'LeakyIV'), alpha = 0.25) +
  geom_ribbon(aes(ymin = bayes_lo, ymax = bayes_hi, fill = 'Bayes'), alpha = 0.25) + 
  #geom_ribbon(aes(ymin = mu - se, ymax = mu + se), alpha = 0.5) +
  geom_line(aes(color = 'Bayes')) + 
  geom_hline(yintercept = 1, linewidth = 0.5, color = 'black') +
  scale_color_d3(guide = 'none') +
  scale_fill_d3() +
  labs(x = expression(paste('Confounding Coefficient ', rho)),
       y = expression(paste('Average Treatment Effect ', theta))) +
  facet_grid(d_z ~ SNR, scales = 'free', 
             labeller = label_bquote(italic(d[Z])==.(d_z))) +
  theme_bw() + 
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = 'bottom')



p1 <- ggplot(tmp2[SNR == 'SNR = 2'], aes(rho, mu, fill = Method)) + 
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = 'LeakyIV'), alpha = 0.25) +
  geom_ribbon(aes(ymin = mu - se, ymax = mu + se), alpha = 0.5) + 
  geom_line(aes(color = Method)) + 
  geom_hline(yintercept = 1, linewidth = 0.5, color = 'black') +
  scale_color_d3(guide = 'none') +
  scale_fill_d3() +
  labs(x = expression(paste('Confounding Coefficient ', rho)),
       y = expression(paste('Average Treatment Effect ', theta))) +
  facet_grid(z_rho ~ d_z, scales = 'free', 
             labeller = label_bquote(cols = italic(d[Z])==.(d_z))) +
  theme_bw() + 
  theme(axis.title = element_text(size = 12),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.position = 'bottom')







