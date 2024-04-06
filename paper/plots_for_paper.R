# Set working directory
setwd('~/Documents/Kings/leakyIV/paper')

# Load libraries
library(leakyIV)
library(data.table)
library(glasso)
library(mvnfast)
library(broom)
library(sisVIVE)
library(ggplot2)
library(ggsci)
library(scales)
library(cowplot)
library(colorspace)
library(ggpattern)
library(rags2ridges)
library(kde1d)
library(corpcor)
library(Matrix)
library(doMC)
registerDoMC(8)

# Load scripts
source('bayesian_baseline.R')
source('MBE.R')
source('simulator.R')

# Set seed
set.seed(123)

################################################################################

### FIG. 3 ###

# Compute theta and gamma norms under a range of values for rho
grid_fn <- function(n, d_z, z_rho, rho, theta, snr_x, snr_y, pr_valid, n_boot) {
  
  # Generate data
  sim <- sim_dat(n, d_z, z_rho, rho, theta, snr_x, snr_y, pr_valid)
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
    Theta_z <- solve(Sigma[3:d, 3:d])
    Sigma_zy <- matrix(Sigma[3:d, 2L], ncol = 1L)
    Sigma_yz <- t(Sigma_zy)
    Sigma_zx <- matrix(Sigma[3:d, 1L], ncol = 1L)
    Sigma_xz <- t(Sigma_zx)
    sigma_xy <- Sigma[1L, 2L]
    var_x <- Sigma[1L, 1L]
    var_y <- Sigma[2L, 2L]
    
    # Gamma is constrained to the hyperplane alpha - theta * beta
    alpha <- as.numeric(Theta_z %*% Sigma_zy)
    beta <- as.numeric(Theta_z %*% Sigma_zx)
    
    # Conditional (co)variances given Z
    k_xx <- var_x - as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
    k_yy <- var_y - as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
    k_xy <- sigma_xy - as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
    
    # Compute theta as a function of rho
    theta_fn <- function(rho) {
      (k_xy - sqrt(k_xx * k_yy - k_xy^2) * tan(asin(rho))) / k_xx
    }
    # Compute gamma norms as a function of rho
    norm_fn <- function(rho, p) {
      theta <- theta_fn(rho)
      gamma <- alpha - theta * beta
      norm <- (sum(abs(gamma)^p))^(1 / p)
      return(norm)
    }
    
    # Loop through a bunch of rho-values
    df <- data.table(rho = seq(-0.95, 0.95, length.out = 1000))
    df[, theta := sapply(rho, theta_fn)]
    df[, l1 := sapply(rho, norm_fn, 1)]
    df[, l2 := sapply(rho, norm_fn, 2)]
    
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

# Plot the rho-theta curve with bootstrap confidence intervals
df <- grid_fn(n = 1000, d_z = 4, z_rho = 0, rho = 1/2, theta = 1, 
              snr_x = 2, snr_y = 2, pr_valid = 0, n_boot = 1000)
df[, theta_mu := mean(theta), by = rho]
df[, theta_se := sd(theta), by = rho]
df[, l2_mu := mean(l2), by = rho]
df[, l2_se := sd(l2), by = rho]
df <- unique(df[, .(rho, theta_mu, theta_se, l2_mu, l2_se)])
p1 <- ggplot(df, aes(rho, theta_mu)) +
  geom_ribbon(aes(ymin = theta_mu - theta_se * qnorm(0.975),
                  ymax = theta_mu + theta_se * qnorm(0.975)),
              color = 'grey90', alpha = 0.25) +
  geom_line(linewidth = 0.75) +
  labs(x = expression(paste('Confounding Coefficient ', rho)),
       y = expression(paste('Average Treatment Effect ', theta))) +
  theme_bw() +
  theme(axis.title = element_text(size = 24))

# Plot the theta-L2 curve with bootstrap confidence intervals
p2 <- ggplot(df, aes(theta_mu, l2_mu)) +
  geom_ribbon(aes(ymin = l2_mu - l2_se * qnorm(0.975),
                  ymax = l2_mu + l2_se * qnorm(0.975)),
              color = 'grey90', alpha = 0.25) +
  geom_line(linewidth = 0.75) +
  labs(x = expression(paste('Average Treatment Effect ', theta)),
       y = expression(paste(L[2], ' Norm, Leakage Weights ', gamma))) +
  theme_bw() +
  theme(axis.title = element_text(size = 24))

plot_grid(p1, p2, labels = c('A', 'B'), label_size = 24)
ggsave2('./plots/lemmas.pdf', height = 7, width = 14)

################################################################################

### FIG. 4 ###

# Illustrate the three-partition of threshold space
set.seed(999)
sim <- sim_cov(d_z = 4, z_rho = 0, rho = 1/2, theta = 1, 
               snr_x = 2, snr_y = 2, pr_valid = 0)
Sigma <- sim$Sigma
l2 <- sqrt(sum(sim$params$gamma^2))
d <- ncol(Sigma)
d_z <- d - 2L

# Find minimum possible tau
Theta_z <- solve(Sigma[3:d, 3:d])
Sigma_zy <- matrix(Sigma[3:d, 2L], ncol = 1L)
Sigma_yz <- t(Sigma_zy)
Sigma_zx <- matrix(Sigma[3:d, 1L], ncol = 1L)
Sigma_xz <- t(Sigma_zx)
sigma_xy <- Sigma[1L, 2L]
var_x <- Sigma[1L, 1L]
var_y <- Sigma[2L, 2L]
alpha <- as.numeric(Theta_z %*% Sigma_zy)
beta <- as.numeric(Theta_z %*% Sigma_zx)
f <- lm(alpha ~ 0 + beta)
tau_check <- sqrt(sum(residuals(f)^2))
min_fctr <- tau_check / l2

# Conditional (co)variances given Z
k_xx <- var_x - as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
k_yy <- var_y - as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
k_xy <- sigma_xy - as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)

# Compute theta as a function of rho
theta_fn <- function(rho) {
  (k_xy - sqrt(k_xx * k_yy - k_xy^2) * tan(asin(rho))) / k_xx
}

# Compute gamma norm a function of theta
norm_fn2 <- function(theta, p = 2) {
  gamma <- alpha - theta * beta
  if (p == Inf) {
    norm <- max(abs(gamma))
  } else {
    norm <- (sum(abs(gamma)^p))^(1 / p)
  }
  return(norm)
}
tau_star <- norm_fn2(2)

df <- data.table(rho = seq(-0.95, 0.95, length.out = 1000))
df[, theta := sapply(rho, theta_fn)]
df[, l2 := sapply(theta, norm_fn2)]
df[, l1 := sapply(theta, norm_fn2, 1)]

# Restrict the space
min_theta <- df[, min(theta)]
tmp <- df[theta <= -min_theta + 1]
tmp[, theta := theta + 2.5]

# Plot regions
p1 <- ggplot(tmp, aes(theta, l2)) + 
  geom_ribbon_pattern(aes(ymin = 4.5, ymax = tau_check), 
                      pattern = 'stripe', fill = 'grey', alpha = 0.5) +
  geom_ribbon(aes(ymin = tau_check, ymax = tau_star), fill = 'red', alpha = 0.5) +
  geom_line(linewidth = 1) +
  geom_vline(xintercept = 4.5, linewidth = 1, color = 'blue') + 
  geom_hline(yintercept = tau_check, linetype = 'dashed', linewidth = 1) + 
  geom_hline(yintercept = tau_star, linetype = 'dashed', linewidth = 1) + 
  ylim(4.5, 9) +
  labs(x = expression(paste('Average Treatment Effect ', theta)),
       y = expression(paste(L[2], ' Norm, Leakage Weights ', gamma))) +
  theme_bw() +
  theme(axis.title = element_text(size = 24))
ggsave('./plots/3partition.pdf')

################################################################################

### FIG. 5 ###

# Benchmark against backdoor adjustment, 2SLS, sisVIVE, and MBE
bnchmrk <- function(d_z, z_rho, rho, snr_x, snr_y, pr_valid, n, n_sim) {
  
  # Update
  cat(paste('Running:', d_z, z_rho, rho, snr_x, snr_y, '\n'))
  
  # Simulate population covariance matrix, treat as ground truth
  sim <- sim_cov(d_z, z_rho, rho, theta = 1, snr_x, snr_y, pr_valid)
  Sigma <- sim$Sigma
  l2 <- sqrt(sum(sim$params$gamma^2))
  d <- d_z + 2
  tau <- 1.1 * l2
  
  # Draw samples from this fixed data generating process
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
    
    # LeakyIV
    suppressWarnings(
      ate_leaky <- leakyIV(tmp, tau = tau, method = 'mle', normalize = FALSE)
    )
    ate_leaky_lo <- ate_leaky$ATE_lo
    ate_leaky_hi <- ate_leaky$ATE_hi
    
    # Export
    out <- data.table(b,
      method = c('Backdoor', 'TSLS', 'sisVIVE', 'MBE', 'leaky_lo', 'leaky_hi'),
      theta = c(ate_bda, ate_2sls, ate_sisvive, ate_mbe, ate_leaky_lo, ate_leaky_hi)
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

### FIG. 6 ###

# Benchmark against Bayesian baseline

bayes_bnchmrk <- function(rho, n) {
  
  # Simulate population covariance matrix, treat as ground truth
  sim <- sim_cov(d_z = 5, z_rho = 0, rho, theta = 1, 
                 snr_x = 2, snr_y = 2, pr_valid = 1/5)
  Sigma <- sim$Sigma
  l2 <- sqrt(sum(sim$params$gamma^2))
  d <- ncol(Sigma)
  d_z <- d - 2
  tau <- 1.1 * l2
  
  # Draw samples from this fixed data generating process
  dat <- as.data.table(rmvn(n, mu = rep(0, d), Sigma))
  colnames(dat) <- c('x', 'y', paste0('z', seq_len(d_z)))
  
  # Bayesian baseline
  bb_res <- baseline_gaussian_sample(
    m = 2000, dat = as.matrix(dat), pos_z = 3:d, pos_x = 1,
    pos_y = 2, prop_var = 0.01, lvar_z_mu = -1, lvar_z_var = 3,
    beta_var = 5, pre_gamma_var = 5, theta_var = 5, 
    lvar_error_mu = -1, lvar_error_var = 5, tau = sqrt(tau), alpha = 0.05
  ) # The sqrt is because this is parameterized for squared Euclidean distance
  
  # LeakyIV
  suppressWarnings(
    ate_bnds <- leakyIV(dat, tau = tau, method = 'mle', normalize = FALSE)
  )
  
  # Export
  out <- data.table(rho,
    Bayes = bb_res$theta, 
    LeakyIV_lo = ate_bnds$ATE_lo, 
    LeakyIV_hi = ate_bnds$ATE_hi
  )
  return(out)
  
}
set.seed(222)
df <- foreach(rr = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75), 
              .combine = rbind) %dopar% bayes_bnchmrk(rr, n = 1000)

# Plot it
setnames(df, 'LeakyIV_lo', 'lo')
setnames(df, 'LeakyIV_hi', 'hi')
tmp <- melt(df, id.vars = c('rho', 'lo', 'hi'), variable.name = 'Method')
pl <- ggplot(tmp, aes(value)) + 
  geom_histogram(bins = 40, color = 'black', alpha = 0.5) + 
  geom_vline(xintercept = 1, linewidth = 0.5, color = 'blue') +
  geom_vline(aes(xintercept = lo), color = 'red', linetype = 'dashed') + 
  geom_vline(aes(xintercept = hi), color = 'red', linetype = 'dashed') + 
  scale_fill_npg() +
  labs(x = expression(paste('MCMC Posterior for the ATE ', theta)), 
       y = 'Count') +
  facet_wrap(~ rho, scales = 'free', ncol = 3, 
             labeller = label_bquote(rho == .(rho))) +
  theme_bw() + 
  theme(axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.position = 'bottom')
ggsave('./plots/bayesian.pdf', width = 10)

################################################################################

### FIG. 7 ###

# Testing exclusion 

# Test statistic function
psi_fn <- function(S) {
  Lambda <- S[3:nrow(S), 1:2]
  psi <- det(crossprod(Lambda))
  return(psi)
}

loop <- function(b, tau_check, n, rho, n_sim = 2000) {
  
  # Simulate population covariance matrix, treat as ground truth
  sim <- sim_pwr(d_z = 5, z_rho = 0, rho, theta = 1, snr_x = 2, snr_y = 2, 
                 tau_check)
  Sigma <- sim$Sigma
  d <- ncol(Sigma)
  d_z <- d - 2
  
  # Simulate data, estimate psi
  dat <- as.data.table(rmvn(n, mu = rep(0, d), Sigma))
  colnames(dat) <- c('x', 'y', paste0('z', seq_len(d_z)))
  Sigma <- Sigma0 <- cov(dat)
  psi <- psi_fn(Sigma)
  
  # Impose H0
  Sigma_zx <- Sigma[3:d, 1L]
  Sigma_zy <- Sigma[3:d, 2L]
  theta_2sls <- as.numeric((Sigma_zx %*% Sigma_zy) / (Sigma_zx %*% Sigma_zx))
  Sigma0[3:d, 2L] <- Sigma0[2L, 3:d] <- theta_2sls * Sigma_zx 
  
  # Positive definite check
  if (!is.positive.definite(Sigma0)) {
    Sigma0 <- as.matrix(nearPD(Sigma0)$mat)
  }
  
  # Draw MC samples
  null <- as.data.table(rmvn(n * n_sim, mu = rep(0, d), Sigma0))
  colnames(null) <- c('x', 'y', paste0('z', seq_len(d_z)))
  
  # Simulate null distro
  psi0 <- sapply(seq_len(n_sim), function(b) {
    tmp <- null[((b - 1) * n + 1):(b * n), ]
    psi_fn(tmp)
  })
  
  # Compute p-value
  p_value <- (sum(psi0 >= psi) + 1L) / (n_sim + 1L)
  
  # Export
  out <- data.table(
    tau_check, n, rho, p_value, b
  )
  return(out)
  
}
res <- foreach(bb = 1:500, .combine = rbind) %:%
  foreach(tt = seq(0, 1, length.out = 21), .combine = rbind) %:%
  foreach(nn = c(500, 1000, 2000), .combine = rbind) %:%
  foreach(rr = c(-3/4, -1/2, -1/4, 1/4, 1/2, 3/4), .combine = rbind) %dopar%
  loop(bb, tt, nn, rr)

# Prepare for plotting
df <- na.omit(res)
df[, hit := fifelse(p_value <= 0.1, 1L, 0L)]
df[, prop_hit := mean(hit), by = .(tau_check, n, rho)]
df[, se_hit := sd(hit) / sqrt(.N), by = .(tau_check, n, rho)]
df[, n := as.factor(n)]
tmp <- unique(df[, .(tau_check, prop_hit, se_hit, n, rho)])

# Plot
pl <- ggplot(tmp, aes(tau_check, prop_hit, fill = n)) + 
  geom_ribbon(aes(ymin = prop_hit - se_hit, ymax = prop_hit + se_hit),
              alpha = 0.5) +
  geom_point(aes(color = n)) + 
  geom_line(aes(color = n), linewidth = 1) + 
  geom_hline(yintercept = 0.1, color = 'red', linetype = 'dashed') + 
  scale_color_discrete_qualitative(palette = 'Dark 3') +
  scale_fill_discrete_qualitative(palette = 'Dark 3') +
  labs(x = 'Minimum Leakage   ', y = 'Rejection Rate') +
  facet_wrap(~ rho, labeller = label_bquote(rho == .(rho))) +
  theme_bw() + 
  theme(legend.position = 'bottom',
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 16),
        axis.title = element_text(size = 16),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16))
ggsave('./plots/power.pdf', width = 10)

################################################################################

### FIG. 8 ###

# Check coverage of estimated bounds across a range of settings
cover <- function(d_z, rho, snr_x, snr_y, alpha, n, b) {
  
  # Simulate population covariance matrix, treat as ground truth
  sim <- sim_cov(d_z, z_rho = 0, rho, theta = 1, snr_x, snr_y, pr_valid = 0)
  d <- d_z + 2
  Sigma <- sim$Sigma
  l2 <- sqrt(sum(sim$params$gamma^2))
  tau <- 1.1 * l2
  S_oracle <- leakyIV(Sigma, tau, normalize = FALSE)
  o_lo <- S_oracle$ATE_lo
  o_hi <- S_oracle$ATE_hi
  
  # Draw samples from this fixed data generating process
  dat <- as.data.table(rmvn(n, mu = rep(0, d), Sigma))
  colnames(dat) <- c('x', 'y', paste0('z', seq_len(d_z)))
  
  # How bad is the estimate?
  m <- colMeans(dat)
  s <- cov(dat)
  kl <- KLdiv(m, rep(0, d), s, Sigma)
  
  # Estimate bounds and quantiles
  pr <- c(alpha/2, 1 - alpha/2)
  suppressWarnings(
    res <- leakyIV(dat, tau, n_boot = 2000, parallel = FALSE, normalize = FALSE)
  )
  qf_lo <- quantile(res$ATE_lo, probs = pr, na.rm = TRUE)
  qf_hi <- quantile(res$ATE_hi, probs = pr, na.rm = TRUE)
  kde_lo <- kde1d(res[!is.na(ATE_lo), ATE_lo])
  kde_hi <- kde1d(res[!is.na(ATE_hi), ATE_hi])
  qd_lo <- qkde1d(pr, kde_lo)
  qd_hi <- qkde1d(pr, kde_hi)
  mu_lo <- res[, mean(ATE_lo, na.rm = TRUE)]
  sd_lo <- res[, sd(ATE_lo, na.rm = TRUE)]
  mu_hi <- res[, mean(ATE_hi, na.rm = TRUE)]
  sd_hi <- res[, sd(ATE_hi, na.rm = TRUE)]
  qn_lo <- qnorm(pr, mu_lo, sd_lo)
  qn_hi <- qnorm(pr, mu_hi, sd_hi)
  
  # Export
  out <- data.table(
    method = rep(c('Standard', 'KDE', 'Gaussian'), each = 2),
    bound = rep(c('Lower', 'Upper'), times = 3),
    hit = c(
      fifelse(o_lo >= qf_lo[1] & o_lo <= qf_lo[2], 1L, 0L),
      fifelse(o_hi >= qf_hi[1] & o_hi <= qf_hi[2], 1L, 0L),
      fifelse(o_lo >= qd_lo[1] & o_lo <= qd_lo[2], 1L, 0L),
      fifelse(o_hi >= qd_hi[1] & o_hi <= qd_hi[2], 1L, 0L),
      fifelse(o_lo >= qn_lo[1] & o_lo <= qn_lo[2], 1L, 0L),
      fifelse(o_hi >= qn_hi[1] & o_hi <= qn_hi[2], 1L, 0L)
    ),
    kldiv = kl, bnd_cor = res[, cor(ATE_lo, ATE_hi, use = 'complete.obs')], tau
  )
  out[, rho := rho]
  return(out)
  
}
set.seed(123)
df <- foreach(rr = c(-0.75, -0.5, -0.25, 0.25, 0.5, 0.75), .combine = rbind) %:%
  foreach(bb = 1:500, .combine = rbind) %dopar%
  cover(d_z = 4, rr, snr_x = 2, snr_y = 2, alpha = 0.1, n = 1000, bb)

# Prepare for plotting
tmp <- df[, mean(hit, na.rm = TRUE), by = .(method, bound, rho)]
setnames(tmp, 'V1', 'cvg')
tmp[, se := df[, sd(hit, na.rm = TRUE) / sqrt(.N), by = .(method, bound, rho)]$V1]
tmp[, method := factor(method, levels = c('Standard', 'KDE', 'Gaussian'))]

# Plot it
pl <- ggplot(tmp, aes(method, cvg, fill = method)) + 
  geom_bar(stat = 'identity') + 
  geom_errorbar(aes(ymin = cvg - se, ymax = cvg + se), width = 0.5) + 
  geom_hline(yintercept = 0.9, linetype = 'dashed', color = 'red') + 
  scale_fill_aaas() +
  scale_y_continuous(limits = c(0.85, 0.95), oob = rescale_none) +
  labs(x = 'Bootstrap Estimator', y = 'Coverage') + 
  facet_grid(bound ~ rho, labeller = label_bquote(cols = rho == .(rho))) +
  theme_bw() + 
  theme(legend.position = 'none',
        axis.title = element_text(size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.text.x = element_text(size = 16),
        strip.text.y = element_text(size = 16))
ggsave('./plots/coverage.pdf', width = 12)  


