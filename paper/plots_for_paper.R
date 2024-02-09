# Set working directory
setwd('~/Documents/Kings/leakyIV/paper')

# Load libraries
library(leakyIV)
library(data.table)
library(glasso)
library(mvnfast)
library(ggplot2)
library(ggsci)
library(cowplot)
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


################################################################################

# This function computes the predictions of a Sigma oracle for a given 
# dimensionality and factor on tau*, as well as n_sim leakyIV estimates at 
# sample size n.

outer_loop <- function(d_z, tau_fctr, n = 2000, n_sim = 200) {
  
  # Generate data, extract "population" covariance matrix
  sim <- sim_dat(n = 1e5, d_z, z_cnt = TRUE, z_rho = 0, rho = 1/4, 
                 theta = 1, r2_x = 3/4, r2_y = 3/4, pr_valid = 0)
  d <- d_z + 2
  l2 <- sqrt(sum(sim$params$gamma^2))
  
  # Treat this as ground truth
  Sigma <- cov(sim$dat)
  S_oracle <- leakyIV(sim$dat$x, sim$dat$y, sim$dat[, -c('x', 'y')],
                      tau = l2 * tau_fctr, method = 'mle')[, b := 0]
  dat <- rmvn(n * n_sim, mu = rep(0, d), Sigma)
  inner_loop <- function(b) {
    tmp <- dat[((b - 1) * n + 1):(b * n), ]
    z <- tmp[, seq_len(d_z)]
    x <- tmp[, d_z + 1]
    y <- tmp[, d_z + 2]
    if (d_z <= 5) {
      out <- leakyIV(x, y, z, tau = l2 * tau_fctr, method = 'mle')
    } else {
      out <- leakyIV(x, y, z, tau = l2 * tau_fctr, method = 'glasso', rho = 0.0015)
    }
    return(out[, b := b])
  }
  res <- foreach(bb = seq_len(n_sim), .combine = rbind) %do% inner_loop(bb)
  
  # Export
  out <- rbind(S_oracle, res)
  out <- melt(out, id.vars = 'b', variable.name = 'Bound', value.name = 'theta')
  out[, Bound := fifelse(grepl('lo', Bound), 'Lower', 'Upper')]
  out[, 'd_z' := d_z][, 'tau_fctr' := tau_fctr]
  return(out)
  
}

################################################################################

# Hold d_z fixed at 5 and run across different tau thresholds
df <- foreach(tt = seq(0.1, 2, by = 0.1), .combine = rbind) %dopar% outer_loop(5, tt)

# Compute means, standard errors
df[, sum(is.na(theta)), by = tau_fctr]
df <- na.omit(df)
mu <- df[b > 0, mean(theta), by = .(tau_fctr, Bound)]
setnames(mu, 'V1', 'mu')
se <- df[b > 0, sd(theta), by = .(tau_fctr, Bound)]
setnames(se, 'V1', 'se')
tmp <- merge(mu, se, by = c('tau_fctr', 'Bound'))
tmp[, Bound := factor(Bound, levels = c('Upper', 'Lower'))]
setnames(df, 'theta', 'mu')
tmp[, obs := 1]
sigma_oracle <- df[b==0, .(tau_fctr, Bound, mu)][, se := NA_real_][, obs := 0]
tmp <- rbind(tmp, sigma_oracle)


# Plot number of candidate instruments against ATE bounds
p1 <- ggplot() + 
  geom_vline(xintercept = 1, linetype = 'dashed', linewidth = 1) +
  geom_hline(yintercept = 1, linewidth = 1, color = 'grey60') +
  geom_ribbon(tmp[obs == 1], 
              mapping = aes(tau_fctr, mu, ymin = mu - se, ymax = mu + se,
                            fill = Bound), alpha = 0.4) + 
  geom_line(tmp[obs == 1], mapping = aes(tau_fctr, mu, color = Bound)) + 
  geom_point(tmp[obs == 0], mapping = aes(tau_fctr, mu, fill = Bound), 
             color = 'black', shape = 21, size = 5) + 
  scale_fill_npg() + 
  scale_color_npg() +
  labs(x = expression(paste('Leakage threshold ', tau)),
       y = expression(paste('Average Treatment Effect ', theta))) +
  theme_bw() + 
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.position = 'bottom')

# Hold tau_fctr fixed at 1.5 and run across different values of d_z
df <- foreach(dd = 1:20, .combine = rbind) %dopar% outer_loop(dd, 1.5)

# Compute means, standard errors
df[, sum(is.na(theta)), by = d_z]
df <- na.omit(df)
mu <- df[b > 0, mean(theta), by = .(d_z, Bound)]
setnames(mu, 'V1', 'mu')
se <- df[b > 0, sd(theta), by = .(d_z, Bound)]
setnames(se, 'V1', 'se')
tmp <- merge(mu, se, by = c('d_z', 'Bound'))
tmp[, Bound := factor(Bound, levels = c('Upper', 'Lower'))]
setnames(df, 'theta', 'mu')
tmp[, obs := 1]
true <- df[b==0, .(d_z, Bound, mu)][, se := NA_real_][, obs := 0]
tmp <- rbind(tmp, true)

# Plot number of candidate instruments against ATE bounds
p2 <- ggplot() + 
  geom_hline(yintercept = 1, linewidth = 1, color = 'grey60') +
  geom_ribbon(tmp[obs == 1], 
              mapping = aes(d_z, mu, ymin = mu - se, ymax = mu + se,
                            fill = Bound), alpha = 0.4) + 
  geom_line(tmp[obs == 1], mapping = aes(d_z, mu, color = Bound)) + 
  scale_fill_npg() + 
  scale_color_npg() +
  geom_point(tmp[obs == 0], mapping = aes(d_z, mu, fill = Bound), 
             color = 'black', shape = 21, size = 5) + 
  geom_hline(yintercept = 1, linewidth = 1, color = 'grey60') +
  labs(x = expression(paste('Number of Candidate Instruments ', d[Z])),
       y = expression(paste('Average Treatment Effect ', theta))) +
  theme_bw() + 
  theme(axis.title = element_text(size = 24),
        legend.text = element_text(size = 24),
        legend.position = 'bottom')

################################################################################

# Save both as a grid
pl <- plot_grid(p1, p2, labels = c('A', 'B'), label_size = 24)
ggsave2('./plots/Dz_and_tau_vs_ate.pdf', height = 7, width = 14)




