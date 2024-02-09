# Set working directory
setwd('C:/Users/k23067841/OneDrive - King\'s College London/Documents/leakyIV2/benchmark_simulations')

# Load libraries, register cores
library(reshape2)
library(data.table)
library(broom)
library(sisVIVE)
library(tidyverse)
library(ggsci)
library(doParallel)
library(LambertW)

cl <- makeCluster(8)
registerDoParallel(cl)

# Load MBE script
source('MBE.R')

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Import sim_idx
sim_idx <- fread('./prop_valid_n50_1000_per_subidx_1/sim_idx.csv')

# Benchmark against sisVIVE and 2SLS:
bnchmrk <- function(sim_idx, sub_idx) {
  
  setwd('C:/Users/k23067841/OneDrive - King\'s College London/Documents/leakyIV2/benchmark_simulations')
  
  source('MBE.R')
  source("leakyIV.R")
  
  # Read the data
  dat <- data.table::fread(paste0('./prop_valid_n50_1000_per_subidx_1/', sim_idx, '.csv'))
  
  d_z <- (ncol(dat) - 3)/2
  
  print(d_z)
  
  # Subset the data
  tmp_with_gamma_theta <- dat[sub_idx:(sub_idx + 999)]
  tmp <- dplyr::select(tmp_with_gamma_theta, 1:(d_z+2))
  theta <- dplyr::select(tmp_with_gamma_theta, (d_z+3))
  gamma <- dplyr::select(tmp_with_gamma_theta, -(1:(d_z+3)))
  
  #print(tmp)
  #print(gamma)
  #print(theta)
  
  # Run 2SLS
  f0 <- stats::lm(x ~ ., data = dplyr::select(tmp, -y))
  tmp$x_hat <- stats::fitted(f0)
  f1 <- stats::lm(y ~ x_hat, data = tmp)
  
  #print(f1)
  
  #print(broom::tidy(f1))
  
  ace_2sls <- broom::tidy(f1)$estimate[2]
  std_err_2sls <- broom::tidy(f1)$std.error[2]
  
  
  
  # Run sisVIVE
  sisvive <- sisVIVE::cv.sisVIVE(tmp$y, tmp$x, as.matrix(dplyr::select(tmp, starts_with('z'))))
  ace_sisVIVE <- sisvive$beta
  
  #print(generics::tidy(f0))
  
  # Run MBE
  beta_hat <- generics::tidy(f0)$estimate[2:d_z+1]
  #print(beta_hat)
  
  se_beta <- generics::tidy(f0)$std.error[2:d_z+1]
  #print(se_beta)
  
  f2 <- stats::lm(y ~ ., data = dplyr::select(tmp, -c(x_hat,x)))
  
  #print(generics::tidy(f2))
  
  gamma_hat <- generics::tidy(f2)$estimate[2:d_z+1]
  #print(gamma_hat)
  
  se_gamma <- generics::tidy(f2)$std.error[2:d_z+1]
  #print(se_gamma)
  
  mbe <- MBE(beta_hat, gamma_hat, se_beta, se_gamma, phi = 1, n_boot = 1)
  
  #print(mbe)
  
  ace_mbe <- mbe$Estimate[2]
  
  print(mbe$SE)
  
  # Run leakyIV
  
  x_dat <- data.matrix(dplyr::select(tmp, x), rownames.force = NA)
  y_dat <- data.matrix(dplyr::select(tmp, y), rownames.force = NA)
  z_dat <- dplyr::select(tmp, -c(x, y, x_hat))
  
  #print(is.numeric(x_dat))
  
  #print(z_dat)
  
  known_theta <- as.numeric(theta[sim_idx])
  print(known_theta)
  
  
  norm_l2 <- LambertW::lp_norm(as.numeric(gamma[1]), 2)
  
  leaky_ace = leakyIV(x_dat, y_dat, z_dat, norm_l2)
  
  leaky_ace_lo = as.numeric(dplyr::select(leaky_ace, 1))
  leaky_ace_hi = as.numeric(dplyr::select(leaky_ace, 2))
  
  print(sim_idx)
  print(sub_idx)
  print(norm_l2)
  print(known_theta)
  print(std_err_2sls)
  print(ace_sisVIVE)
  print(ace_mbe)
  print(leaky_ace_lo)
  print(leaky_ace_hi)
  
  
  # Export
  out <- data.table(
    sim_idx = sim_idx, sub_idx = sub_idx, norm_gamma = norm_l2, tau_div_norm = 1, norm = 'L2', regularization = 'none (MLE)', 
    theta_true = known_theta,
    'TSLS' = ace_2sls, 'sd_err_TSLS' = std_err_2sls, 'sisVIVE' = ace_sisVIVE, 'MBE' = ace_mbe,
    ATE_lo = leaky_ace_lo, ATE_hi = leaky_ace_hi
  )
  
  print(paste0('sim_idx: ', toString(sim_idx), '  sub_idx: ', toString(sub_idx)))
  
  # tau = d ||gamma||_2
  
  #  for(d in c(1.01, 1.05, 1.1, 1.5, 2)){
  
  #    tau_l2 <- d*norm_l2
  
  #    leaky_ace = leakyIV(x_dat, y_dat, z_dat, tau_l2)
  
  #    leaky_ace_lo = as.numeric(dplyr::select(leaky_ace, 1))
  #    leaky_ace_hi = as.numeric(dplyr::select(leaky_ace, 2))
  
  #    out <- rbind(out, list(sim_idx, sub_idx, norm_l2,  d, 'L2', 'none (MLE)',
  #                           known_theta, ace_2sls, std_err_2sls,  ace_sisVIVE, ace_mbe, std_err_mbe,
  #                           leaky_ace_lo, leaky_ace_hi))
  
  #  }
  
  
  return(out)
}

bnchmrk(1, 4001)

# Execute in parallel
res <- foreach(aa = seq_len(84), .combine = rbind, .packages = c("data.table", "broom", "sisVIVE", 
                                                                "tidyverse", "ggsci", "LambertW", "reshape2", "tidyr")) %:%
  foreach(bb = seq(1, 49001, 1000), .combine = rbind) %dopar%
  #prev_time <- Sys.time()
  #curr_time <- Sys.time()
  bnchmrk(aa, bb)

print(res)

fwrite(res, './prop_valid_n50_1000_per_subidx_1/results.csv')

# Plot results
df <- res %>%
  select(res, -sub_idx) %>%
  rename(idx = sim_idx) %>%
  pivot_longer(TSLS:sisVIVE, names_to = 'Method', values_to = 'alpha') %>%
  inner_join(sim_idx, by = 'idx')


df %>%
  filter(z_rho == 0) %>%
  ggplot(aes(Method, alpha, color = Method)) + 
  geom_jitter(size = 0.75) + 
  geom_boxplot() + 
  scale_color_d3() + 
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') + 
  theme_bw() + 
  facet_grid(rho ~ prop_valid)

df %>%
  filter(z_rho != 0) %>%
  ggplot(aes(Method, alpha, fill = Method)) + 
  geom_boxplot() + 
  scale_fill_d3() + 
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') + 
  theme_bw() + 
  facet_grid(rho ~ prop_valid)


df <- fread('./grid2largeB/sim_idx.csv')[1:24]
colnames(df) <- c('sim_idx', 'lo', 'hi')
df[is.na(lo), lo := hi]
df[is.na(hi), hi := lo]
df <- df %>% pivot_longer(-sim_idx, names_to = 'bound', values_to = 'soft_iv')
df$sim_idx <- as.numeric(df$sim_idx)
res2 <- res %>% 
  group_by(sim_idx) %>% 
  summarize(TSLS_mu = mean(TSLS), sisVIVE_mu = mean(sisVIVE), MBE_mu = mean(MBE),
            TSLS_sd = sd(TSLS), sisVIVE_sd = sd(sisVIVE), MBE_sd = sd(MBE))
colnames(res2) <- gsub('_mu', '', colnames(res2))
res2 <- res2 %>%
  pivot_longer(TSLS:MBE, names_to = 'Method', values_to = 'alpha')
colnames(res2) <- gsub('_sd', '', colnames(res2))
res3 <- res2 %>%
  select(sim_idx:MBE) %>%
  pivot_longer(-sim_idx, names_to = 'Method', values_to = 'se') %>%
  unique(.)
res2 <- res2 %>% 
  select(sim_idx, Method, alpha) %>%
  unique(.) %>%
  inner_join(res3, by = c('sim_idx', 'Method'))


colnames(sim_idx)[1] <- 'sim_idx'



df <- df %>% 
  select(-bound) %>%
  inner_join(res2, by = 'sim_idx') %>%
  inner_join(sim_idx, by = 'sim_idx') %>%
  pivot_longer(soft_iv:MBE, names_to = 'Method', values_to = 'alpha')

df <- df %>% 
  select(-bound) %>%
  inner_join(select(res, -sub_idx), by = 'sim_idx') %>%
  inner_join(sim_idx, by = 'sim_idx') %>%
  pivot_longer(soft_iv:MBE, names_to = 'Method', values_to = 'alpha') %>%
  mutate(Sigma_z = ifelse(z_rho == 0, 'identity', 'toeplitz'),
         prop_valid = paste('pr_valid =', prop_valid))

df <- as.data.table(df)
df[Method == 'TSLS', Method := '2SLS']
df[, Method := factor(Method, levels = c('soft_iv', '2SLS', 'sisVIVE', 'MBE'))]

df %>%
  ggplot(aes(as.factor(rho), alpha, fill = Method)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') + 
  scale_fill_d3() + 
  labs(x = 'Unobserved Confounding', y = 'Average Treatment Effect') +
  theme_bw() + 
  facet_grid(Sigma_z ~ prop_valid)


df %>%
  ggplot(aes(as.factor(rho), alpha, fill = Method)) +
  geom_bar(stat = 'identity', color = 'black', position = position_dodge()) +
  geom_errorbar(aes(ymin = alpha - se, ymax = alpha + se), width = .2,
                position = position_dodge(.9))
geom_hline(yintercept = 1, color = 'red', linetype = 'dashed') + 
  scale_fill_d3() + 
  labs(x = 'Unobserved Confounding', y = 'Average Treatment Effect') +
  theme_bw() + 
  facet_grid(Sigma_z ~ prop_valid)