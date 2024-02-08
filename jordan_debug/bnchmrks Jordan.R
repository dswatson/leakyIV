# Set working directory
setwd('C:/Users/k23067841/Downloads/LeakyIV simulations again')

# Load libraries, register cores
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
sim_idx <- fread('./grid2_largeB/sim_idx.csv')

# Benchmark against sisVIVE and 2SLS:
bnchmrk <- function(sim_idx, sub_idx) {
  
  source('MBE.R')
  source("leakyIV.R")
  
  # Read the data
  dat <- data.table::fread(paste0('./grid2largeB/', sim_idx, '.csv'))
  
  d_z <- (ncol(dat) - 3)/2
  
  # Subset the data
  tmp_with_gamma_theta <- dat[sub_idx:(sub_idx + 9999)]
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
  
  ace_2sls <- broom::tidy(f1)$estimate[2]
  
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
  
  # Run leakyIV
  
  x_dat <- data.matrix(dplyr::select(tmp, x), rownames.force = NA)
  y_dat <- data.matrix(dplyr::select(tmp, y), rownames.force = NA)
  z_dat <- dplyr::select(tmp, -c(x, y, x_hat))
  
  #print(is.numeric(x_dat))
  
  #print(z_dat)
  
  tau_l2 <- LambertW::lp_norm(as.numeric(gamma[1]), 2)
  
  leaky_ace_0 = leakyIV(x_dat, y_dat, z_dat, tau_l2)
  
  leaky_ace_0_lo = dplyr::select(leaky_ace_0, 1)
  leaky_ace_0_hi = dplyr::select(leaky_ace_0, 2)
  
  leaky_ace_1 = leakyIV(x_dat, y_dat, z_dat, 1.01*tau_l2)
  
  leaky_ace_1_lo = dplyr::select(leaky_ace_1, 1)
  leaky_ace_1_hi = dplyr::select(leaky_ace_1, 2)
  
  leaky_ace_2 = leakyIV(x_dat, y_dat, z_dat, 1.1*tau_l2)
  
  leaky_ace_2_lo = dplyr::select(leaky_ace_2, 1)
  leaky_ace_2_hi = dplyr::select(leaky_ace_2, 2)
  
  leaky_ace_3 = leakyIV(x_dat, y_dat, z_dat, 2*tau_l2)
  
  leaky_ace_3_lo = dplyr::select(leaky_ace_3, 1)
  leaky_ace_3_hi = dplyr::select(leaky_ace_3, 2)
  
  leaky_ace_4 = leakyIV(x_dat, y_dat, z_dat, 5*tau_l2)
  
  leaky_ace_4_lo = dplyr::select(leaky_ace_4, 1)
  leaky_ace_4_hi = dplyr::select(leaky_ace_4, 2)
  
  leaky_ace_5 = leakyIV(x_dat, y_dat, z_dat, 10*tau_l2)
  
  leaky_ace_5_lo = dplyr::select(leaky_ace_5, 1)
  leaky_ace_5_hi = dplyr::select(leaky_ace_5, 2)
  
  leaky_ace_6 = leakyIV(x_dat, y_dat, z_dat, 10*tau_l2)
  
  leaky_ace_6_lo = dplyr::select(leaky_ace_6, 1)
  leaky_ace_6_hi = dplyr::select(leaky_ace_6, 2)
  
  print(paste0('sim_idx: ', toString(sim_idx), '  sub_idx: ', toString(sub_idx)))
  
  # Export
  out <- data.frame(
    sim_idx = sim_idx, sub_idx = sub_idx, 
    'TSLS' = ace_2sls, 'sisVIVE' = ace_sisVIVE, 'MBE' = ace_mbe,
    'D is 0 low' = leaky_ace_0_lo, 'D is 0 high' = leaky_ace_0_hi,
    'D is 0.01 low' = leaky_ace_1_lo, 'D is 0.01 high' = leaky_ace_1_hi,
    'D is 0.1 low' = leaky_ace_2_lo, 'D is 0.1 high' = leaky_ace_2_hi,
    'D is 1 low' = leaky_ace_3_lo, 'D is 1 high' = leaky_ace_3_hi,
    'D is 4 low' = leaky_ace_4_lo, 'D is 4 high' = leaky_ace_4_hi,
    'D is 9 low' = leaky_ace_5_lo, 'D is 9 high' = leaky_ace_5_hi,
    'D is 99 low' = leaky_ace_6_lo, 'D is 99 high' = leaky_ace_6_hi
    )
  return(out)
}

bnchmrk(1, 20001)

# Execute in parallel
res <- foreach(aa = seq_len(72), .combine = rbind, .packages = c("data.table", "broom", "sisVIVE", 
                                                                "tidyverse", "ggsci", "LambertW")) %:%
  foreach(bb = seq(1, 90001, 10000), .combine = rbind) %dopar%
  #prev_time <- Sys.time()
  #curr_time <- Sys.time()
  bnchmrk(aa, bb)

print(res)

fwrite(res, './grid2largeB/results.csv')

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