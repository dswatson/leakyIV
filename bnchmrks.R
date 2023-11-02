# Set working directory
setwd('~/Documents/UCL/soft_instruments')

# Load libraries, register cores
library(data.table)
library(broom)
library(sisVIVE)
library(tidyverse)
library(ggsci)
library(doMC)
registerDoMC(8)

# Load MBE script
source('MBE.R')

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Import sim_idx
sim_idx <- fread('./simulations/sim_idx.csv')

# Benchmark against sisVIVE and 2SLS:
bnchmrk <- function(sim_idx, sub_idx) {
  # Read the data
  dat <- fread(paste0('./simulations/sim', sim_idx, '.csv'))
  # Subset the data
  tmp <- dat[sub_idx:(sub_idx + 999)]
  # Run 2SLS
  f0 <- lm(x ~ ., data = select(tmp, -y))
  tmp$x_hat <- fitted(f0)
  f1 <- lm(y ~ x_hat, data = tmp)
  ace_2sls <- tidy(f1)$estimate[2]
  # Run sisVIVE
  sisvive <- cv.sisVIVE(tmp$y, tmp$x, as.matrix(select(tmp, starts_with('z'))))
  ace_sisVIVE <- sisvive$beta
  # Run MBE
  beta_hat <- tidy(f0)$estimate[2:9]
  se_beta <- tidy(f0)$std.error[2:9]
  f2 <- lm(y ~ ., data = select(tmp, -x_hat))
  gamma_hat <- tidy(f2)$estimate[2:9]
  se_gamma <- tidy(f2)$std.error[2:9]
  mbe <- MBE(beta_hat, gamma_hat, se_beta, se_gamma, phi = 1, n_boot = 1)
  ace_mbe <- mbe$Estimate[2]
  # Export
  out <- data.frame(
    sim_idx = sim_idx, sub_idx = sub_idx, 
    'TSLS' = ace_2sls, 'sisVIVE' = ace_sisVIVE, 'MBE' = ace_mbe
  )
  return(out)
}

# Execute in parallel
res <- foreach(aa = seq_len(24), .combine = rbind) %:%
  foreach(bb = seq(1, 9001, 1000), .combine = rbind) %dopar%
  bnchmrk(aa, bb)

# Plot results
df <- res %>%
  select(-sub_idx) %>%
  rename(idx = sim_idx) %>%
  pivot_longer(TSLS:MBE, names_to = 'Method', values_to = 'alpha') %>%
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


df <- fread('./simulations/experiment_simulations_3.csv')[1:24]
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












