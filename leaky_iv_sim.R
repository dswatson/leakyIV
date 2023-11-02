# Set working directory
setwd('~/Documents/UCL/soft_instruments')

# Load libraries, register cores
library(data.table)
library(matrixStats)
library(doRNG)
library(doMC)
registerDoMC(10)

# Simulation script
source('simulation.R')

# Leaky IV script
source('leaky_iv.R')

# Set seed
set.seed(123, "L'Ecuyer-CMRG")

# Loop over different data configurations
loop_fn <- function(b, theta_b, z_rhob, rho_b, r2_xb, r2_yb, prop_b) {
  tmp <- sim_dat(
    n = 1e4, d_z = 10, z_cnt = TRUE, z_rho = z_rhob, rho = rho_b, 
    theta = theta_b, r2_x = r2_xb, r2_y = r2_yb, pr_valid = prop_b, idx_b
  )
  out <- foreach(k = 1:10, .combine = rbind) %do% {
    dat <- tmp$dat
  }
  
  
  # Assume oracle access to tau
  out <- leaky_iv(tmp$dat$z, tmp$dat$x, tmp$dat$y, sum(tmp$params$gamma^2),
                  parallel = FALSE)
  out[, ACE := theta_b][, rho := rho_b][, r2_x := r2_xb][, r2_y := r2_yb][, pr_valid := prop_b]
  return(out)
}
# Execute in parallel
df <- foreach(thetas = c(-2, -1, -1/2, -1/4, 1/4, 1/2, 1, 2), .combine = rbind) %:%
  foreach(z_rhos = c(0, 1/4, 1/2, 3/4), .combine = rbind) %:%
  foreach(rhos = c(-3/4, -1/2, -1/4, 1/4, 1/2, 3/4), .combine = rbind) %:%
  foreach(r2_xs = c(1/4, 1/2, 3/4), .combine = rbind) %:%
  foreach(r2_ys = c(1/4, 1/2, 3/4), .combine = rbind) %:%
  foreach(props = c(0, 0.1, 0.2, 0.3, 0.4, 0.5), .combine = rbind) %dopar%
  loop_fn(theta_b = thetas, z_rhob = z_rhos, rho_b = rhos, r2_xb = r2_xs, r2_yb = r2_ys, 
          prop_b = props)

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


