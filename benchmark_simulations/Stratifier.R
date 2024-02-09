# Set working directory
setwd('C:/Users/k23067841/Downloads/LeakyIV simulations again')

# Load libraries, register cores
library(data.table)
library(broom)
library(sisVIVE)
library(tidyverse)
library(ggsci)
library(LambertW)
library(dplyr)
library(tidyselect)

# Set seed
set.seed(123, kind = "L'Ecuyer-CMRG")

# Import sim_idx
sim_idx_metadata <- fread('./sisVIVE test/sim_idx.csv')

results <- fread('./sisVIVE test/results.csv')


print(results)

diag_sims <- dplyr::select(sim_idx_metadata[sim_idx_metadata$z_rho == 0], sim_idx, pr_valid)

toeplitz_sims <- dplyr::select(sim_idx_metadata[sim_idx_metadata$z_rho == 0.5], sim_idx, pr_valid)

print(diag_sims)

benchmark_estimates <- dplyr::select(results, sim_idx, tau_div_norm, TSLS:MBE, theta_true)

tight_benchmark_estimates <- benchmark_estimates[benchmark_estimates$tau_div_norm == 1]

print(tight_benchmark_estimates)

ATE_bounds <- dplyr::select(results, sim_idx, tau_div_norm, ATE_lo:ATE_hi)

tight_ATE_bounds <- ATE_bounds[ATE_bounds$tau_div_norm == 1]

print(tight_ATE_bounds)

tight_benchmark_estimates <- tight_benchmark_estimates %>%
  group_by(sim_idx) %>%
  summarize(sim_idx       = unique(sim_idx),
            TSLS_mean     = mean(TSLS),
            sisVIVE_mean  = mean(sisVIVE),
            MBE_mean      = mean(MBE),
            theta_true    = unique(theta_true) )

tight_ATE_bounds <- tight_ATE_bounds %>%
  group_by(sim_idx) %>%
  summarize(sim_idx       = unique(sim_idx),
            ATE_lo_mean   = mean(ATE_lo, na.rm=TRUE),
            ATE_hi_mean   = mean(ATE_hi, na.rm=TRUE))




tight_avg_results <- dplyr::left_join(tight_benchmark_estimates, tight_ATE_bounds, by="sim_idx")

print(tight_benchmark_estimates)
print(tight_ATE_bounds)
print(tight_avg_results)


tight_avg_results_diag <- dplyr::left_join(diag_sims, tight_avg_results, by="sim_idx")
print(tight_avg_results_diag)

tight_avg_results_toeplitz <- dplyr::left_join(toeplitz_sims, tight_avg_results, by="sim_idx")
print(tight_avg_results_toeplitz)

polygon(c(tight_avg_results_diag$pr_valid, rev(tight_avg_results_diag$pr_valid)), c(abs(tight_avg_results_diag$ATE_lo_mean), abs(rev(tight_avg_results_diag$ATE_hi_mean))),
        col = "#6BD7AF")

lines(tight_avg_results_diag$pr_valid, abs(tight_avg_results_diag$theta_true), type = "l", ylab = "y")

points(tight_avg_results_diag$pr_valid, abs(tight_avg_results_diag$MBE_mean), col = "red")
points(tight_avg_results_diag$pr_valid, abs(tight_avg_results_diag$TSLS_mean), col = "blue")
points(tight_avg_results_diag$pr_valid, abs(tight_avg_results_diag$sisVIVE), col = "yellow")

#tight_benchmark_estimates[tight_benchmark_estimates$sim_idx == 1] <- 

#tight_avg_benchmark <- 
#tight_avg_bounds <- 


plot()





