# Load libraries
library(Deriv)

# Define functions
norm_mn <- function(rho) {
  aa <- Sigma_xz %*% Theta_z %*% Sigma_zx - var_x
  bb <- Sigma_xz %*% Theta_z %*% Sigma_zy - sigma_xy
  cc <- Sigma_yz %*% Theta_z %*% Sigma_zy - var_y
  theta_r <- bb / aa - sign(rho) * sqrt((bb^2 - aa * cc) * (1 + aa / (eta_x^2 * rho^2))) / 
    (aa * (1 + aa / (eta_x^2 * rho^2)))
  gamma_r <- Theta_z %*% (Sigma_zy - theta_r * Sigma_zx)
  (sum(abs(gamma_r)^p))^(1 / p)
}
Deriv(norm_mn)

# And the corresponding loss
loss_fn <- function(rho) {
  aa <- Sigma_xz %*% Theta_z %*% Sigma_zx - var_x
  bb <- Sigma_xz %*% Theta_z %*% Sigma_zy - sigma_xy
  cc <- Sigma_yz %*% Theta_z %*% Sigma_zy - var_y
  theta_r <- bb / aa - sign(rho) * sqrt((bb^2 - aa * cc) * (1 + aa / (eta_x^2 * rho^2))) / 
    (aa * (1 + aa / (eta_x^2 * rho^2)))
  gamma_r <- Theta_z %*% (Sigma_zy - theta_r * Sigma_zx)
  norm_r <- (sum(abs(gamma_r)^p))^(1 / p)
  (tau - norm_r)^2
}
Deriv(loss_fn)

# Then we set these derivatives to zero and solve for rho


