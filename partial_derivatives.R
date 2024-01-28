# Load libraries
library(Deriv)

# Define functions
norm_mn <- function(rho) {
  gg <- (ee - var_x) * (1 + ((ee - var_x) / (eta_x^2 * rho^2)))
  hh <- -(dd - sigma_xy) * (1 + ((ee - var_x) / (eta_x^2 * rho^2)))
  ii <- ((dd - sigma_xy)^2 / (eta_x^2 * rho^2)) + ff - var_y
  # TECHNICALLY THE SOLUTION INVOLVES THIS IF/ELSE STATEMENT:
  #if (rho < 0) {
  #  theta_r <- (-hh + sqrt(hh^2 - gg * ii)) / gg
  #} else {
  #  theta_r <- (-hh - sqrt(hh^2 - gg * ii)) / gg
  #}
  # BUT I'M GUESSING THAT WON'T PLAY NICE WITH AUTODIFF, SO 
  # LET'S JUST FOLLOW ONE ROOT FOR NOW
  theta_r <- (-hh + sqrt(hh^2 - gg * ii)) / gg
  gamma_r <- Theta_z %*% (Sigma_zy - theta_r * Sigma_zx)
  (sum(abs(gamma_r)^p))^(1 / p)
}
Deriv(norm_mn)

# And the corresponding loss
loss_fn <- function(rho) {
  gg <- (ee - var_x) * (1 + ((ee - var_x) / (eta_x^2 * rho^2)))
  hh <- -(dd - sigma_xy) * (1 + ((ee - var_x) / (eta_x^2 * rho^2)))
  ii <- ((dd - sigma_xy)^2 / (eta_x^2 * rho^2)) + ff - var_y
  # TECHNICALLY THE SOLUTION INVOLVES THIS IF/ELSE STATEMENT:
  #if (rho < 0) {
  #  theta_r <- (-hh + sqrt(hh^2 - gg * ii)) / gg
  #} else {
  #  theta_r <- (-hh - sqrt(hh^2 - gg * ii)) / gg
  #}
  # BUT I'M GUESSING THAT WON'T PLAY NICE WITH AUTODIFF, SO 
  # LET'S JUST FOLLOW ONE ROOT FOR NOW
  theta_r <- (-hh + sqrt(hh^2 - gg * ii)) / gg
  gamma_r <- Theta_z %*% (Sigma_zy - theta_r * Sigma_zx)
  norm_r <- (sum(abs(gamma_r)^p))^(1 / p)
  (tau - norm_r)^2
}
Deriv(loss_fn)

# Then we set these partial derivatives to zero and solve for rho