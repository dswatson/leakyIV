# Load libraries
library(Deriv)

# Define functions
norm_mn <- function(rho) {
  theta <- (psi / eta_x2) - sign(rho) * 
    ((sqrt((1 - 1/rho^2) * (psi^2 - phi2 * eta_x2)))) / 
    (-eta_x2 * (1 - 1/rho^2))
  gamma <- Theta_z %*% (Sigma_zy - theta * Sigma_zx)
  (sum(abs(gamma)^p))^(1 / p)
}
Deriv(norm_mn)

# And the corresponding loss
loss_fn <- function(rho) {
  theta <- bb / aa - sign(rho) * sqrt((bb^2 - aa * cc) * (1 + aa / (-aa * rho^2))) / 
    (aa * (1 + aa / (-aa * rho^2)))
  gamma <- Theta_z %*% (Sigma_zy - theta * Sigma_zx)
  norm <- (sum(abs(gamma)^p))^(1 / p)
  (tau - norm)^2
}
Deriv(loss_fn)

# Then we set these derivatives to zero and solve for rho!


