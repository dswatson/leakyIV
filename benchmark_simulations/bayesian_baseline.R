################################################################################
# baseline.R
#
# A heuristic baseline approach that does Bayesian inference with Gaussian
# likelihoods and return bounds based on tail intervals. Here, we assume that
# instruments Z are mutually independent to speed things up.

library(MASS)
library(matrixsampling)
library(progress)

################################################################################
# baseline_gaussian_sample:
#
# Inputs:
#
# - m: number of MCMC samples
# - dat: data
# - pos_z, pos_x, pos_y: positions of instruments, treatment and outcome
#                        within dat
# - lvar_z_mu, lvar_z_var: hyper prior (log normal) for variances of Z
# - beta_var, pre_gamma_var: prior variance for coefficients beta and gamma
# - theta_var: prior variance for coefficient theta
# - lvar_error_mu, lvar_error_var: hyper prior (log normal) for error variances
# - tau: bound on the L2 norm of gamma
# - alpha: probability level for intervals

baseline_gaussian_sample <- function(m, dat, pos_z, pos_x, pos_y, prop_var,
                                     lvar_z_mu, lvar_z_var,
                                     beta_var, pre_gamma_var, theta_var,
                                     lvar_error_mu, lvar_error_var, tau,
                                     alpha)
{
  
  ##############################################################################
  # Basic information
  
  n <- nrow(dat)
  S <- t(dat[, c(pos_z, pos_x, pos_y)]) %*% dat[, c(pos_z, pos_x, pos_y)]
  num_z <- length(pos_z)
  
  prior_info <- list(lvar_z_mu=lvar_z_mu, lvar_z_var=lvar_z_var,
                     beta_var=beta_var, pre_gamma_var=pre_gamma_var, 
                     theta_var=theta_var,
                     lvar_error_mu=lvar_error_mu, lvar_error_var=lvar_error_var)
  
  ##############################################################################
  # Initialization
  
  lvar_z <- array(0, dim=c(num_z, m))
  beta_c <- array(0, dim=c(num_z, m))
  pre_gamma_c <- array(0, dim=c(num_z, m))
  leta_x2 <- rep(0, m)
  leta_y2 <- rep(0, m)
  zrho <- rep(0, m)
  theta <- rep(0, m)
  zkappa <- rep(0, m)
  
  lvar_z[, 1] <- rnorm(num_z) / 10
  beta_c[, 1] <- rnorm(num_z) / 10
  pre_gamma_c[, 1] <- rnorm(num_z) / 10
  leta_x2[1] <- rnorm(1) / 10
  leta_y2[1] <- rnorm(1) / 10
  zrho[1] <- rnorm(1) / 10
  theta[1] <- rnorm(1) / 10
  zkappa[1] <- rnorm(1) / 10
  
  ##############################################################################
  # Main loop
  
  pb <- progress_bar$new(total = m)
  
  for (i in 2:m) {
    
    pb$tick()
    
    params <- list(lvar_z=lvar_z[, i - 1], beta_c=beta_c[, i - 1], 
                   pre_gamma_c=pre_gamma_c[, i - 1], leta_x2=leta_x2[i - 1],
                   leta_y2=leta_y2[i - 1], zrho=zrho[i - 1], 
                   theta=theta[i - 1], zkappa=zkappa[i - 1])
    
    # Variances of Z
    for (j in 1:num_z) 
      params$lvar_z[j] <- do_mh_step_baseline(n, S, num_z, prop_var, params, prior_info, tau, j, 1)
    
    # beta coefficients  
    for (j in 1:num_z) 
      params$beta_c[j] <- do_mh_step_baseline(n, S, num_z, prop_var, params, prior_info, tau, j, 2)
    
    # pre-gamma coefficients  
    #for (j in 1:num_z) 
      params$pre_gamma_c[j] <- do_mh_step_baseline(n, S, num_z, prop_var, params, prior_info, tau, j, 3)
    
    # Other
    params$leta_x2 <- do_mh_step_baseline(n, S, num_z, prop_var, params, prior_info, tau, j, 4)
    params$leta_y2 <- do_mh_step_baseline(n, S, num_z, prop_var, params, prior_info, tau, j, 5)
    params$zrho <- do_mh_step_baseline(n, S, num_z, prop_var, params, prior_info, tau, j, 6)
    params$theta <- do_mh_step_baseline(n, S, num_z, prop_var, params, prior_info, tau, j, 7)
    params$zkappa <- do_mh_step_baseline(n, S, num_z, prop_var, params, prior_info, tau, j, 8)
    
    # Move on
    lvar_z[, i] <- params$lvar_z
    beta_c[, i] <- params$beta_c
    pre_gamma_c[, i] <- params$pre_gamma_c
    leta_x2[i] <- params$leta_x2
    leta_y2[i] <- params$leta_y2
    zrho[i] <- params$zrho
    theta[i] <- params$theta
    zkappa[i] <- params$zkappa
    
  }
  
  ##############################################################################
  # Wrap-up
  
  sort_theta <- sort(theta)
  alpha_idx <- round(m * alpha)
  lb_interval <- c(-Inf, sort_theta[alpha_idx])
  ub_interval <- c(sort_theta[m - alpha_idx + 1], Inf)
  
  result <- list(lvar_z=lvar_z, beta_c=beta_c, pre_gamma_c=pre_gamma_c, 
                 leta_x2=leta_x2, leta_y2=leta_y2, zrho=zrho, 
                 theta=theta, zkappa=zkappa, 
                 lb_interval=lb_interval, ub_interval=ub_interval)
  return(result)
  
}

################################################################################
# baseline_gaussian_model_function:
#
# The log-likelihood + log-prior (up to an additive constant) of a model
# where Z are independently distributed with their own variances,
# X = Z * beta + e_x, Y = Z * gamma + X * theta + e_y.
#
# Inputs:
# Inputs:
#
# - n, S: sufficient statistics
# - num_z: number of instruments. We assume that the first num_z columns of
#          S are the instrument, num_z + 1 is the position of the treatment and
#          num_z + 2 is the position of the outcome
# - lvar_z: log vector of variances of Z
# - bet: coefficients beta
# - pre_g: pre-normalized coefficients gamma
# - zkappa: unconstrained representation of redundant parameter to define gamma
# - th: coefficient theta
# - lve_x2, lve_y2: log variances of error terms e_x and e_y
# - zrho: unconstrained representation of correlation coefficient
# - lvar_z_mu, lvar_z_var: hyper prior (log normal) for variances of Z
# - beta_var, pre_gamma_var: prior variance for coefficients beta and gamma
# - theta_var: prior variance for coefficient theta
# - lvar_error_mu, lvar_error_var: hyper prior (log normal) for error variances
# - tau: bound on the L2 norm of gamma
# - new_v, new_v_pos, type_param: do we need to change a parameter? Where and
#                                 by which value?

baseline_gaussian_model_function <- function(n, S, num_z,
                                     params, prior_info, tau,
                                     new_v=NA, new_v_pos=NA, type_param=NA)
{
  
  if (!is.na(new_v)) {
    if (type_param == 1) {
      params$lvar_z[new_v_pos] <- new_v
    } else if (type_param == 2) {
      params$beta_c[new_v_pos] <- new_v
    } else if (type_param == 3) {
      params$pre_gamma_c[new_v_pos] <- new_v
    } else if (type_param == 4) {
      params$leta_x2 <- new_v
    } else if (type_param == 5) {
      params$leta_y2 <- new_v
    } else if (type_param == 6) {
      params$zrho <- new_v
    } else if (type_param == 7) {
      params$theta <- new_v
    } else if (type_param == 8) {
      params$zkappa <- new_v
    }
  }
  
  var_z <- exp(params$lvar_z)
  kappa <- pnorm(params$zkappa)
  gamma_c <- params$pre_gamma_c * sqrt(kappa * tau / sum(params$pre_gamma_c^2))
  eta_x2 <- exp(params$leta_x2)
  eta_y2 <- exp(params$leta_y2)
  rho <- pnorm(params$zrho) * 2 - 1
  eta_xy <- sqrt(eta_x2 * eta_y2) * rho
  
  Sigma <- matrix(0, num_z + 2, num_z + 2)
  
  # Covariance of Z
  for (i in 1:num_z)
    Sigma[i, i] <- var_z[i]
  
  # Cross-covariance (Z, X)
  Sigma[1:num_z, num_z + 1] <- var_z * params$beta_c
  Sigma[num_z + 1, 1:num_z] <- Sigma[1:num_z, num_z + 1]
  
  # Variance of X
  Sigma[num_z + 1, num_z + 1] <- sum(params$beta_c^2 * var_z) + eta_x2
  
  # Cross-covariance (Z, Y)
  Sigma[1:num_z, num_z + 2] <- var_z * gamma_c + 
                               Sigma[1:num_z, num_z + 1] * params$theta
  Sigma[num_z + 2, 1:num_z] <- Sigma[1:num_z, num_z + 2]
    
  # Cross-covariance (X, Y)
  Sigma[num_z + 1, num_z + 2] <- Sigma[num_z + 1, 1:num_z] %*% gamma_c + 
                                 eta_x2 * params$theta + eta_xy
  Sigma[num_z + 2, num_z + 1] <- Sigma[num_z + 1, num_z + 2]
  
  # Variance Y
  Sigma[num_z + 2, num_z + 2] <- sum(gamma_c^2 * var_z) + 
                                 2 * params$theta * sum(gamma_c * Sigma[1:num_z, num_z + 1]) +
                                 params$theta^2 * Sigma[num_z + 1, num_z + 1] +
                                 2 * eta_xy + eta_y2

  # Log-likelihood
  llik <- -0.5 * (n * determinant(Sigma)$modulus[1] + sum(diag(solve(Sigma, S))))

  # Log-prior
  lprior <- -0.5 * (sum((params$lvar_z - prior_info$lvar_z_mu)^2) / prior_info$lvar_z_var +
                    sum(params$beta_c^2) / prior_info$beta_var +
                    sum(params$pre_gamma_c^2) / prior_info$pre_gamma_var +
                    params$zkappa^2 +
                    params$theta^2 / prior_info$theta_var + 
                    (params$leta_x2 - prior_info$lvar_error_mu)^2 / prior_info$lvar_error_var +
                    (params$leta_y2 - prior_info$lvar_error_mu)^2 / prior_info$lvar_error_var +
                    params$zrho^2)

  # Result
  return(llik + lprior)
  
}

get_baseline_value <- function(params, new_v_pos, type_param)
{
  if (type_param == 1) {
    return(params$lvar_z[new_v_pos])
  } else if (type_param == 2) {
    return(params$beta_c[new_v_pos])
  } else if (type_param == 3) {
    return(params$pre_gamma_c[new_v_pos])
  } else if (type_param == 4) {
    return(params$leta_x2)
  } else if (type_param == 5) {
    return(params$leta_y2)
  } else if (type_param == 6) {
    return(params$zrho)
  } else if (type_param == 7) {
    return(params$theta)
  } else if (type_param == 8) {
    return(params$zkappa)
  }
  return(NA)
}

do_mh_step_baseline <- function(n, S, num_z, prop_var,
                                params, prior_info, tau,
                                new_v_pos, type_param)
{
  
  current <- baseline_gaussian_model_function(n, S, num_z, params, prior_info, tau)
  old_v <- get_baseline_value(params, new_v_pos, type_param)
  new_v <- old_v + rnorm(1) * sqrt(prop_var)
  attempt <- baseline_gaussian_model_function(n, S, num_z, params, prior_info, tau,
                                              new_v, new_v_pos, type_param)
  if (runif(1) <= exp(attempt - current))
    return(new_v)
  
  return(old_v)

}

################################################################################
# debug_baseline:
#
# Totally artificial example just for a reality check.

debug_baseline <- function()
{
  # Model
  
  num_z <- 5
  beta_true <- rnorm(num_z) / sqrt(num_z)
  gamma_true <- rnorm(num_z) / sqrt(3 * num_z)
  theta_true <- 1.
  e_x2_true <- 0.2
  e_y2_true <- 0.2
  rho_true <- 0.5
  var_z_true <- rep(1, num_z)
  
  tau <- sum(gamma_true^2)
  
  # Data

  n <- 1000
  
  sigma_error <- matrix(c(e_x2_true, rho_true * sqrt(e_x2_true * e_y2_true), 
                          rho_true * sqrt(e_x2_true * e_y2_true), e_y2_true), ncol=2)
  errors <- mvrnorm(n, c(0, 0), sigma_error)
  
  z <- matrix(nrow=n, ncol=num_z)
  for (i in 1:num_z)
    z[, i] <- rnorm(n) * sqrt(var_z_true[i])
  x <- z %*% beta_true + errors[, 1]
  y <- z %*% gamma_true + x * theta_true + errors[, 2]
  dat <- cbind(z, x, y)
  pos_z <- 1:num_z
  pos_x <- num_z + 1
  pos_y <- num_z + 2
    
  # Prior
  
  lvar_z_mu <- -1
  lvar_z_var <- 3
  beta_var <- 5
  pre_gamma_var <- 1
  theta_var <- 5
  lvar_error_mu <- -1
  lvar_error_var <- 1
  
  # Call MCMC
  
  m <- 5000
  prop_var <- 0.01
  alpha <- 0.05
  result <- baseline_gaussian_sample(m, dat, pos_z, pos_x, pos_y, prop_var, 
                                     lvar_z_mu, lvar_z_var,
                                     beta_var, pre_gamma_var, theta_var,
                                     lvar_error_mu, lvar_error_var, tau, alpha)
}
