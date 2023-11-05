# Attempt at finite data leaky_IV model with certain confidence bounds 
# (up to deletion of any uninformative/weakly informative eigen-instruments)
# Assumes different samples are stored as different rows

install.packages("comprehenr")
library("comprehenr")

get_bounded_theta <- function(x, y, z, N, confidence){
  
  N <- length(x)
  
  x_bounds <- get_bounded_expectation(x, N, confidence)
  x_sq_bounds <- get_bounded_expectation(x*x, N, confidence)
  
  y_bounds <- get_bounded_expectation(y, N, confidence)
  y_sq_bounds <- get_bounded_expectation(y*y, N, confidence)
  
  xy_bounds <- get_bounded_expectation(x*y, N, confidence)
  
  zi_lbounds <- vector(for (i in 1:N) get_bounded_expectation(z[ ,i]*z[ ,j])[1], N, confidence)
  zi_ubounds <- vector(for (i in 1:N) get_bounded_expectation(z[ ,i]*z[ ,j])[2], N, confidence)

  # For each z_i, z_j find expectation values relevant in cov(z_i, z_j)
  
  zi_zj_lbounds <- matrix(for (i in 1:N) for (j in 1:N) get_bounded_expectation(z[ ,i]*z[ ,j])[1], N, confidence)
  zi_zj_ubounds <- matrix(for (i in 1:N) for (j in 1:N) get_bounded_expectation(z[ ,i]*z[ ,j])[2], N, confidence)

  zi_times_zj_lbounds <- matrix(for (i in 1:N) for (j in 1:N) get_bounded_multiplication(c(zi_lbounds[i], zi_ubounds[i]),
                                                                                         c(zi_lbounds[i], zi_ubounds[i]))[1])
  
  zi_times_zj_ubounds <- matrix(for (i in 1:N) for (j in 1:N) get_bounded_multiplication(c(zi_lbounds[i], zi_ubounds[i]),
                                                                                         c(zi_lbounds[i], zi_ubounds[i]))[2])
  
  cov_z_lbounds <- zi_zj_lbounds - zi_times_zj_ubounds
  cov_z_ubounds <- zi_zj_ubounds - zi_times_zj_lbounds
  
  
  #FACTORIALLY HARD in number of instruments?? Want bounds in eigenvalues and eigenvectors on cov_z
  #Need to grid search through all allowed cov values?
  #Could it be cheap to sample random allowed cov matrices (maybe always including extreme points) and update bounds on e_vals, e_vecs until bored?
  #Also could seed a set of gradient ascent/descent procedures to avoid the curse of dimensionality?

  
  
  
  }



get_bounded_expectation <- function(data, N, confidence){

  mean <- mean(data)
  t_score <- qt(p=confidence/2, df=N,lower.tail=F)
  sd_err <- sd(data) / sqrt(N)
  delta <- t_score * sd_err
  lower_bound <- mean - delta
  upper_bound <- mean + delta

  return(c(lower_bound, upper_bound))

  }

get_bounded_multiplication <- function(factor_1_bounds, factor_2_bounds){

  combinations <- c(factor_1_bounds[1]*factor_2_bounds[1], factor_1_bounds[1]*factor_2_bounds[2], 
                    factor_1_bounds[2]*factor_2_bounds[1], factor_1_bounds[2]*factor_2_bounds[2])

  return( c(min(combinations), max(combinations)) )
  
  }
                          
