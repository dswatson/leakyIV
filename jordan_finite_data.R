# Attempt at finite data leaky_IV model with certain confidence bounds 
# (up to deletion of any uninformative/weakly informative eigen-instruments)
# Assumes different samples are stored as different rows

library(comprehenr)

get_bounded_theta <- function(x, y, z, N, confidence_level){
  
  N <- length(x)
  
  x_bounds <- get_bounded_expectation(x, N)
  x_sq_bounds <- get_bounded_expectation(x*x, N)
  
  y_bounds <- get_bounded_expectation(y, N)
  y_sq_bounds <- get_bounded_expectation(y*y, N)
  
  xy_bounds <- get_bounded_expectation(x*y, N)
  
  cov_z_lbounds <- matrix(0, N, N)
  cov_z_ubounds <- matrix(0, N, N)
  
  zi_lbounds <- vector(for (i in 1:N) get_bounded_expectation(z[ ,i]*z[ ,j])[1])
  zi_ubounds <- vector(for (i in 1:N) get_bounded_expectation(z[ ,i]*z[ ,j])[2])
  
  # For each z_i, z_j find expectation values relevant in cov(z_i, z_j)
  
  zi_zj_lbounds <- matrix(for (i in 1:N) for (j in 1:N) get_bounded_expectation(z[ ,i]*z[ ,j])[1], N)
  zi_zj_ubounds <- matrix(for (i in 1:N) for (j in 1:N) get_bounded_expectation(z[ ,i]*z[ ,j])[2], N)
  
  
  
  
  #FACTORIALLY HARD in number of instruments?? Want bounds in eigenvalues and eigenvectors on cov_z -> 
  #Need to grid search through all allowed cov values?
  #Could be cheap to sample random allowed cov matrices (maybe always including extreme points) and update bounds on e_vals, e_vecs until bored?
  
  

  }



get_bounded_expectation <- function(data, N){

  mean <- mean(data)
  t_score <- qt(p=alpha/2, df=N,lower.tail=F)
  sd_err <- sd(data) / sqrt(N)
  delta <- t_score * sd_err
  lower_bound <- mean - delta
  upper_bound <- mean + delta

  return(c(lower_bound, upper_bound))
  
  }
