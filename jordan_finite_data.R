# Attempt at finite data leaky_IV model with certain confidence bounds 
# (up to deletion of any uninformative/weakly informative eigen-instruments)

get_bounded_theta <- function(x, y, z, N, confidence_level){

  x_bounds <-
  x_sq_bounds <- 

  y_bounds <- 
  y_sq_bounds <- 

  # For each z_i, z_j find expectation values relevant in cov(z_i, z_j)
  cov_z_bounded 
  
  #FACTORIALLY HARD in number of instruments?? Want bounds in eigenvalues and eigenvectors on cov_z -> 
  #Need to grid search through all allowed cov values?
  #Could be cheap and sample random allowed cov matrices (maybe always including extreme points) and update bounds on e_vals, e_vecs until bored?
  
  inv_cov_z_bounded

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
