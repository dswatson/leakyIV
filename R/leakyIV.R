#' Bounding Causal Effects with Leaky Instruments
#' 
#' Estimates bounds on average treatment effects in linear IV models under 
#' limited violations of the exclusion criterion.
#' 
#' @param dat Input data. Either (a) an \eqn{n \times d} data frame or matrix of 
#'   observations with columns for treatment, outcome, and candidate instruments; 
#'   or (b) a \eqn{d \times d} covariance matrix over such variables. The latter 
#'   is incompatible with bootstrapping. Note that in either case, the order of
#'   variables is presumed to be treatment (\eqn{X}), outcome (\eqn{Y}), leaky 
#'   instruments (\eqn{Z}).
#' @param tau Either (a) a scalar representing the upper bound on the p-norm of 
#'   linear weights on \eqn{Z} in the structural equation for \eqn{Y}; or (b) a 
#'   vector representing upper bounds on the absolute value of each such 
#'   coefficient. See details.
#' @param p Power of the norm for the \code{tau} threshold.
#' @param normalize Scale candidate instruments to unit variance? 
#' @param method Estimator for the covariance matrix, if one is not supplied by
#'   \code{dat}. Options include (a) \code{"mle"}, the default; (b) \code{"shrink"}, 
#'   an analytic empirical Bayes solution; or (c) \code{"glasso"}, the graphical 
#'   lasso. See details.
#' @param approx Use nearest positive definite approximation if the estimated 
#'   covariance matrix is singular? See details.
#' @param n_boot Optional number of bootstrap replicates.
#' @param bayes Use Bayesian bootstrap? 
#' @param parallel Compute bootstrap estimates in parallel? Must register 
#'   backend beforehand, e.g. via \code{doParallel}.
#' @param ... Extra arguments to be passed to graphical lasso estimator if
#'   \code{method = "glasso"}. Note that the regularization parameter \code{rho}
#'   is required as input, with no default. 
#' 
#' 
#' @details 
#' Instrumental variables are defined by three structural assumptions: they must
#' be (A1) \emph{relevant}, i.e. associated with the treatment; (A2) 
#' \emph{unconfounded}, i.e. independent of common causes between treatment and 
#' outcome; and (A3) \emph{exclusive}, i.e. only affect outcomes through the 
#' treatment. The \code{leakyIV} algorithm (Watson et al., 2024) relaxes (A3), 
#' allowing some information leakage from IVs \eqn{Z} to outcomes \eqn{Y} in 
#' linear systems. While the average treatment effect (ATE) is no longer 
#' identifiable in this setting, sharp bounds can be computed exactly. 
#' 
#' We assume the following structural equation for the treatment: 
#' \eqn{X := Z \beta + \epsilon_X}, where the final summand is a noise term that
#' correlates with the additive noise in the structural equation for the outcome: 
#' \eqn{Y := Z \gamma + X \theta + \epsilon_Y}. The ATE is given by the 
#' parameter \eqn{\theta}. Whereas classical IV models require each \eqn{\gamma} 
#' coefficient to be zero, we permit some direct signal from \eqn{Z} to 
#' \eqn{Y}. Specifically, \code{leakyIV} provides support for two types of
#' information leakage: (a) thresholding the \emph{p}-norm of linear weights 
#' \eqn{\gamma} (scalar \code{tau}); and (b) thresholding the absolute value of 
#' each \eqn{\gamma} coefficient one by one (vector \code{tau}). 
#' 
#' Numerous methods exist for estimating covariance matrices. \code{leakyIV}
#' provides support for maximum likelihood estimation (the default), as well as
#' empirical Bayes shrinkage via \code{corpcor::\link[corpcor]{cov.shrink}} 
#' (Schäfer & Strimmer, 2005) and the graphical lasso via 
#' \code{glasso::\link[glasso]{glasso}} (Friedman et al., 2007). These latter 
#' methods are preferable in high-dimensional settings where sample covariance 
#' matrices may be unstable or singular. Alternatively, users can pass a 
#' pre-computed covariance matrix directly as \code{dat}.
#' 
#' Estimated covariance matrices may be singular for some datasets or bootstrap
#' samples. Behavior in this case is determined by the \code{approx} argument. 
#' If \code{TRUE}, \code{leakyIV} proceeds with the nearest positive definite 
#' approximation, computed via Higham's (2002) algorithm (with a warning). If 
#' \code{FALSE}, bounds are NA (also with a warning).
#' 
#' Uncertainty can be evaluated in leaky IV models using the bootstrap, provided
#' that covariances are estimated internally and not passed directly. 
#' Bootstrapping provides a nonparametric sampling distribution for min/max 
#' values of the ATE. Set \code{bayes = TRUE} to replace the classical bootstrap 
#' with a Bayesian bootstrap for approximate posterior inference (Rubin, 1981).
#' 
#' 
#' @return  
#' A data frame with columns for \code{ATE_lo} and \code{ATE_hi}, representing
#' lower and upper bounds of the partial identification interval for the 
#' causal effect of \eqn{X} on \eqn{Y}. When bootstrapping, the output data 
#' frame contains \code{n_boot} rows, one for each bootstrap replicate. 
#' 
#' 
#' @references  
#' Watson, D., Penn, J., Gunderson, L., Bravo-Hermsdorff, G., Mastouri, A., and
#' Silva, R. (2024). Bounding causal effects with leaky instruments. \emph{arXiv}
#' preprint, 2404.04446.
#' 
#' Friedman, J., Hastie, T., and Tibshirani, R. (2007). Sparse inverse 
#' covariance estimation with the lasso. \emph{Biostatistics}, 9:432-441.
#' 
#' Schäfer, J., and Strimmer, K. (2005). A shrinkage approach to large-scale 
#' covariance estimation and implications for functional genomics. 
#' \emph{Statist. Appl. Genet. Mol. Biol.}, 4:32.
#' 
#' Higham, N. (2002). Computing the nearest correlation matrix: A problem from 
#' finance. \emph{IMA J. Numer. Anal.}, 22:329–343.
#' 
#' Rubin, D.R. (1981). The Bayesian bootstrap. \emph{Ann. Statist.}, 
#' \emph{9}(1): 130-134. 
#' 
#' 
#' @examples  
#' set.seed(123)
#' 
#' # Hyperparameters
#' n <- 200
#' d_z <- 4
#' beta <- rep(1, d_z)
#' gamma <- rep(0.1, d_z)
#' theta <- 2
#' rho <- 0.5
#' 
#' # Simulate correlated residuals
#' S_eps <- matrix(c(1, rho, rho, 1), ncol = 2)
#' eps <- matrix(rnorm(n * 2), ncol = 2)
#' eps <- eps %*% chol(S_eps)
#' 
#' # Simulate observables from a leaky IV model
#' z <- matrix(rnorm(n * d_z), ncol = d_z)
#' x <- z %*% beta + eps[, 1]
#' y <- z %*% gamma + x * theta + eps[, 2]
#' obs <- cbind(x, y, z)
#' 
#' # Run the algorithm
#' leakyIV(obs, tau = 1)
#' 
#' # With bootstrapping
#' leakyIV(obs, tau = 1, n_boot = 10)
#' 
#' # With covariance matrix input
#' S <- cov(obs)
#' leakyIV(S, tau = 1)
#' 
#' 
#' @export 
#' @import data.table
#' @importFrom corpcor is.positive.definite cov.shrink 
#' @importFrom glasso glasso
#' @importFrom stats rexp cov.wt optim lm coef residuals
#' @importFrom Matrix nearPD
#' @importFrom foreach foreach %do% %dopar%

leakyIV <- function(
    dat,
    tau, 
    p = 2, 
    normalize = TRUE,
    method = "mle",
    approx = TRUE,
    n_boot = NULL, 
    bayes = FALSE, 
    parallel = TRUE, 
    ...) {
  
  # To avoid data.table check issues
  bb <- rho <- sat <- NULL
  
  # Checks, warnings
  if (nrow(dat) == ncol(dat)) {
    Sigma <- dat
    if (!is.positive.definite(Sigma)) {
      stop('Pre-computed covariance matrix must be positive definite.')
    }
    if (!is.null(n_boot)) {
      if (n_boot > 0L) {
        warning('Bootstrapping cannot be performed with covariance matrix input. ',
                'Setting n_boot = 0.')
        n_boot <- 0L
      }
    }
    Sigma_input <- TRUE
  } else {
    dat <- as.data.frame(dat)
    n <- nrow(dat)
    w <- rep(1L, n)
    Sigma_input <- FALSE
  }
  if (!length(tau) %in% c(1L, ncol(dat) - 2L)) {
    stop('tau must be either a scalar or a vector of length ncol(dat) - 2.')
  }
  if (any(tau < 0)) {
    stop('tau must be >= 0.')
  }
  if (max(tau) == 0) {
    stop('Some tau must be positive or else there is no information leakage. ',
         'Consider using 2SLS instead.')
  }
  if (p < 0) {
    stop('p must be >= 0.')
  }
  if (p > 0 & p < 1) {
    warning('Norms with p < 1 are improper and non-convex. ',
            'Interpret results with caution.')
  }
  if (!isTRUE(Sigma_input) & !method %in% c('mle', 'shrink', 'glasso')) {
    stop('method not recognized. Must be one of "mle", "shrink", or "glasso".')
  }
  if (is.null(n_boot)) {
    n_boot <- 0L
    parallel <- FALSE
  }
  
  # Optionally prepare transition matrix
  d <- ncol(dat)
  d_z <- d - 2L
  s0 <- NULL
  t_matrix <- NULL
  if (length(tau) > 1L) {
    if (any(tau == 0)) {
      # Partition z into valid and leaky instruments
      s0 <- which(tau == 0)
      s1 <- which(tau != 0)
      tau[s0] <- 1
    }
    tau <- c(1, 1, tau)
    t_matrix <- matrix(nrow = d, ncol = d)
    diag(t_matrix) <- 1 / tau^2
    for (i in 2:d) {
      for (j in 1:(i - 1L)) {
        t_matrix[i, j] <- t_matrix[j, i] <- 1 / (tau[i] * tau[j])
      }
    }
    tau <- 1L
    p <- Inf
  }
  
  # Bootstrap samples or weights
  if (n_boot > 0L) {
    if (isTRUE(bayes)) {
      # Draw Dirichlet weights
      w_mat <- matrix(stats::rexp(n * n_boot), ncol = n_boot)
      w_mat <- (w_mat / rowSums(w_mat)) * n
    } else {
      # Draw bootstrap samples
      b_mat <- matrix(sample(n, n * n_boot, replace = TRUE), ncol = n_boot)
    }
  }

  # Compute bounds
  loop <- function(b) {
    
    # Estimate covariance parameters
    if (b > 0L) {
      if (isTRUE(bayes)) {
        w <- w_mat[, b]
      } else {
        dat <- dat[b_mat[, b], ]
      }
    }
    if (!isTRUE(Sigma_input)) {
      if (method == 'mle') {
        Sigma <- stats::cov.wt(dat, wt = w)$cov
      } else if (method == 'shrink') {
        Sigma <- cov.shrink(dat, w = w, verbose = FALSE)[seq_len(d), seq_len(d)]
        #Theta_z <- invcov.shrink(dat, verbose = FALSE)[seq_len(d_z), seq_len(d_z)]
      } else if (method == 'glasso') {
        s <- glasso(stats::cov.wt(dat, wt = w)$cov, ...)
        Sigma <- s$w
        #Theta_z <- s$wi[seq_len(d_z), seq_len(d_z)]
      }
    }
    # Optionally normalize 
    if (isTRUE(normalize)) {
      fctr <- c(1, 1, 1/sqrt(diag(Sigma[3:d, 3:d])))
      f_matrix <- matrix(nrow = d, ncol = d)
      for (i in seq_len(d)) {
        for (j in 1:i) {
          f_matrix[i, j] <- f_matrix[j, i] <- fctr[i] * fctr[j]
        }
      }
      Sigma <- Sigma * f_matrix
    }
    if (!is.null(t_matrix)) { # For vector tau
      Sigma <- Sigma * t_matrix
    }
    # Check for singularity
    if (!is.positive.definite(Sigma) & !isTRUE(approx)) {
      warning('Covariance estimator results in a singular matrix. ',
              'Consider rerunning with another method or setting approx = TRUE.')
      ATE_lo <- ATE_hi <- NA_real_
    } else {
      if (!is.positive.definite(Sigma) & isTRUE(approx)) {
        warning('Covariance estimator results in a singular matrix. ',
                'Substituting nearest positive definite approximation.')
        Sigma <- as.matrix(nearPD(Sigma)$mat)
      }
      # Moving on...
      Theta_z <- solve(Sigma[3:d, 3:d])
      Sigma_zy <- matrix(Sigma[3:d, 2L], ncol = 1L)
      Sigma_yz <- t(Sigma_zy)
      Sigma_zx <- matrix(Sigma[3:d, 1L], ncol = 1L)
      Sigma_xz <- t(Sigma_zx)
      sigma_xy <- Sigma[1L, 2L]
      var_x <- Sigma[1L, 1L]
      var_y <- Sigma[2L, 2L]
      # Gamma is constrained to the hyperplane alpha - theta * beta
      alpha <- as.numeric(Theta_z %*% Sigma_zy)
      # alpha <- as.numeric(coef(lm(y ~ 0 + z)))
      beta <- as.numeric(Theta_z %*% Sigma_zx)
      # beta <- as.numeric(coef(lm(x ~ 0 + z)))
      if (!is.null(s0)) {
        alpha <- alpha[s1] 
        beta <- beta[s1]
      }
      
      if (p == 2) {  # Exact solution in the L2 case
        
        f <- stats::lm(alpha ~ 0 + beta)
        theta_check <- as.numeric(stats::coef(f))
        tau_check <- sqrt(sum(stats::residuals(f)^2))
        if (tau <= tau_check) {
          warning('tau is too low, resulting in an empty feasible region. ',
                  'Consider rerunning with a higher threshold.')
          ATE_lo <- ATE_hi <- NA_real_
        } else {
          delta <- as.numeric(
            sqrt(crossprod(beta) * (tau^2 - crossprod(alpha)) + (alpha %*% beta)^2) / 
              crossprod(beta)
          )
          ATE_lo <- theta_check - delta
          ATE_hi <- theta_check + delta
        }
        
      } else {       # Numerical solution otherwise
        
        # Conditional (co)variances given Z
        k_xx <- var_x - as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
        k_yy <- var_y - as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
        k_xy <- sigma_xy - as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
        if (any(c(k_xx, k_yy) < 0) | k_xx * k_yy < k_xy^2) {
          warning('Covariance estimator implies inconsistent parameters. ',
                  'Consider rerunning with another method.')
          ATE_lo <- ATE_hi <- NA_real_
        } else {     # Define some helper functions
          
          # Compute theta as a function of rho
          theta_fn <- function(rho) {
            (k_xy - sqrt(k_xx * k_yy - k_xy^2) * tan(asin(rho))) / k_xx
          }
          # Compute gamma norms as a function of rho
          norm_fn <- function(rho) {
            theta <- theta_fn(rho)
            gamma <- alpha - theta * beta
            if (p == 0L) {
              norm <- sum(gamma != 0)
            } else if (p == Inf) {
              norm <- max(abs(gamma))
            } else {
              norm <- (sum(abs(gamma)^p))^(1 / p)
            }
            return(norm)
          }
          # Compute loss as a function of rho
          loss_fn <- function(rho) {
            norm <- norm_fn(rho)
            loss <- (tau - norm)^2
            return(loss)
          }
          
          # Find the split point: upper and lower bounds lie on either side of
          # rho_check, assuming tau > tau_check = min_norm$value
          min_norm <- stats::optim(0, norm_fn, method = 'Brent', lower = -1, upper = 1)
          tau_check <- min_norm$value
          if (tau <= tau_check) {
            warning('tau is too low, resulting in an empty feasible region. ',
                    'Consider rerunning with a higher threshold.')
            ATE_lo <- ATE_hi <- NA_real_
          } else {
            rho_check <- min_norm$par
            rho_lo <- stats::optim(mean(c(-1, rho_check)), loss_fn, method = 'Brent', 
                                   lower = -1, upper = rho_check)$par
            rho_hi <- stats::optim(mean(c(rho_check, 1)), loss_fn, method = 'Brent', 
                                   lower = rho_check, upper = 1)$par
            ATE_lo <- theta_fn(rho_hi)
            ATE_hi <- theta_fn(rho_lo)
            
          }
        }
      }
    }
    
    # Export
    out <- data.table(ATE_lo, ATE_hi)
    return(out)
  }
  
  # Optionally compute in parallel
  if (n_boot == 0L) {
    out <- loop(0L)
  } else {
    if (isTRUE(parallel)) {
      out <- foreach(bb = seq_len(n_boot), .combine = rbind) %dopar% loop(bb)
    } else {
      out <- foreach(bb = seq_len(n_boot), .combine = rbind) %do% loop(bb)
    }
  }
  
  # Export
  return(out)
  
}


