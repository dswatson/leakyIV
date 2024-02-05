#' Bounding Causal Effects with Leaky Instruments
#' 
#' This function estimates bounds on average treatment effects in linear IV 
#' models under limited violations of the exclusion criterion.
#' 
#' @param x Treatment variable. 
#' @param y Outcome variable.
#' @param z One or more leaky instrumental variable(s). 
#' @param tau Either (a) a scalar representing the upper bound on the p-norm of 
#'   linear weights from \code{z} to \code{y}; or (b) a vector of length 
#'   \code{ncol(z)} representing upper bounds on the absolute value of each such 
#'   coefficient. See details.
#' @param p Power of the norm on linear weights from \code{z} to \code{y} (only 
#'   relevant if \code{tau} is scalar).
#' @param n_rho Number of rho-values to sweep through in a grid search (only 
#'   relevant if \code{tau} is a vector).
#' @param method Estimator for the covariance matrix. Options include 
#'   (a) \code{"mle"}, the default; (b) \code{"shrink"}, an analytic empirical 
#'   Bayes solution; or (c) \code{"glasso"}, the graphical lasso. See details.
#' @param Sigma Optional pre-computed covariance matrix for \code{x, y, z}.
#'   If non-\code{NULL}, then \code{Sigma} overrides \code{method}. This is
#'   incompatible with bootstrapping.
#' @param n_boot Optional number of bootstrap replicates.
#' @param bayes Use Bayesian bootstrap? 
#' @param parallel Compute bootstrap estimates in parallel?
#' @param ... Extra arguments to be passed to graphical lasso estimator if
#'   \code{method = "glasso"}. Note that the regularization parameter \code{rho}
#'   is required as input, with no default. 
#' 
#' @details 
#' Instrumental variables are defined by three structural assumptions: they must
#' be (A1) \emph{relevant}, i.e. associated with the treatment; (A2) 
#' \emph{unconfounded}, i.e. independent of common causes between treatment and 
#' outcome; and (A3) \emph{exclusive}, i.e. only affect outcomes through the 
#' treatment. The \code{leakyIV} algorithm relaxes (A3), allowing some 
#' information leakage from IVs \code{z} to outcomes \code{y} in linear 
#' structural equation models. While the average treatment effect (ATE) is no 
#' longer identifiable in this setting, tight bounds can be computed exactly for 
#' scalar \code{tau} and approximately for vector \code{tau}. 
#' 
#' We assume the following structural equation for \code{x}: 
#' \eqn{X := Z \beta + \epsilon_X}, where the final summand is a noise term that
#' correlates with the additive noise in the structural equation for \code{y}: 
#' \eqn{Y := Z \gamma + X \theta + \epsilon_Y}. The ATE is given by the 
#' parameter \eqn{\theta}. Whereas classical IV models require each \eqn{gamma} 
#' coefficient to be zero, we permit some direct signal from \code{z} to 
#' \code{y}. Specifically, \code{leakyIV} provides native support for two types 
#' of information leakage: (a) thresholding the p-norm of linear weights 
#' \eqn{gamma} (scalar \code{tau}); and (b) thresholding the absolute value of 
#' each \eqn{gamma} coefficient one by one (vector \code{tau}). In the latter 
#' case, we perform a grid search over a range of \code{n_rho} candidate values 
#' for the correlation coefficient between \eqn{\epsilon_X} and \eqn{\epsilon_Y}, 
#' recording the min/max ATE consistent with assumptions and data as we vary the 
#' magnitude and direction of latent confounding.
#' 
#' Numerous methods exist for estimating covariance matrices. \code{leakyIV}
#' provides support for maximum likelihood estimation (the default), as well as
#' empirical Bayes shrinkage via \code{corpcor::\link[corpcor]{cov.shrink}} 
#' (Schäfer & Strimmer, 2005) and the graphical lasso via 
#' \code{glasso::\link[glasso]{glasso}} (Friedman et al., 2007). These latter 
#' methods are preferable in high-dimensional settings where sample covariance 
#' matrices may be unstable or singular. Alternatively, users can pass a 
#' pre-computed covariance matrix via the \code{Sigma} argument.
#' 
#' Uncertainty can be evaluated in leaky IV models using the bootstrap, provided
#' that covariances are estimated internally and not passed directly via the 
#' \code{Sigma} argument. Bootstrapping provides a nonparametric approximate
#' posterior distribution for min/max values of the average treatment effect of 
#' X on Y. Set \code{bayes = TRUE} to replace the classical bootstrap with a 
#' Bayesian bootstrap (Rubin, 1981).
#' 
#' @return  
#' A data frame with columns for \code{ATE_lo} and \code{ATE_hi}, representing
#' lower and upper bounds of the partial identification interval for the 
#' causal effect of \code{x} on \code{y}. When bootstrapping, the output data 
#' frame contains \code{n_boot} rows, one for each bootstrap replicate. 
#' 
#' @references  
#' Schäfer, J., and Strimmer, K. (2005). A shrinkage approach to large-scale 
#' covariance estimation and implications for functional genomics. 
#' \emph{Statist. Appl. Genet. Mol. Biol.}, 4:32.
#' 
#' Friedman, J., Hastie, T., and Tibshirani, R. (2007). Sparse inverse 
#' covariance estimation with the lasso. \emph{Biostatistics}, 9:432-441.
#' 
#' Rubin, D.R. (1981). The Bayesian bootstrap. \emph{Ann. Statist.}, 
#' \emph{9}(1): 130-134. 
#' 
#' @examples  
#' set.seed(123)
#' 
#' # Hyperparameters
#' n <- 200
#' d_z <- 5
#' beta <- rep(1, d_z)
#' gamma <- rep(0.1, d_z)
#' theta <- 2
#' 
#' # Simulate correlated residuals
#' S_eps <- matrix(c(1, 0.5, 0.5, 1), ncol = 2)
#' eps <- matrix(rnorm(n * 2), ncol = 2)
#' eps <- eps %*% chol(S_eps)
#' 
#' # Simulate observables from a leaky IV model
#' z <- matrix(rnorm(n * d_z), ncol = d_z)
#' x <- as.numeric(z %*% beta) + eps[, 1]
#' y <- as.numeric(z %*% gamma) + x * theta + eps[, 2]
#' 
#' # Run the algorithm
#' leakyIV(x, y, z, tau = 1)
#' 
#' @export 
#' @import data.table
#' @importFrom corpcor cov.shrink invcov.shrink
#' @importFrom glasso glasso
#' @importFrom foreach foreach %do% %dopar%

leakyIV <- function(
    x,
    y, 
    z,
    tau, 
    p = 2, 
    n_rho = 1999L,
    method = "mle",
    Sigma = NULL,
    n_boot = NULL, 
    bayes = FALSE, 
    parallel = TRUE, 
    ...) {
  
  # To avoid data.table check issues
  bb <- rho <- sat <- NULL
  
  # Preliminaries
  if (is.matrix(z) || is.data.frame(z)) {
    n_z <- nrow(z)
    d_z <- ncol(z)
  } else {
    n_z <- length(z)
    d_z <- 1L
  }
  if (is.null(n_boot)) {
    n_boot <- 0L
    parallel <- FALSE
  }
  stopifnot(
    is.numeric(z) || is.matrix(z) || is.data.frame(z),
    is.numeric(x), is.numeric(y), is.numeric(tau), is.numeric(p), 
    is.character(method), is.numeric(n_boot), is.logical(bayes), 
    is.logical(parallel)
  )
  if (length(x) != n_z) {
    stop('x and z must have the same number of samples.')
  }
  if (length(y) != n_z) {
    stop('y and z must have the same number of samples.')
  }
  if (!length(tau) %in% c(1, d_z)) {
    stop('tau must be either a scalar or a vector of length ncol(z).')
  }
  if (any(tau < 0)) {
    stop('tau must be strictly positive.')
  }
  if (length(tau) == 1) {
    if (p < 0) {
      stop('p must be >= 0.')
    } else if (p < 1) {
      warning('Exact solutions are only possible for p >= 1, using approximate ',
              'methods instead.')
    }
  }
  if (!method %in% c('mle', 'shrink', 'glasso')) {
    stop('method not recognized. Must be one of "mle", "shrink", or "glasso".')
  }
  if (!is.null(Sigma)) {
    if (ncol(Sigma) != d_z | nrow(Sigma) != d_z) {
      stop('Pre-computed covariance matrix Sigma must have ncol(z) rows and ', 
           'ncol(z) columns.')
    }
  }
  dat <- cbind(z, x, y)
  n <- nrow(dat)
  d <- ncol(dat)
  w <- rep(1, n)
  
  # Bootstrap samples or weights
  if (n_boot > 0) {
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
    
    # Estimate covariance and precision matrices
    if (is.null(Sigma)) {
      if (b > 0) {
        if (isTRUE(bayes)) {
          w <- w_mat[, b]
        } else {
          dat <- dat[b_mat[, b], ]
        }
      }
      if (method == 'mle') {
        Sigma <- stats::cov.wt(dat, wt = w)$cov
        Theta_z <- solve(Sigma[seq_len(d_z), seq_len(d_z)])
      } else if (method == 'shrink') {
        Sigma <- cov.shrink(dat, w = w, verbose = FALSE)[seq_len(d), seq_len(d)]
        Theta_z <- invcov.shrink(dat, verbose = FALSE)[seq_len(d_z), seq_len(d_z)]
      } else if (method == 'glasso') {
        s <- glasso(stats::cov.wt(dat, wt = w)$cov, ...)
        Sigma <- s$w
        Theta_z <- s$wi[seq_len(d_z), seq_len(d_z)]
      }
    } else {
      Theta_z <- solve(Sigma[seq_len(d_z), seq_len(d_z)])
    }
    
    # Extract other covariance parameters
    Sigma_z <- Sigma[seq_len(d_z), seq_len(d_z)]
    Sigma_zy <- matrix(Sigma[seq_len(d_z), d], ncol = 1L)
    Sigma_yz <- t(Sigma_zy)
    Sigma_zx <- matrix(Sigma[seq_len(d_z), d - 1L], ncol = 1L)
    Sigma_xz <- t(Sigma_zx)
    sigma_xy <- Sigma[d, d - 1L]
    var_x <- Sigma[d - 1L, d - 1L]
    var_y <- Sigma[d, d]
    
    # Our magic variables
    eta_x2 <- var_x - as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
    phi2 <- var_y - as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
    psi <- sigma_xy - as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
    
    # Compute theta as a function of rho 
    theta_fn <- function(rho) {
      theta <- (psi / eta_x2) - sign(rho) * 
        ((sqrt((1 - 1/rho^2) * (psi^2 - phi2 * eta_x2)))) / 
        (-eta_x2 * (1 - 1/rho^2))
      return(theta)
    }
    
    if (length(tau) == 1) {
      # Compute norms as a function of rho
      norm_fn <- function(rho) {
        theta <- theta_fn(rho)
        gamma <- as.numeric(Theta_z %*% (Sigma_zy - theta * Sigma_zx))
        norm <- (sum(abs(gamma)^p))^(1 / p)
        return(norm)
      }
      # Compute loss as a function of rho
      loss_fn <- function(rho) {
        norm <- norm_fn(rho)
        loss <- (tau - norm)^2
        return(loss)
      }
      # Find the split point: upper and lower bounds lie on either side of
      # rho_star, assuming tau > min_norm$value
      min_norm <- stats::optim(0, norm_fn, method = 'Brent', lower = -1, upper = 1)
      if (tau < min_norm$value) {
        warning('tau is too low, resulting in an empty feasible region. ',
                'Consider rerunning with a higher threshold.')
        ATE_lo <- ATE_hi <- NA_real_
      } else {
        rho_star <- min_norm$par
        rho_lo <- stats::optim(mean(c(-1, rho_star)), loss_fn, method = 'Brent', 
                               lower = -1, upper = rho_star)$par
        rho_hi <- stats::optim(mean(c(rho_star, 1)), loss_fn, method = 'Brent', 
                               lower = rho_star, upper = 1)$par
        ATE_lo <- theta_fn(rho_hi)
        ATE_hi <- theta_fn(rho_lo)
      }
    } else {
      # Check criterion
      tmp <- data.table(rho = seq(-0.999, 0.999, length.out = n_rho))
      tmp <- tmp[rho != 0]
      tmp[, sat := sapply(rho, function(r) {
        theta <- theta_fn(r)
        gamma <- as.numeric(Theta_z %*% (Sigma_zy - theta * Sigma_zx))
        out <- all(abs(gamma) <= tau)
        return(out)
      })]
      if (tmp[, sum(sat)] == 0) {
        warning('tau is too restrictive, resulting in an empty feasible region. ',
                'Consider rerunning with different thresholds.')
        ATE_lo <- ATE_hi <- NA_real_
      } else if (tmp[, sum(sat)] == 1) {
        ATE_lo <- ATE_hi <- tmp[sat == TRUE, rho]
      } else {
        rho_lo <- tmp[sat == TRUE, min(rho)]
        rho_hi <- tmp[sat == TRUE, max(rho)]
        ATE_lo <- theta_fn(rho_hi)
        ATE_hi <- theta_fn(rho_lo)
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


