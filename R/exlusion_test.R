#' Testing Exclusion
#' 
#' Performs a Monte Carlo test against the null hypothesis that minimum leakage
#' is zero, a necessary but insufficient condition for exclusion.
#' 
#' @param dat Input data. Either (a) an \eqn{n \times d} data frame or matrix of 
#'   observations with columns for treatment, outcome, and candidate instruments; 
#'   or (b) a \eqn{d \times d} covariance matrix over such variables. Note that 
#'   in either case, the order of variables is presumed to be treatment 
#'   (\eqn{X}), outcome (\eqn{Y}), leaky instruments (\eqn{Z}). 
#'   \code{exclusion_test} requires at least two candidate instruments \eqn{Z}.
#' @param normalize Scale candidate instruments to unit variance?
#' @param method Estimator for the covariance matrix. Options include 
#'   (a) \code{"mle"}, the default; (b) \code{"shrink"}, an analytic empirical 
#'   Bayes solution; or (c) \code{"glasso"}, the graphical lasso. See details.
#' @param approx Use nearest positive definite approximation if the estimated 
#'   covariance matrix is singular? See details.
#' @param n_sim Number of Monte Carlo replicates.
#' @param parallel Run Monte Carlo simulations in parallel? Must register 
#'   backend beforehand, e.g. via \code{doParallel}.
#' @param return_stats Return observed statistic and simulated null distribution?
#' @param ... Extra arguments to be passed to graphical lasso estimator if
#'   \code{method = "glasso"}. Note that the regularization parameter \code{rho}
#'   is required as input, with no default. 
#' 
#' 
#' @details 
#' The classic linear instrumental variable (IV) model relies on the 
#' \emph{exclusion criterion}, which states that instruments \eqn{Z} have no 
#' direct effect on the outcome \eqn{Y}, but can only influence it through the
#' treatment \eqn{X}. This implies a series of tetrad constraints that can be
#' directly tested, given a model for sampling data from the covariance matrix 
#' of the observable variables (Watson et al., 2024). 
#' 
#' We assume that data are multivariate normal and impose the null hypothesis
#' by modifying the estimated covariance matrix to induce a linear dependence
#' between the vectors for Cov(\eqn{Z, X}) and Cov(\eqn{Z, Y}). Our test 
#' statistic is the determinant of the cross product of these vectors, which 
#' equals zero if and only if the null hypothesis is true. We generate a null
#' distribution by simulating from the null covariance matrix and compute a
#' \emph{p}-value by estimating the proportion of statistics that exceed the 
#' observed value. Future releases will provide support for a wider range of 
#' data generating processes.
#' 
#' Numerous methods exist for estimating covariance matrices. \code{exclusion_test}
#' provides support for maximum likelihood estimation (the default), as well as
#' empirical Bayes shrinkage via \code{corpcor::\link[corpcor]{cov.shrink}} 
#' (Schäfer & Strimmer, 2005) and the graphical lasso via 
#' \code{glasso::\link[glasso]{glasso}} (Friedman et al., 2007). These latter 
#' methods are preferable in high-dimensional settings where sample covariance 
#' matrices may be unstable or singular. Alternatively, users can pass a 
#' pre-computed covariance matrix directly as \code{dat}.
#' 
#' Estimated covariance matrices may be singular for some datasets or Monte 
#' Carlo samples. Behavior in this case is determined by the \code{approx} 
#' argument. If \code{TRUE}, the test proceeds with the nearest positive 
#' definite approximation, computed via Higham's (2002) algorithm (with a 
#' warning). If \code{FALSE}, the sampler will attempt to use the singular 
#' covariance matrix (also with a warning), but results may be invalid.
#' 
#' 
#' @return  
#' Either a scalar representing the Monte Carlo \emph{p}-value of the exclusion 
#' test (default) or, if \code{return_stats = TRUE}, a named list with three 
#' entries: \code{psi}, the observed statistic; \code{psi0}, a vector of length 
#' \code{n_sim} with simulated null statistics; and \code{p_value}, the 
#' resulting \emph{p}-value.
#' 
#' 
#' @references  
#' Watson, D., Penn, J., Gunderson, L., Bravo-Hermsdorff, G., Mastouri, A., and
#' Silva, R. (2024). Bounding causal effects with leaky instruments. \emph{arXiv}
#' preprint, 2404.04446.
#' 
#' Spirtes, P. Calculation of entailed rank constraints in partially 
#' non-linear and cyclic models. In \emph{Proceedings of the 29th Conference on 
#' Uncertainty in Artificial Intelligence}, 606–615, 2013.
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
#' 
#' @examples  
#' set.seed(123)
#' 
#' # Hyperparameters
#' n <- 200
#' d_z <- 4
#' beta <- rep(1, d_z)
#' theta <- 2
#' rho <- 0.5
#' 
#' # Simulate correlated residuals
#' S_eps <- matrix(c(1, rho, rho, 1), ncol = 2)
#' eps <- matrix(rnorm(n * 2), ncol = 2)
#' eps <- eps %*% chol(S_eps)
#' 
#' # Simulate observables from the linear IV model
#' z <- matrix(rnorm(n * d_z), ncol = d_z)
#' x <- z %*% beta + eps[, 1]
#' y <- x * theta + eps[, 2]
#' obs <- cbind(x, y, z)
#' 
#' # Compute p-value of the test
#' exclusion_test(obs, parallel = FALSE)
#' 
#' @export 
#' @importFrom corpcor is.positive.definite cov.shrink 
#' @importFrom glasso glasso
#' @importFrom stats cov
#' @importFrom Matrix nearPD
#' @importFrom mvnfast rmvn
#' @importFrom foreach foreach %dopar%

exclusion_test <- function(
    dat,
    normalize = TRUE,
    method = "mle", 
    approx = TRUE,
    n_sim = 1999L, 
    parallel = TRUE, 
    return_stats = FALSE,
    ...) {
  
  # To avoid data.table check issues
  bb <- NULL
  
  # Checks, warnings
  if (nrow(dat) == ncol(dat)) {
    Sigma <- dat
    if (!is.positive.definite(Sigma)) {
      stop('Pre-computed covariance matrix must be positive definite.')
    }
    Sigma_input <- TRUE
  } else {
    dat <- as.data.frame(dat)
    n <- nrow(dat)
    Sigma_input <- FALSE
  }
  if (!isTRUE(Sigma_input) & !method %in% c('mle', 'shrink', 'glasso')) {
    stop('method not recognized. Must be one of "mle", "shrink", or "glasso".')
  }
  d <- ncol(dat)
  if (d <= 3L) {
    stop('exclusion_test requires at least two candidate instruments.')
  }
  
  # Compute covariance
  cov_fn <- function(input_data) {
    if (method == 'mle') {
      Sigma <- stats::cov(input_data)
    } else if (method == 'shrink') {
      Sigma <- cov.shrink(input_data, verbose = FALSE)[seq_len(d), seq_len(d)]
    } else if (method == 'glasso') {
      Sigma <- glasso(stats::cov(input_data)$cov, ...)$w
    }
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
    if (!is.positive.definite(Sigma)) {
      if (isTRUE(approx)) {
        warning('Covariance estimator results in a singular matrix. ',
                'Substituting nearest positive definite approximation.')
        Sigma <- as.matrix(nearPD(Sigma)$mat)
      } else {
        warning('Covariance estimator results in a singular matrix. ',
                'Consider rerunning with another method or setting approx = TRUE.')
      }
    }
    return(Sigma)
  }
  if (!isTRUE(Sigma_input)) {
    Sigma <- cov_fn(dat)
  }
  
  # Compute test statistic
  psi_fn <- function(Sigma) {
    Lambda <- Sigma[3:d, 1:2]
    psi <- det(crossprod(Lambda))
    return(psi)
  }
  psi <- psi_fn(Sigma)
  
  # Impose H0
  Sigma0 <- Sigma
  Sigma_zx <- Sigma[3:d, 1L]
  Sigma_zy <- Sigma[3:d, 2L]
  theta_2sls <- as.numeric((Sigma_zx %*% Sigma_zy) / (Sigma_zx %*% Sigma_zx))
  Sigma0[3:d, 2L] <- Sigma0[2L, 3:d] <- theta_2sls * Sigma_zx 
  
  # Draw Monte Carlo samples, compute null statistics and p-value
  null <- rmvn(n * n_sim, mu = rep(0, d), Sigma0)
  mc_loop <- function(b) {
    tmp <- null[((b - 1) * n + 1):(b * n), ]
    S_tmp <- cov_fn(tmp)
    return(psi_fn(S_tmp))
  }
  if (isTRUE(parallel)) {
    psi0 <- foreach(bb = seq_len(n_sim), .combine = c) %dopar% mc_loop(bb)
  } else {
    psi0 <- sapply(seq_len(n_sim), mc_loop)
  }
  p_value <- (sum(psi0 >= psi) + 1L) / (n_sim + 1L)
  
  if (isTRUE(return_stats)) {
    out <- list('psi' = psi, 'psi0' = psi0, 'p_value' = p_value)
  } else {
    out <- p_value
  }
  return(out)
  
}



