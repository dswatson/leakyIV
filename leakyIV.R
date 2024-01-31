#' Bounding Causal Effects with Leaky Instruments
#' 
#' This function estimates bounds on average treatment effects in linear IV 
#' models under limited violations of the exclusion criterion.
#' 
#' @param x Treatment variable. 
#' @param y Outcome variable.
#' @param z One or more candidate instrumental variable(s). 
#' @param tau Upper bound on the p-norm of linear weights \code{gamma} from 
#'   \code{z} to \code{y}. See details.
#' @param p Power of the norm on \code{gamma}.
#' @param beta Method for regressing \code{x} on \code{z}. Options include 
#'   \code{"ols"} (default), \code{"ridge"}, and \code{"lasso"}. For regularized
#'   methods, the value of the Lagrange multiplier is learned via 10-fold CV. 
#'   Alternatively, a numeric vector of length \code{ncol(z)} specifying linear 
#'   coefficients. 
#' @param n_boot Number of bootstrap replicates.
#' @param bayes Use Bayesian bootstrap?
#' @param parallel Compute bootstrap estimates in parallel?
#' 
#' @details 
#' Instrumental variables are defined by three structural assumptions: they must
#' be (A1) \emph{relevant}, i.e. associated with the treatment; (A2) 
#' \emph{unconfounded}, i.e. independent of common causes between treatment and 
#' outcome; and (A3) \emph{exclusive}, i.e. only affect outcomes through the 
#' treatment. The \code{leaky_iv} algorithm relaxes (A3), allowing some 
#' information leakage from IVs \code{z} to outcomes \code{y} in linear 
#' structural equation models (SEMs). While causal effects are no longer 
#' identifiable in this setting, tight bounds can be computed with precision.
#' 
#' Violations of the exclusion restriction may come in many forms. 
#' \code{leakyIV} provides native support for thresholding the p-norm of
#' linear weights \code{gamma} in the structural equation for \code{y}:
#' \eqn{Y := Z \gamma + X \theta + \epsilon_Y}. Both the threshold \code{tau}
#' and the power of the norm \code{p} are up to the user. 
#' 
#' The uncertainty of associated bounds is estimated via the bootstrap, with 
#' optional support for Bayesian bootstrapping (Rubin, 1981).
#' 
#' @return  
#' A 2x3 data.frame with the following columns: \code{bound} (lower or upper),
#' \code{ate} (average treatment effect), and \code{se} (standard error). 
#' 
#' @references  
#' Rubin, D.R. (1981). The Bayesian bootstrap. \emph{Ann. Statist.}, 
#' \emph{9}(1): 130-134. 
#' 
#' @examples  
#' # Load data
#' BLAH
#' 
#' # Run leaky IV
#' df <- leakyIV(dat$z, dat$x, dat$y, tau = 1)
#' 
#' # Compute 95% confidence interval
#' 
#' @export 
#' @import data.table 
#' @import glmnet
#' @importFrom matrixStats cov.wt
#' @importFrom foreach foreach

leakyIV <- function(
    x,
    y, 
    z,
    tau, 
    p = 2, 
    beta = 'ols',
    n_boot = 1999, 
    bayes = FALSE, 
    parallel = TRUE) {
  
  # Preliminaries
  if (is.matrix(z) || is.data.frame(z)) {
    n_z <- nrow(z)
    d_z <- ncol(z)
  } else {
    n_z <- length(z)
    d_z <- 1L
  }
  stopifnot(
    is.numeric(z) || is.matrix(z) || is.data.frame(z),
    is.numeric(x), is.numeric(y), is.numeric(p), # SOMETHING FOR BETA?
    is.numeric(n_boot), is.logical(bayes), 
    is.logical(parallel),
    length(x) == n_z,
    length(y) == n_z
  )
  dat <- cbind(z, x, y)
  n <- nrow(dat)
  d <- ncol(dat)
  d_z <- d - 2L
  if (beta == 'ridge') {
    alpha_glmnet <- 1
  } else if (beta == 'lasso') {
    alpha_glmnet <- 0
  }
  
  # Bootstrap samples or weights
  if (n_boot > 0) {
    if (isTRUE(bayes)) {
      # Draw Dirichlet weights
      w_mat <- matrix(rexp(n * n_boot), ncol = n_boot)
      w_mat <- (w_mat / rowSums(w_mat)) * n
    } else {
      # Draw bootstrap samples
      b_mat <- matrix(sample(n, n * n_boot, replace = TRUE), ncol = n_boot)
    }
  }

  # Compute bounds
  loop <- function(b) {
    
    # Compute Sigma and eta_x
    if (n_boot > 1 & isTRUE(bayes)) {
      wts <- w_mat[, b]
      # Estimate eta_x
      if (beta == 'ols') {
        f1 <- lm(x ~ ., data = dat[, 1:(d - 1)], weights = wts)
        eps_x <- residuals(f1)
      } else if (beta %in% c('ridge', 'lasso')) {
        f1 <- cv.glmnet(z, x, alpha = alpha_glmnet, weights = wts)
        eps_x <- x - predict(f1, z, s = 'lambda.min')
      } else {
        eps_x <- x - beta %*% z # WITH WEIGHTS THO?
      }
      eta_x <- sqrt(weighted.mean(x = eps_x^2, w = wts))
      # Estimate data covariance
      Sigma <- cov.wt(dat, wt = wts)$cov
    } else {
      if (n_boot == 0) {
        tmp <- dat
      } else {
        # Draw bootstrap sample
        tmp <- dat[b_mat[, b], ]
      }
      # Estimate eta_x
      if (beta == 'ols') {
        f1 <- lm(x ~ ., data = tmp[, 1:(d - 1)])
        eps_x <- residuals(f1)
      } else if (beta %in% c('ridge', 'lasso')) {
        f1 <- cv.glmnet(tmp[, seq_len(d_z)], tmp$x, alpha = alpha_glmnet)
        eps_x <- tmp$x - predict(f1, tmp[, seq_len(d_z)], s = 'lambda.min')
      } else {
        eps_x <- x - beta %*% z
      }
      eta_x <- sqrt(mean(eps_x^2))
      # Estimate data covariance
      Sigma <- cov(tmp)
    }
    
    # Extract elements of covariance matrix
    Sigma_z <- Sigma[seq_len(d_z), seq_len(d_z)]
    Theta_z <- solve(Sigma_z)
    Sigma_zy <- matrix(Sigma[seq_len(d_z), d], ncol = 1L)
    Sigma_yz <- t(Sigma_zy)
    Sigma_zx <- matrix(Sigma[seq_len(d_z), d - 1L], ncol = 1L)
    Sigma_xz <- t(Sigma_zx)
    sigma_xy <- Sigma[d, d - 1L]
    var_x <- Sigma[d - 1L, d - 1L]
    var_y <- Sigma[d, d]
    
    # Bounding
    aa <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx) - var_x
    bb <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy) - sigma_xy
    cc <- as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy) - var_y
    theta_fn <- function(rho) {
      theta <- bb / aa - sign(rho) * sqrt((bb^2 - aa * cc) * (1 + aa / (eta_x^2 * rho^2))) / 
        (aa * (1 + aa / (eta_x^2 * rho^2)))
      return(theta)
    }
    norm_fn <- function(rho) {
      theta <- theta_fn(rho)
      gamma <- as.numeric(Theta_z %*% (Sigma_zy - theta * Sigma_zx))
      norm <- (sum(abs(gamma)^p))^(1 / p)
      return(norm)
    }
    loss_fn <- function(rho) {
      norm <- norm_fn(rho)
      loss <- (tau - norm)^2
      return(loss)
    }
    min_norm <- optim(0, norm_fn, method = 'Brent', lower = -1, upper = 1)
    if (tau < min_norm$value) {
      warning('tau is too low, resulting in an empty feasible region. ',
              'Consider rerunning with a higher threshold.')
      theta_lo <- theta_hi <- NA_real_
    } else {
      rho_star <- min_norm$par
      rho_lo <- optim(mean(c(-1, rho_star)), loss_fn, method = 'Brent', 
                      lower = -1, upper = rho_star)$par
      rho_hi <- optim(mean(c(rho_star, 1)), loss_fn, method = 'Brent', 
                      lower = rho_star, upper = 1)$par
      theta_lo <- theta_fn(rho_hi)
      theta_hi <- theta_fn(rho_lo)
    }
    
    # Export
    out <- data.table(b, theta_lo, theta_hi)
    return(out)
  }
  
  # Optionally compute in parallel
  if (n_boot == 0) {
    out <- loop(0)
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


