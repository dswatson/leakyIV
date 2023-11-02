#' Bounding Causal Effects with Leaky Instruments
#' 
#' This function estimates bounds on average treatment effects in linear IV 
#' models under limited violations of the exclusion criterion.
#' 
#' @param z One or more instrumental variable(s). 
#' @param x Treatment variable. 
#' @param y Outcome variable.
#' @param tau Upper bound on the p-norm of linear weights \code{gamma} from 
#'   \code{z} to \code{y}. See details.
#' @param p Power of the norm on \code{gamma}.
#' @param tau_fn Alternatively, an indicator function on \code{gamma} that 
#'   evaluates to \code{TRUE} if the information leakage from \code{z} to
#'   \code{y} is not too severe, else \code{FALSE}. See examples.
#' @param rho_min Minimum correlation between residuals for \code{x} and \code{y}.
#' @param rho_max Maximum correlation between residuals for \code{x} and \code{y}.
#' @param n_rho Number of candidate rhos to evaluate. By default, these will be 
#'   evenly spaced across \code{[rho_min, rho_max]}. Alternatively, a sequence 
#'   of values for \code{rho}.  
#' @param n_boot Number of bootstrap replicates for standard error estimates.
#' @param bayes Use Bayesian bootstrap?
#' @param parallel Compute bootstrap estimates in parallel? Must register 
#'   backend beforehand, e.g. via \code{doMC}. 
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
#' \code{leaky_iv} provides native support for thresholding the p-norm of
#' linear weights \code{gamma} in the structural equation for \code{y}:
#' \eqn{Y := Z \gamma + X \theta + \epsilon_Y}. Both the threshold \code{tau}
#' and the power of the norm \code{p} are up to the user. Alternatively, users
#' may supply any Boolean condition on \code{gamma} via the \code{tau_fn}
#' argument (see Examples).
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
#' df <- leaky_iv(dat$z, dat$x, dat$y, tau = 1)
#' 
#' # Compute 95% confidence interval
#' alpha <- 0.05
#' ci_hi <- df$ate + df$se * qnorm(1 - alpha/2)
#' ci_lo <- df$ate - df$se * qnorm(1 - alpha/2)
#' 
#' # Compute 95% credible interval
#' df <- leaky_iv(dat$z, dat$x, dat$y, tau = 1, bayes = TRUE)
#' ci_hi <- df$ate + df$se * qnorm(1 - alpha/2)
#' ci_lo <- df$ate - df$se * qnorm(1 - alpha/2)
#' 
#' # Supply your own tau_fn
#' 
#' # Restrict the range of rho
#' r <- seq(0.01, 0.99, length.out = 100)
#' df <- leaky_iv(dat$z, dat$x, dat$y, tau = 1, n_rho = r)
#' 
#' @export 
#' @import data.table 
#' @importFrom matrixStats cov.wt
#' @importFrom foreach foreach
#' @importFrom doRNG dorng

leaky_iv <- function(
    z,
    x,
    y,
    tau, 
    p = 2L,
    tau_fn = NULL,
    rho_min = -0.99,
    rho_max = 0.99,
    n_rho = 199L, 
    n_boot = 199L,
    bayes = FALSE, 
    parallel = TRUE
  ) {
  
  # Prelimz
  if (is.matrix(z) || is.data.frame(z)) {
    n_z <- nrow(z)
    d_z <- ncol(z)
  } else {
    n_z <- length(z)
    d_z <- 1L
  }
  if (!is.null(tau_fn)) {
    stopifnot(is.logical(tau_fn(rnorm(d_z))))
  }
  stopifnot(
    is.numeric(z) || is.matrix(z) || is.data.frame(z),
    is.numeric(x), is.numeric(y), !is.null(tau), is.numeric(p),
    is.numeric(n_rho), is.numeric(n_boot), is.logical(bayes), 
    is.logical(parallel),
    length(x) == n_z,
    length(y) == n_z
  )
  dat <- cbind(z, x, y)
  n <- nrow(dat)
  d <- ncol(dat)
  if (length(n_rho) == 1L) {
    rhos <- seq(rho_min, rho_max, length.out = n_rho)
  } else {
    rhos <- n_rho
  }
  if (any(rhos == 0)) {
    zero_idx <- which(rhos == 0)
    rhos <- rhos[-zero_idx]
  }
  
  # Define bootstrap loop
  boot_loop <- function(b) {
    if (n_boot == 1L) {
      Sigma <- cov(dat)
    } else {
      if (bayes) {
        # Draw Dirichlet weights
        wts <- rexp(n)
        wts <- (wts / sum(wts)) * n
        # Estimate eta_x
        f1 <- lm(x ~ ., data = dat[, 1:(d - 1)], weights = wts)
        eta_x <- sqrt(weighted.mean(x = residuals(f1)^2, w = wts))
        # Estimate data covariance
        Sigma <- cov.wt(dat, wt = wts)$cov
      } else {
        # Draw bootstrap sample
        tmp <- dat[sample.int(n, replace = TRUE)]
        # Estimate eta_x
        f1 <- lm(x ~ ., data = tmp[, 1:(d - 1)])
        eta_x <- sqrt(mean((residuals(f1)^2)))
        # Estimate data covariance
        Sigma <- cov(tmp)
      }
    }
    # Extract elements of covariance matrix
    Sigma_z <- Sigma[seq_len(d_z), seq_len(d_z)]
    Theta_z <- solve(Sigma_z)
    Theta_z2 <- Theta_z %*% Theta_z
    Sigma_zy <- matrix(Sigma[seq_len(d_z), d], ncol = 1L)
    Sigma_yz <- t(Sigma_zy)
    Sigma_zx <- matrix(Sigma[seq_len(d_z), d - 1L], ncol = 1L)
    Sigma_xz <- t(Sigma_zx)
    sigma_xy <- Sigma[d, d - 1L]
    var_x <- Sigma[d - 1L, d - 1L]
    var_y <- Sigma[d, d]
    # Define
    dd <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zy)
    ee <- as.numeric(Sigma_xz %*% Theta_z %*% Sigma_zx)
    ff <- as.numeric(Sigma_yz %*% Theta_z %*% Sigma_zy)
    if (!is.null(tau_fn) & p == 2L) {
      if (tau < ff - ee^2 / dd) {
        stop('tau threshold too low, select larger value.')
      }
      # Compute feasible region according to rho
      theta_from_rho <- function(rho_i) {
        gg <- (ee - var_x) * (1 + ((ee - var_x) / (eta_x^2 * rho_i^2)))
        hh <- -(dd - sigma_xy) * (1 + ((ee - var_x) / (eta_x^2 * rho_i^2)))
        ii <- ((dd - sigma_xy)^2 / (eta_x^2 * rho_i^2)) + ff - var_y
        if (rho_i < 0) {
          theta_r <- as.numeric((-hh + sqrt(hh^2 - gg * ii)) / gg)
        } else {
          theta_r <- as.numeric((-hh - sqrt(hh^2 - gg * ii)) / gg)
        }
        return(theta_r)
      }
      theta_min <- theta_from_rho(rho_max)
      theta_max <- theta_from_rho(rho_min)
      # Define
      aa <- as.numeric(Sigma_xz %*% Theta_z2 %*% Sigma_zx)
      bb <- as.numeric(Sigma_xz %*% Theta_z2 %*% Sigma_zy)
      cc <- as.numeric(Sigma_yz %*% Theta_z2 %*% Sigma_zy)
      out <- data.table(
        'bound' = c('lo', 'hi'), 
        'theta' = c((bb - sqrt(aa * (tau - cc) + bb^2)) / aa,
                    (bb + sqrt(aa * (tau - cc) + bb^2)) / aa)
      )
      if (out[bound == 'lo', theta] < theta_min | 
          out[bound == 'hi', theta] > theta_max) {
        warning('tau-feasible bounds exceed rho-feasible bounds.
                This is only possible under *extreme* confounding. 
                Consider reducing tau.')
      }
    #} else if (!is.null(tau_fn) & p == 1L) {
      # BLAH
    } else {
      # Compute theta and gamma as functions of rho
      rho_loop <- function(rho_i) {
        gg <- (ee - var_x) * (1 + ((ee - var_x) / (eta_x^2 * rho_i^2)))
        hh <- -(dd - sigma_xy) * (1 + ((ee - var_x) / (eta_x^2 * rho_i^2)))
        ii <- ((dd - sigma_xy)^2 / (eta_x^2 * rho_i^2)) + ff - var_y
        if (rho_i < 0) {
          theta_r <- as.numeric((-hh + sqrt(hh^2 - gg * ii)) / gg)
        } else {
          theta_r <- as.numeric((-hh - sqrt(hh^2 - gg * ii)) / gg)
        }
        gamma_r <- as.numeric(Theta_z %*% (Sigma_zy - theta_r * Sigma_zx))
        if (is.null(tau_fn)) {
          norm <- (sum(abs(gamma_r)^p))^(1 / p)
          sat <- NA
        } else {
          norm <- NA_real_
          sat <- tau_fn(gamma_r)
        }
        # Export
        out <- data.table(
          'rho' = rho_i, 'theta' = theta_r, 'norm' = norm, 'sat' = sat
        )
        return(out)
      }
      rho_grid <- foreach(r = rhos, .combine = rbind) %do% rho_loop(r) 
      if (is.null(tau_fn)) {
        rho_grid[, sat := ifelse(norm <= tau, TRUE, FALSE)]
      }
      rho_grid <- na.omit(rho_grid[sat == TRUE])
      if (nrow(rho_grid) == 0L) {
        out <- NULL
      } else {
        out <- data.table(
          'bound' = c('lo', 'hi'), 
          'theta' = c(rho_grid[, min(theta)], rho_grid[, max(theta)])
        )
      }
    }
    return(out)
  }
  
  # Run bootstrap
  if (isTRUE(parallel)) {
    boots <- foreach(i = seq_len(n_boot), .combine = rbind) %dorng% 
      boot_loop(i)
  } else {
    boots <- foreach(i = seq_len(n_boot), .combine = rbind) %do% 
      boot_loop(i)
  }
  out <- data.table(
    'bound' = c('lo', 'hi'),
    'ate' = c(boots[bound == 'lo', mean(theta)], 
              boots[bound == 'hi', mean(theta)]),
    'se' = c(boots[bound == 'lo', sd(theta)],
             boots[bound == 'hi', sd(theta)])
  )
  return(out)
}



