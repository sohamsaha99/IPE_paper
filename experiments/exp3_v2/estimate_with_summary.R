library(dplyr)
library(tidyr)
library(purrr)
library(mvtnorm)
source("./experiments/exp3_v2/weighted_rmst.R")
source("./experiments/exp3_v2/impute_censored.R")

#' Function for estimating theta_min, theta_max
#' Input df need not be normalized, but the
#' predictor columns are named age, sex, kps, mgmt, isGTR and
#' everything numeric.
#' Returns a list with alpha_min, theta_min, wts_min,
#' alpha_max, theta_max, wts_max,
#' grid_points, L, tau
estimate_IPE <- function(df,
                         target_con,
                         control_rmst_value,
                         tau,
                         grid_points,
                         L) {
  # Normalize columns and adjust target accordingly
  X_cols <- c("age", "sex", "kps", "mgmt", "isGTR")
  # Store the original moments
  original_moments <- target_con
  for (x in X_cols) {
    m <- mean(df[, x])
    s <- sd(df[, x])
    df[, x] <- (df[, x] - m) / s
    # t1, t2 might be empty, but not a problem, checked later
    t1 <- original_moments %>% filter(xname == x, type == 1) %>% pull(value)
    t2 <- original_moments %>% filter(xname == x, type == 2) %>% pull(value)
    # Adjust target
    for (i in seq_len(nrow(target_con))) {
      if (target_con[i, "xname"] != x) next
      if (target_con[i, "type"] == 1) {
        # only update if an original t1 exists (and is unique)
        if (length(t1) == 1) {
          target_con[i, "value"] <- (original_moments[i, "value"] - m) / s
        } else {
          stop("t1 is problematic. Check.")
        }
      } else if (target_con[i, "type"] == 2) {
        # need t1 to transform E[X^2] -> E[Z^2]
        if (length(t2) == 1 && length(t1) == 1) {
          target_con[i, "value"] <- (original_moments[i, "value"] - 2 * m * t1 + m * m) / (s ^ 2)
        } else {
          stop("t1 or t2 problematic. Check.")
        }
      }
    }
  }
  # Now we have target_con and df normalized
  # Impute censored patients, stored in column emin_tau
  df_trt_imputed <- impute_censored(df = df, method = "emin_tau", tau = tau)
  # Compute dmvnorm
  normal_den <- function(x, v) {
    stopifnot(ncol(x) == length(v))
    dmvnorm(x, mean = v, sigma = diag(length(v)) / (L^2))
  }
  # Y column contains E(min{T, tau}|X, Y, delta)
  df_trt_imputed$Y <- df_trt_imputed$emin_tau
  # If "prob" does not exist, create equal probabilities
  if (!"prob" %in% colnames(df_trt_imputed)) {
    df_trt_imputed$prob <- rep(1 / nrow(df_trt_imputed), nrow(df_trt_imputed))
  }
  # Ensure probabilities sum to 1
  df_trt_imputed$prob <- df_trt_imputed$prob / sum(df_trt_imputed$prob)
  df_trt <- df_trt_imputed
  # Set linear programming:
  # min c^\top alpha
  # subject to A alpha = target_con$value
  K <- nrow(grid_points)
  cvec <- rep(NA, K)
  for (k in seq_len(K)) {
    cvec[k] <- sum(df_trt$Y * df_trt$prob *
                   normal_den(df_trt[, X_cols], as.numeric(grid_points[k, ])))
  }
  rhsvec <- c(1, target_con$value)
  dirvec <- rep("==", length(rhsvec))
  Amat <- matrix(NA, nrow = length(rhsvec), ncol = K)
  # Set up each row
  for (k in seq_len(K)) {
    den_k <- normal_den(df_trt[, X_cols], as.numeric(grid_points[k, ]))
    stopifnot(length(den_k) == nrow(df_trt))
    Amat[1, k] <- sum(df_trt$prob  * den_k)
    Amat[2, k] <- sum(df_trt$age   * df_trt$prob * den_k)
    Amat[3, k] <- sum(df_trt$age^2 * df_trt$prob * den_k)
    Amat[4, k] <- sum(df_trt$sex   * df_trt$prob * den_k)
    Amat[5, k] <- sum(df_trt$kps   * df_trt$prob * den_k)
    Amat[6, k] <- sum(df_trt$mgmt  * df_trt$prob * den_k)
    Amat[7, k] <- sum(df_trt$isGTR * df_trt$prob * den_k)
  }

  # Solve linear program
  # Solve LP: min c^T alpha  s.t. Amat %*% alpha = rhsvec, alpha >= 0
  lp_out <- lpSolve::lp(direction    = "min",
                        objective.in = cvec,
                        const.mat    = Amat,
                        const.dir    = dirvec,
                        const.rhs    = rhsvec
  )

  if (lp_out$status != 0) {
    stop("LP did not solve to optimality. lpSolve status code: ", lp_out$status)
  }

  alpha_min <- lp_out$solution
  EY1_min <- as.numeric(lp_out$objval)
  theta_min <- EY1_min - control_rmst_value
  wts_min <- rep(NA, nrow(df_trt))
  for (i in seq_len(nrow(df_trt))) {
    wts_min[i] <- density_ratio_wf(x = as.numeric(df_trt[i, X_cols]),
                                   alpha = alpha_min,
                                   grid_points = grid_points,
                                   L = L
    )
  }

  # Solve linear program
  # Solve LP: max c^T alpha  s.t. Amat %*% alpha = rhsvec, alpha >= 0
  lp_out <- lpSolve::lp(direction    = "max",
                        objective.in = cvec,
                        const.mat    = Amat,
                        const.dir    = dirvec,
                        const.rhs    = rhsvec
  )

  if (lp_out$status != 0) {
    stop("LP did not solve to optimality. lpSolve status code: ", lp_out$status)
  }

  alpha_max <- lp_out$solution
  EY1_max <- as.numeric(lp_out$objval)
  theta_max <- EY1_max - control_rmst_value
  wts_max <- rep(NA, nrow(df_trt))
  for (i in seq_len(nrow(df_trt))) {
    wts_max[i] <- density_ratio_wf(x = as.numeric(df_trt[i, X_cols]),
                                   alpha = alpha_max,
                                   grid_points = grid_points,
                                   L = L
    )
  }

  list(
       theta_min   = theta_min,
       alpha_min   = alpha_hat,
       wts_min     = wts_min,
       theta_max   = theta_max,
       alpha_max   = alpha_max,
       wts_max     = wts_max,
       prob        = df_trt$prob,
       grid_points = grid_points,
       L           = L,
       tau         = tau
  )
}

# Visualize the optimizer weights
# w_f(x) = sum_k alpha_k * N(x | v_k, I_p / L^2)
density_ratio_wf <- function(x, alpha, grid_points, L) {
  # x: numeric vector in R^p
  # alpha: numeric vector length K
  # grid_points: data.frame / matrix with K rows, p cols (v_1,...,v_K)
  # L: positive scalar

  x <- as.numeric(x)
  V <- as.matrix(grid_points)
  alpha <- as.numeric(alpha)

  p <- length(x)
  K <- nrow(V)

  stopifnot(ncol(V) == p)
  stopifnot(length(alpha) == K)
  stopifnot(is.finite(L), L > 0)

  # Evaluate phi(x; v_k, I_p / L^2) for all k
  Sigma <- diag(p) / (L^2)
  dens <- mvtnorm::dmvnorm(V, mean = x, sigma = Sigma)  # length K

  sum(alpha * dens)
}

