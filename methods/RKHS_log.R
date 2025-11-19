source("../utils/tabulate.R", chdir = TRUE)
source("../utils/my_ROI_solver.R", chdir = TRUE)
library(ROI)
library(ROI.plugin.nloptr)
library(nloptr)

# Function to compute theta_min and wts_min based on categorized data
# Input contains:
# Experimental: for each level of X, the proportions and mean(Y) within the level
# Control: Summary statistics of X and mean(Y)
# Centers of Gaussian kernels: grid_points,
# Regularization: L, the upper-bound on RKHS norm of log w
# Regularization: sigm, the standard deviation used in Gaussian kernel
find_theta_min_table_rkhs_log <- function(categorized, control_moments, control_mean, L, sigm) {
  # ---- Prep & ordering ----
  sample_space <- categorized$X
  K <- length(sample_space)

  ord <- seq_along(sample_space)
  if (!is.ordered(sample_space)) {
    ord <- order(sample_space)
    sample_space <- sample_space[ord]
    categorized <- categorized[ord, ]
  }

  # ---- Linear equality constraints: Amat %*% w == rhsvec ----
  Amat <- rbind(categorized$prop,
                categorized$prop * sample_space,
                categorized$prop * sample_space^2
  )
  rhsvec <- c(1, control_moments)  # length 3

  # ---- Kernel + inverse (with small ridge) ----
  Sigma <- outer(sample_space, sample_space, function(a, b) exp(-0.5 * ((a - b) / sigm)^2))
  lambda <- 1e-8
  Kinv <- solve(Sigma + diag(lambda, K))

  # ---- Reparam: u = log w - 1  ==>  w = exp(u + 1)
  # This ensures w > 0 and makes the inequality constraint simply u' Kinv u <= L^2.
  to_w <- function(u) exp(u + 1.0)

  # ---- Objective: minimize sum(cvec * w) with w = exp(u+1) ----
  cvec <- categorized$prop * categorized$mean_Y

  eval_f <- function(u) {
    w <- to_w(u)
    sum(cvec * w)
  }
  eval_grad_f <- function(u) {
    w <- to_w(u)
    # d/du_j [ sum_i c_i w_i ] = c_j * w_j
    cvec * w
  }

  # ---- Equality constraints: g_eq(u) = Amat %*% w - rhsvec = 0 ----
  eval_g_eq <- function(u) {
    w <- to_w(u)
    as.vector(Amat %*% w - rhsvec)
  }
  eval_jac_g_eq <- function(u) {
    w <- to_w(u)
    # Jacobian rows i, columns j: Amat[i, j] * w_j
    Amat %*% diag(w, nrow = K, ncol = K)
  }

  # ---- Inequality constraint: g_ineq(u) = u' Kinv u - L^2 <= 0 ----
  eval_g_ineq <- function(u) {
    as.numeric(t(u) %*% Kinv %*% u - L^2)
  }
  eval_jac_g_ineq <- function(u) {
    # 1 x K row
    matrix(2 * as.vector(Kinv %*% u), nrow = 1)
  }

  # ---- Bounds on u (optional for numerical stability) ----
  # Using moderately wide bounds to avoid overflow in exp(u+1)
  max_abs_u <- 40
  lb <- rep(-max_abs_u, K)
  ub <- rep( max_abs_u, K)

  # ---- Solve with SLSQP on u ----
  # Start at u = 0 (i.e., w = exp(1)), which satisfies the inequality (u=0 ⇒ u'Kinv u = 0 <= L^2)
  x0 <- rep(0, K)
  opts <- list(algorithm  = "NLOPT_LD_SLSQP",
               xtol_rel   = 1e-8,
               ftol_rel   = 1e-8,
               maxeval    = 10000,
               print_level = 0,
               check_derivatives = FALSE
  )

  res <- try(
             nloptr::nloptr(
                            x0              = x0,
                            eval_f          = eval_f,
                            eval_grad_f     = eval_grad_f,
                            lb              = lb,
                            ub              = ub,
                            eval_g_ineq     = eval_g_ineq,
                            eval_jac_g_ineq = eval_jac_g_ineq,
                            eval_g_eq       = eval_g_eq,
                            eval_jac_g_eq   = eval_jac_g_eq,
                            opts            = opts
                            ),
             silent = TRUE
  )

  # Optional fallback to derivative-free COBYLA if SLSQP fails
  use_cobyla_if_needed <- TRUE
  if (inherits(res, "try-error") && isTRUE(use_cobyla_if_needed)) {
    message("Using COBYLA fallback")
    res <- nloptr::nloptr(
                          x0              = x0,
                          eval_f          = eval_f,
                          lb              = lb,
                          ub              = ub,
                          eval_g_ineq     = eval_g_ineq,
                          eval_g_eq       = eval_g_eq,
                          # COBYLA ignores analytic jacobians
                          opts            = list(algorithm   = "NLOPT_LN_COBYLA",
                                                 xtol_rel    = 1e-9,
                                                 maxeval     = 10000,
                                                 print_level = 0)
    )
  }

  # ---- Parse result ----
  is_valid_solution <- is.list(res) && !is.null(res$status) && (as.numeric(res$status) %in% 3:4)

  if (!is_valid_solution) {
    print("Failure.")
    message("nloptr did not report success. Status: ", if (!is.null(res$status)) res$status,
            " | Message: ", if (!is.null(res$message)) res$message)
  }

  u_hat   <- if (is_valid_solution) res$solution else rep(NA_real_, K)
  wts_min <- if (is_valid_solution) to_w(u_hat) else rep(NA_real_, K)
  objval  <- if (is_valid_solution) res$objective else NA_real_

  if (is_valid_solution) {
    # Zero-out tiny weights for neatness
    wts_min[abs(wts_min) < 1e-8] <- 0
  }

  thetahat_min <- if (is_valid_solution) objval - control_mean else NA_real_

  list(theta_min         = thetahat_min,
       wts_min           = wts_min,
       dual_lambda       = NA,     # not computed in this nloptr version
       dual_slack        = NA,     # not computed in this nloptr version
       basis_indices     = NA,     # not computed in this nloptr version
       X                 = sample_space,
       prop              = categorized$prop,
       mean_Y            = categorized$mean_Y,
       is_valid_solution = is_valid_solution,
       L                 = L,
       sigm              = sigm,
       nloptr_status     = if (!is.null(res$status)) res$status else NA_integer_,
       nloptr_message    = if (!is.null(res$message)) res$message else NA_character_
  )
}


# Find theta_min estimate based on data
estimate_theta_min_rkhs_log <- function(dat, L, sigm) {
  trt <- dat$trt
  target <- dat$target
  control_mean <- dat$control_mean
  stopifnot(unique(target$xname) == "X")
  stopifnot(identical(target$type, c("first", "second")))
  # Construct categorized data
  categorized <- make_table(trt)
  # Estimate theta_min
  sol_categorized <- find_theta_min_table_rkhs_log(categorized = categorized,
                                                   control_moments = target$value,
                                                   control_mean = control_mean,
                                                   L = L,
                                                   sigm = sigm
  )
  sol_categorized
}

# Find true theta_min based on population information
true_theta_min_rkhs_log <- function(population, L, sigm) {
  find_theta_min_table_rkhs_log(categorized = population$categorized,
                                control_moments = population$control_moments,
                                control_mean = population$control_mean,
                                L = L,
                                sigm = sigm
  )
}

#' Utility function to summarize output from theta computation using Gaussian
summarize_theta_rkhs_log <- function(result) {
  cat("Minimum: ", result$theta_min, "\n")
  par(mfrow = c(2, 1))
  xvalues <- result$X
  # Compute the weights
  wts_min <- result$wts_min
  plot(xvalues, wts_min)
  # Also plot the optimizer distributions
  dens0min <- wts_min * result$prop
  plot(xvalues, dens0min)
  # Check if the constraints are satisfied
  for (j in 0:2) {
    cat(sum(wts_min * result$prop * xvalues ^ j))
    cat("\n")
  }
}

