source("../utils/tabulate.R", chdir = TRUE)
source("../utils/my_ROI_solver.R", chdir = TRUE)

# Function to create A, b, c, direction based on categorized data
create_lp_vars_gaussian <- function(categorized, control_moments, control_mean, grid_points, L) {
  sample_space <- categorized$X
  K <- length(grid_points)
  # Sort data based on sample_space
  ord <- seq_along(sample_space)
  if (!is.ordered(sample_space)) {
    ord <- order(sample_space)
    sample_space <- sample_space[ord]
    categorized <- categorized[ord, ]
  }
  # Construct LP s.t. Ax <= rhs
  nconstr <- 3
  Amat <- matrix(NA, ncol = K, nrow = nconstr)
  rhsvec <- rep(NA, nconstr)
  dirvec <- rep(NA, nconstr)
  # Constraints: \sum_x x^j w(x) f_1(x) = mu_j
  for (j in 0:2) {
    for (k in 1:K) {
      vec <- sample_space^j * categorized$prop *
        dnorm(sample_space,
              mean = grid_points[k],
              sd = 1 / L)
      Amat[j + 1, k] <- sum(vec)
    }
  }
  rhsvec[1] <- 1
  dirvec[1] <- "=="
  rhsvec[2] <- control_moments[1]
  dirvec[2] <- "=="
  rhsvec[3] <- control_moments[2]
  dirvec[3] <- "=="
  s <- 3
  stopifnot(s == nconstr)
  stopifnot(!any(is.na(Amat)))
  stopifnot(!any(is.na(dirvec)))
  stopifnot(!any(is.na(rhsvec)))
  # Condition number of Amat
  # cat(kappa(Amat), "\n")
  # Construct the objective function
  # \sum_x w(x) E(Y | X = x) f_1(x)
  cvec <- rep(NA, K)
  for (k in seq_len(K)) {
    vec <- categorized$mean_Y * categorized$prop *
      dnorm(sample_space,
            mean = grid_points[k],
            sd = 1 / L)
    cvec[k] <- sum(vec)
  }
  list(Amat = Amat, rhsvec = rhsvec, dirvec = dirvec, cvec = cvec)
}

# Function to compute theta_min and wts_min based on categorized data
# Input contains:
# Experimental: for each level of X, the proportions and mean(Y) within the level
# Control: Summary statistics of X and mean(Y)
# Centers of Gaussian kernels: grid_points,
# Regularization: L, the inverse of sd of the Gausiian components
find_theta_min_table_gaussian <- function(categorized, control_moments, control_mean, grid_points, L) {
  sample_space <- categorized$X
  K <- length(grid_points)
  # Sort data based on sample_space
  ord <- seq_along(sample_space)
  if (!is.ordered(sample_space)) {
    ord <- order(sample_space)
    sample_space <- sample_space[ord]
    categorized <- categorized[ord, ]
  }
  # Construct LP s.t. Ax <= rhs
  nconstr <- 3
  lp_vars <- create_lp_vars_gaussian(categorized, control_moments, control_mean, grid_points, L)
  Amat <- lp_vars$Amat
  dirvec <- lp_vars$dirvec
  rhsvec <- lp_vars$rhsvec
  cvec <- lp_vars$cvec
  # Solve the optimization
  lp <- OP(objective = cvec,
           constraints = L_constraint(L = Amat,
                                      dir = dirvec,
                                      rhs = rhsvec
           ),
           maximum = FALSE
  )
  model <- my_ROI_solver(lp)
  is_valid_solution <- (model$status$code == 0)
  if (!is_valid_solution) {print(model$status$msg)}

  # Store solution
  alpha_min <- model$solution
  thetahat_min <- model$objval - control_mean
  wts_min <- rep(NA, length(sample_space))
  for (i in seq_along(sample_space)) {
    wts_min[i] <- sum(alpha_min * dnorm(sample_space[i],
                                        mean = grid_points,
                                        sd = 1 / L))
  }
  if (is_valid_solution) {
    # Make small entries zero
    alpha_min[abs(alpha_min) < 1e-7] <- 0
    if (sum(alpha_min != 0) > nconstr) {
      print("Warning: Solution might not be a vertex.")
      print(sum(alpha_min != 0))
    }
    # Make and solve dual LP
    # dual_lp <- make_dual_lp(Amat, rhsvec, cvec, dirvec)
    # dual_lambda <- my_ROI_solver(dual_lp)$solution
    # In the dual problem, we typically have 3 positive entries
    # idx <- order(abs(alpha_min), decreasing = TRUE)[1:nconstr] # indices of top-3 values
    idx <- (alpha_min != 0) # indices of non-zero values
    dual_lambda <- qr.solve(t(Amat[, idx]), cvec[idx])
    dual_lambda <- as.numeric(dual_lambda)
    dual_slack <- cvec - t(Amat) %*% as.numeric(dual_lambda)
    dual_slack <- as.numeric(dual_slack)
  } else {
    dual_lambda <- NA
    dual_slack <- NA
    idx <- NA
  }
  list(theta_min = thetahat_min,
       alpha_min = alpha_min,
       wts_min = wts_min,
       dual_lambda = dual_lambda,
       dual_slack = dual_slack,
       basis_indices = idx,
       X = sample_space,
       prop = categorized$prop,
       mean_Y = categorized$mean_Y,
       is_valid_solution = is_valid_solution,
       grid_points = grid_points,
       L = L
  )
}

# Find theta_min estimate based on data
estimate_theta_min_gaussian <- function(dat, grid_points, L) {
  trt <- dat$trt
  target <- dat$target
  control_mean <- dat$control_mean
  stopifnot(unique(target$xname) == "X")
  stopifnot(identical(target$type, c("first", "second")))
  # Construct categorized data
  categorized <- make_table(trt)
  # Estimate theta_min
  sol_categorized <- find_theta_min_table_gaussian(categorized = categorized,
                                                   control_moments = target$value,
                                                   control_mean = control_mean,
                                                   grid_points = grid_points,
                                                   L = L)
  sol_categorized
}

# Find true theta_min based on population information
true_theta_min_gaussian <- function(population, grid_points, L) {
  find_theta_min_table_gaussian(categorized = population$categorized,
                                control_moments = population$control_moments,
                                control_mean = population$control_mean,
                                grid_points = grid_points,
                                L = L
  )
}

#' Utility function to summarize output from theta computation using Gaussian
summarize_theta_gaussian <- function(result) {
  cat("Minimum: ", result$theta_min, "\n")
  par(mfrow = c(3, 1))
  xvalues <- result$X
  grid_points <- result$grid_points
  alpha_min <- result$alpha_min
  plot(grid_points, alpha_min)
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

###### DATA SPLITTING ######
#' Usual data-splitting without extra term
splitting_basic_gaussian <- function(dat, grid_points, L) {
  trt <- dat$trt
  control_mean <- dat$control_mean
  K <- 3 # Number of folds
  n1 <- nrow(trt)
  # make (roughly) balanced fold assignments
  fold_id <- sample(rep(seq_len(K), length.out = n1))
  extend_weights <- function(grid_points, xtest, alpha, L) {
    wts <- rep(NA, length(xtest))
    for (j in seq_along(xtest)) {
      vec <- alpha * dnorm(xtest[j], mean = grid_points, sd = 1 / L)
      wts[j] <- sum(vec)
    }
    wts
  }
  theta_min_cv <- rep(NA, K)
  for (k in seq_len(K)) {
    ind_A <- which(fold_id != k)
    trt_A <- trt[ind_A, , drop = FALSE]
    trt_B <- trt[-ind_A, , drop = FALSE]
    # Get minimizer from data A
    dat_A <- dat
    dat_A$trt <- trt_A
    est_A <- estimate_theta_min_gaussian(dat_A, grid_points = grid_points, L = L)
    # Compute weighted average based on data B
    wts_B <- extend_weights(grid_points, trt_B$X, est_A$alpha_min, L)
    theta_min_cv[k] <- mean(wts_B * trt_B$Y) # / mean(wts_B)
  }
  mean(theta_min_cv) - control_mean
}

#' Correction with lagrangian
splitting_lagrange_gaussian <- function(dat, grid_points, L) {
  trt <- dat$trt
  control_mean <- dat$control_mean
  K <- 3 # Number of folds
  n1 <- nrow(trt)
  # make (roughly) balanced fold assignments
  fold_id <- sample(rep(seq_len(K), length.out = n1))
  theta_min_cv <- rep(NA, K)
  for (k in seq_len(K)) {
    ind_A <- which(fold_id != k)
    trt_A <- trt[ind_A, , drop = FALSE]
    trt_B <- trt[-ind_A, , drop = FALSE]
    # Get minimizer from data A
    dat_A <- dat
    dat_A$trt <- trt_A
    est_A <- estimate_theta_min_gaussian(dat_A, grid_points = grid_points, L = L)
    # Compute LP variables based on data B
    lp_vars_B <- create_lp_vars_gaussian(categorized = make_table(trt_B),
                                         control_moments = dat$target$value,
                                         control_mean = dat$target$value,
                                         grid_points = grid_points,
                                         L = L
    )
    prelim <- sum(lp_vars_B$cvec * est_A$alpha_min)
    constraint_mismatch <- as.numeric(lp_vars_B$Amat %*% est_A$alpha_min) - lp_vars_B$rhsvec
    stopifnot(length(est_A$dual_lambda) == length(constraint_mismatch))
    theta_min_cv[k] <- prelim - sum(est_A$dual_lambda * constraint_mismatch)
  }
  mean(theta_min_cv) - control_mean
}

#' Correction with derivative
splitting_derivative_gaussian <- function(dat, grid_points, L) {
  trt <- dat$trt
  control_mean <- dat$control_mean
  K <- 3 # Number of folds
  n1 <- nrow(trt)
  # make (roughly) balanced fold assignments
  fold_id <- sample(rep(seq_len(K), length.out = n1))
  theta_min_cv <- rep(NA, K)
  for (k in seq_len(K)) {
    ind_A <- which(fold_id != k)
    trt_A <- trt[ind_A, , drop = FALSE]
    trt_B <- trt[-ind_A, , drop = FALSE]
    # Get minimizer from data A
    dat_A <- dat
    dat_A$trt <- trt_A
    est_A <- estimate_theta_min_gaussian(dat_A, grid_points = grid_points, L = L)
    lp_vars_A <- create_lp_vars_gaussian(categorized = make_table(trt_A),
                                         control_moments = dat$target$value,
                                         control_mean = dat$control_mean,
                                         grid_points = grid_points,
                                         L = L
    )
    Amat_A <- lp_vars_A$Amat
    cvec_A <- lp_vars_A$cvec
    idx_A <- est_A$basis_indices
    # Compute LP variables based on data B
    lp_vars_B <- create_lp_vars_gaussian(categorized = make_table(trt_B),
                                         control_moments = dat$target$value,
                                         control_mean = dat$control_mean,
                                         grid_points = grid_points,
                                         L = L
    )
    if (sum(idx_A) > length(lp_vars_B$rhsvec)) {
      print("WARNING: Solution supported on more than 3 points (NA return)")
      print(length(idx_A))
      print(length(lp_vars_B$rhsvec))
      return(NA)
    }
    prelim <- sum(lp_vars_B$cvec * est_A$alpha_min)
    constraint_mismatch <- as.numeric(lp_vars_B$Amat %*% est_A$alpha_min) - lp_vars_B$rhsvec
    # mismatch_coef <- as.numeric(qr.solve(t(Amat_A[, idx_A]), lp_vars_B$cvec[idx_A]))
    mismatch_coef <- as.numeric(qr.solve(t(Amat_A[, idx_A]), cvec_A[idx_A]))
    stopifnot(length(mismatch_coef) == length(constraint_mismatch))
    theta_min_cv[k] <- prelim - sum(mismatch_coef * constraint_mismatch)
  }
  mean(theta_min_cv) - control_mean
}

#' Correction with 3*3 system of equations
splitting_inverse_gaussian <- function(dat, grid_points, L) {
  trt <- dat$trt
  control_mean <- dat$control_mean
  K <- 3 # Number of folds
  n1 <- nrow(trt)
  # make (roughly) balanced fold assignments
  fold_id <- sample(rep(seq_len(K), length.out = n1))
  theta_min_cv <- rep(NA, K)
  for (k in seq_len(K)) {
    ind_A <- which(fold_id != k)
    trt_A <- trt[ind_A, , drop = FALSE]
    trt_B <- trt[-ind_A, , drop = FALSE]
    # Get minimizer from data A
    dat_A <- dat
    dat_A$trt <- trt_A
    est_A <- estimate_theta_min_gaussian(dat_A, grid_points = grid_points, L = L)
    idx_A <- est_A$basis_indices
    # Compute LP variables based on data B
    lp_vars_B <- create_lp_vars_gaussian(categorized = make_table(trt_B),
                                         control_moments = dat$target$value,
                                         control_mean = dat$control_mean,
                                         grid_points = grid_points,
                                         L = L
    )
    alpha_min_B <- est_A$alpha_min * 0
    alpha_min_B[idx_A] <- qr.solve(lp_vars_B$Amat[, idx_A], lp_vars_B$rhsvec)
    prelim <- sum(lp_vars_B$cvec * alpha_min_B)
    theta_min_cv[k] <- prelim
  }
  mean(theta_min_cv) - control_mean
}

