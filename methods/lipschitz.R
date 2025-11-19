source("../utils/tabulate.R", chdir = TRUE)
source("../utils/my_ROI_solver.R", chdir = TRUE)
source("./MAIC.R", chdir = TRUE, local = TRUE)

# Function to create A, b, c, direction based on categorized data
create_lp_vars_lipschitz <- function(categorized, control_moments, control_mean, baseline, L) {
  sample_space <- categorized$X
  K <- length(sample_space)
  # Check if baseline length is valid
  stopifnot(length(baseline) == K)
  # Sort data based on sample_space
  ord <- seq_along(sample_space)
  if (!is.ordered(sample_space)) {
    ord <- order(sample_space)
    sample_space <- sample_space[ord]
    categorized <- categorized[ord, ]
    baseline <- baseline[ord]
  }
  # Construct LP s.t. Ax <= rhs
  nconstr <- 3 + 2 * (K - 1)
  Amat <- matrix(NA, ncol = K, nrow = nconstr)
  rhsvec <- rep(NA, nconstr)
  dirvec <- rep(NA, nconstr)
  # Constraints: \sum_x x^j w(x) f_1(x) = mu_j
  Amat[1, ] <- categorized$prop
  rhsvec[1] <- 1
  dirvec[1] <- "=="
  Amat[2, ] <- categorized$prop * sample_space
  rhsvec[2] <- control_moments[1]
  dirvec[2] <- "=="
  Amat[3, ] <- categorized$prop * sample_space ^ 2
  rhsvec[3] <- control_moments[2]
  dirvec[3] <- "=="
  s <- 3
  # Regularity constraints
  for (k in seq_len(K - 1)) {
    Amat[s + k, ] <- 0
    dist <- abs(categorized$X[k + 1] - categorized$X[k])
    Amat[s + k, c(k, k + 1)] <- c(baseline[k + 1],
                                  - baseline[k] * exp(L * dist))
    rhsvec[s + k] <- 0
    dirvec[s + k] <- "<="
  }
  s <- s + K - 1
  for (k in seq_len(K - 1)) {
    Amat[s + k, ] <- 0
    dist <- abs(categorized$X[k + 1] - categorized$X[k])
    Amat[s + k, c(k, k + 1)] <- c(- baseline[k + 1],
                                  baseline[k] * exp(- L * dist))
    rhsvec[s + k] <- 0
    dirvec[s + k] <- "<="
  }
  s <- s + K - 1
  stopifnot(s == nconstr)
  stopifnot(!any(is.na(Amat)))
  stopifnot(!any(is.na(dirvec)))
  stopifnot(!any(is.na(rhsvec)))
  # Condition number of Amat
  # cat(kappa(Amat), "\n")
  # Construct the objective function
  # \sum_x w(x) E(Y | X = x) f_1(x)
  cvec <- categorized$mean_Y * categorized$prop
  list(Amat = Amat, rhsvec = rhsvec, dirvec = dirvec, cvec = cvec)
}

# Function to compute theta_min and wts_min based on categorized data
# Input contains:
# Experimental: for each level of X, the proportions and mean(Y) within the level
# Control: Summary statistics of X and mean(Y)
# Regularization: L, the Lipschitz constant for log (w / w_{MAIC})
find_theta_min_table_lipschitz <- function(categorized, control_moments, control_mean, baseline, L) {
  sample_space <- categorized$X
  # Sort data based on sample_space
  ord <- seq_along(sample_space)
  if (!is.ordered(sample_space)) {
    ord <- order(sample_space)
    sample_space <- sample_space[ord]
    categorized <- categorized[ord, ]
  }
  # Construct LP s.t. Ax <= rhs
  nconstr <- 3
  lp_vars <- create_lp_vars_lipschitz(categorized, control_moments, control_mean, baseline, L)
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
  wts_min <- model$solution
  thetahat_min <- model$objval - control_mean
  if (is_valid_solution) {
    # Make and solve dual LP
    # dual_lp <- make_dual_lp(Amat, rhsvec, cvec, dirvec)
    # dual_lambda <- my_ROI_solver(dual_lp)$solution
    # In the dual problem, we typically have 3 positive entries
    # idx <- (alpha_min != 0) # indices of non-zero values
    # dual_lambda <- qr.solve(t(Amat[, idx]), cvec[idx])
    # dual_lambda <- as.numeric(dual_lambda)
    # dual_slack <- cvec - t(Amat) %*% as.numeric(dual_lambda)
    # dual_slack <- as.numeric(dual_slack)
    dual_lambda <- NA
    dual_slack <- NA
    idx <- NA
  } else {
    dual_lambda <- NA
    dual_slack <- NA
    idx <- NA
  }
  list(theta_min = thetahat_min,
       wts_min = wts_min,
       dual_lambda = dual_lambda,
       dual_slack = dual_slack,
       basis_indices = idx,
       X = sample_space,
       prop = categorized$prop,
       mean_Y = categorized$mean_Y,
       is_valid_solution = is_valid_solution,
       L = L
  )
}

# Find theta_min estimate based on data
estimate_theta_min_lipschitz <- function(dat, L) {
  trt <- dat$trt
  target <- dat$target
  control_mean <- dat$control_mean
  stopifnot(unique(target$xname) == "X")
  stopifnot(identical(target$type, c("first", "second")))
  # Construct categorized data
  categorized <- make_table(trt)
  # Find baseline from MAIC
  eta_hat <- MAIC_categorized(categorized, target$value)
  baseline <- exp(eta_hat[1] * categorized$X + eta_hat[2] * categorized$X^2)
  baseline <- baseline / sum(baseline * categorized$prop) # \sum w_i p_i = 1
  # Estimate theta_min
  sol_categorized <- find_theta_min_table_lipschitz(categorized = categorized,
                                                    control_moments = target$value,
                                                    control_mean = control_mean,
                                                    baseline = baseline,
                                                    L = L)
  sol_categorized
}

# Find true theta_min based on population information
true_theta_min_lipschitz <- function(population, L) {
  categorized <- population$categorized
  # Find baseline from MAIC
  eta_star <- MAIC_categorized(categorized, population$control_moments)
  baseline <- exp(eta_star[1] * categorized$X + eta_star[2] * categorized$X^2)
  baseline <- baseline / sum(baseline * categorized$prop) # \sum w_i p_i = 1
  find_theta_min_table_lipschitz(categorized = categorized,
                                 control_moments = population$control_moments,
                                 control_mean = population$control_mean,
                                 baseline = baseline,
                                 L = L
  )
}

#' Utility function to summarize output from theta computation using Gaussian
summarize_theta_lipschitz <- function(result) {
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

