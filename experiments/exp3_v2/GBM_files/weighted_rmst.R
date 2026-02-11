library(survival)

# Calculate Restricted Mean Survival Time (RMST) using weights for patients
# Modified from R function simtrial:::rmst_single_arm to include weights
weighted_rmst <- function(time_var, event_var, time_cutoff, weights = NULL) {
  # Input type check
  # check_args(time_var, type = c("integer", "numeric"))
  # check_args(event_var, type = c("integer", "numeric"))
  # check_args(time_cutoff, type = c("integer", "numeric"), length = 1)
  
  # Input value check
  stopifnot(time_var >= 0)
  stopifnot(event_var %in% c(0, 1))
  stopifnot(time_cutoff >= 0)
  
  # Set weights to 1 if NULL
  if (is.null(weights)) {
    weights <- rep(1, length(time_var))
  }
  # Normalize weights so total of weights is equal to sample size
  weights <- length(time_var) * weights / sum(weights)
  stopifnot(length(weights) == length(time_var))
  
  # Fit a single Kaplan-Meier curve
  fit <- survival::survfit(survival::Surv(time_var, event_var) ~ 1,
                           weights = weights)
  
  # Extract survival probability, number of event, number at risk,
  # and number of censored along with observed time from the fitted model
  # as a new data frame.
  df <- data.frame(
    time = fit$time,
    n_risk = fit$n.risk,
    n_event = fit$n.event,
    n_censor = fit$n.censor,
    surv = fit$surv, # Estimated survival probability
    stringsAsFactors = FALSE
  )
  
  # Filter df by survival time less or equal to the pre-specified cut-off time point tau
  df_fit1 <- df[df$time <= time_cutoff, ]
  
  # Add initial value of (time, survival) = (0,1) for calculating time difference
  df_fit2 <- rbind(df_fit1, c(0, NA, NA, NA, 1))
  
  # Add cut-off time if no records observed on the pre-specified time point
  if (max(df_fit1$time) != time_cutoff) {
    df_fit2 <- rbind(df_fit2, c(time_cutoff, NA, NA, NA, NA))
  }
  
  # Sort the data frame by time
  df_fit2 <- df_fit2[order(df_fit2$time), ]
  n_event <- df_fit2$n_event
  n_risk <- df_fit2$n_risk
  
  # Calculate the time difference and set the last value as NA
  time_diff <- c(diff((sort(df_fit2$time))), NA)
  
  # Calculate the area under the curve per time interval
  area <- time_diff * df_fit2$surv
  
  # Calculate the inverse cumulated area under the curve per observed time point A_i
  big_a <- rev(c(0, cumsum(rev(area)[-1])))
  
  # Calculation of dev refers to di / Yi * (Yi - di)
  dev <- (n_event / (n_risk * (n_risk - n_event))) * (big_a^2)
  
  # Based on the calculation, create a data frame with below items:
  # rmst is the estimated RMST
  rmst <- sum(area, na.rm = TRUE)
  # std is the standard error of the estimated RMST
  variance <- sum(dev, na.rm = TRUE) * sum(n_event, na.rm = TRUE) / (sum(n_event, na.rm = TRUE) - 1)
  std <- sqrt(variance)
  # lcl and ucl are lower/upper control limit of CIs for RMST
  # lcl <- rmst - stats::qnorm(1 - alpha / 2) * std
  # ucl <- rmst + stats::qnorm(1 - alpha / 2) * std
  event <- sum(n_event, na.rm = TRUE)
  
  ans <- data.frame(
    time_cutoff, rmst, variance, std, event,
    # time_cutoff, rmst, variance, std, lcl, ucl, event,
    stringsAsFactors = FALSE
  )
  
  return(ans)
}


# # Simulate survival data
# set.seed(20250323)
# N <- 200
# event_times <- rexp(N, rate = 1/30)
# censoring_times <- rexp(N, rate = 1/60)
# observed_times <- pmin(event_times, censoring_times)
# observed_status <- as.numeric(event_times <= censoring_times)
# cat("Censoring proportion:", mean(1 - observed_status), "\n")
# 
# # Generate random weights
# wts <- runif(N)
# 
# # Determine RMST cutoff as the 80th percentile
# tau <- qexp(0.80, rate = 1/30)
# 
# # Plot the KM curve
# fit <- survfit(Surv(observed_times, observed_status) ~ 1)
# plot(fit); abline(v = tau, lty = "dashed", col = "red")
# 
# # Find RMST
# rmst_value <- weighted_rmst(observed_times, observed_status, tau, wts)[1, "rmst"]
# cat("RMST obtained by function:", rmst_value, "\n")
# 
# # Now compute weighted averages
# truncated_observed_times <- pmin(observed_times, tau)
# weighted.mean(truncated_observed_times, wts)
# 
# observed_times <- pmin(event_times, censoring_times, tau)
# observed_status <- as.numeric(event_times <= censoring_times | observed_times == tau)
# weighted_rmst(observed_times, observed_status, tau+100000, wts)[1, "rmst"]
# 
# # Example: n = 3
# wts <- runif(3); wts <- wts / sum(wts)
# obs <- sort(runif(3))
# ind <- c(0, 1, 1)
# weighted_rmst(obs, ind, 2, wts)
# # Do it manually
# wts[1]*(wts[2]*obs[2] + wts[3]*obs[3])/(wts[2]+wts[3]) + (wts[2]*obs[2] + wts[3]*obs[3])
# 
# # Example: n = 4
# wts <- runif(4); wts <- wts / sum(wts)
# obs <- sort(runif(4))
# ind <- c(0, 1, 0, 1)
# weighted_rmst(obs, ind, 2, wts)
# # Do it manually
# wts[2]*obs[2] + wts[4]*obs[4] + wts[3]*obs[4] + wts[1]*(wts[2]*obs[2] + wts[4]*obs[4] + wts[3]*obs[4])/(1 - wts[1])
# 
# # Example: n = 4. Tied censored times
# wts <- runif(4); wts <- wts / sum(wts)
# obs <- sort(runif(3)); obs <- obs[c(1, 1, 2, 3)]
# ind <- c(0, 0, 1, 1)
# weighted_rmst(obs, ind, 2, wts)
# # Do it manually
# wts[3]*obs[3] + wts[4]*obs[4] + (wts[1] + wts[2]) * (wts[3]*obs[3] + wts[4]*obs[4]) / (wts[3] + wts[4])
# weighted_rmst_recursion(obs, ind, 2, wts)$rmst
# 
# obs <- c(4, 4, 3, 2, 1)
# ind <- c(0, 1, 0, 0, 0)
# ord <- order(obs, 1 - ind)
# obs[ord]; ind[ord]

# Function to compute weighted RMST based on recursive formula
weighted_rmst_recursion <- function(time_var, event_var, time_cutoff, weights = NULL) {
  # Input value check
  stopifnot(time_var >= 0)
  stopifnot(event_var %in% c(0, 1))
  stopifnot(time_cutoff >= 0)
  
  N <- length(time_var)
  # Set weights to 1 if NULL
  if (is.null(weights)) {
    weights <- rep(1, N)
  }
  # Normalize weights so total of weights is equal to sample size
  weights <- length(time_var) * weights / sum(weights)
  stopifnot(length(weights) == N)
  
  # Apply truncation at cutoff and change event status to 1 if truncation is applied
  event_var <- ifelse(time_var >= time_cutoff, 1, event_var)
  time_var <- pmin(time_var, time_cutoff)
  
  # If there are no censored observations, simply return the weighted average
  if (all(event_var == 1)) {
    return(list(
      rmst = weighted.mean(x = time_var, w = weights),
      T_hat = time_var
    ))
  }
  
  # Sort indices so that U_i is in increasing order
  # In case of ties, death events are placed before censored events
  ord <- order(time_var, 1 - event_var)
  time_var <- time_var[ord]
  event_var <- event_var[ord]
  weights <- weights[ord]
  # Make sure that the final observation has delta = 1
  stopifnot(event_var[N] == 1)
  
  # Compute the conditional expectations \hat{T}_i = \hat{E}(U_i | U_i, \delta_i),
  # where \hat{E} refers to expectation w.r.t. the weighted Kaplan-Meier.
  #
  # For death events, \hat{T}_i = U_i
  T_hat <- ifelse(event_var == 1, time_var, NA)
  # Find all i such that delta_i = 0
  cens_ind <- which(is.na(T_hat))
  # Reverse the list so we can traverse from right
  cens_ind <- cens_ind[length(cens_ind):1]
  # Start from the right
  for (ind in cens_ind) {
    # Weighted average of subsequent \hat{T}_i's
    T_hat[ind] <- weighted.mean(T_hat[(ind + 1):N], weights[(ind + 1):N])
  }
  stopifnot(all(!is.na(T_hat)))
  # Find T_hat in the original order
  T_hat_original_order <- rep(NA, N)
  T_hat_original_order[ord] <- T_hat
  return(list(
    rmst = weighted.mean(x = T_hat, w = weights),
    T_hat = T_hat_original_order
  ))
}

# ##################################
# # Simulate survival data
# set.seed(20250323)
# N <- 200
# event_times <- rexp(N, rate = 1/30)
# censoring_times <- rexp(N, rate = 1/60)
# observed_times <- pmin(event_times, censoring_times)
# observed_status <- as.numeric(event_times <= censoring_times)
# cat("Censoring proportion:", mean(1 - observed_status), "\n")
# 
# # Generate random weights
# wts <- runif(N)
# 
# # Determine RMST cutoff as the 80th percentile
# tau <- qexp(0.80, rate = 1/30)
# 
# # Plot the KM curve
# fit <- survfit(Surv(observed_times, observed_status) ~ 1)
# plot(fit); abline(v = tau, lty = "dashed", col = "red")
# 
# # Find RMST
# rmst_value <- weighted_rmst(observed_times, observed_status, tau, wts)[1, "rmst"]
# cat("RMST obtained by function:", rmst_value, "\n")
# 
# # Now compute RMST by my function
# rmst_value_2 <- weighted_rmst_recursion(observed_times, observed_status, tau, wts)$rmst
# cat("RMST obtained by my function:", rmst_value_2, "\n")
