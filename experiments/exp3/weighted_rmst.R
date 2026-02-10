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

