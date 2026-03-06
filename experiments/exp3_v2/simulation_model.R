library(survival)

#' Input: dat with columns age, sex, kps, isGTR, mgmt
#' n: Number of rows in new data
#' Procedure: First fit Weibull for C and T | X.
#' Resample predictors. Generate C, T independently given X.
resample_data <- function(dat, n) {
  required <- c("age", "sex", "kps", "mgmt", "isGTR", "os", "os_status")
  missing_cols <- setdiff(required, names(dat))
  if (length(missing_cols) > 0) {
    stop("Missing columns in dat: ", paste(missing_cols, collapse = ", "))
  }

  # 1) Resample X rows
  newdat <- dat[sample(nrow(dat), n, replace = TRUE), , drop = FALSE]

  # 2) Parametric model for event times: Weibull AFT via survreg
  fit_T <- survreg(Surv(os, os_status) ~ age + sex + kps + isGTR + mgmt,
                   data = dat, dist = "weibull")

  sigma_T <- fit_T$scale
  mu_T    <- predict(fit_T, newdata = newdat, type = "lp")  # log(scale) in Weibull AFT

  shape_T <- 1 / sigma_T
  scale_T <- exp(mu_T)

  uT <- runif(n)
  T  <- scale_T * (-log(uT))^(1 / shape_T)

  # 3) Parametric model for censoring times: fit on "censoring events"
  #    Treat censoring as the event indicator (1=censored, 0=observed event)
  dat$cen_status <- 1L - dat$os_status

  fit_C <- survreg(Surv(os, cen_status) ~ age + sex + kps + isGTR + mgmt,
                   data = dat, dist = "weibull")

  sigma_C <- fit_C$scale
  mu_C    <- predict(fit_C, newdata = newdat, type = "lp")

  shape_C <- 1 / sigma_C
  scale_C <- exp(mu_C)

  uC <- runif(n)
  C  <- scale_C * (-log(uC))^(1 / shape_C)

  # 4) Observed time/status
  time   <- pmin(T, C)
  status <- as.integer(T <= C)  # 1=event observed, 0=censored

  # cbind(newdat, T_event = T, C = C, time = time, status = status)
  newdat$os <- time
  newdat$os_status <- status
  newdat
}
