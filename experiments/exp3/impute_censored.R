# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(survival)

impute_censored <- function(df,
                            method = c("simulate", "median", "emin_tau"),
                            tau = NULL,
                            seed = NULL) {
  method <- match.arg(method)

  if (method == "emin_tau" && is.null(tau)) {
    stop("For method = 'emin_tau', you must provide a numeric tau.")
  }
  if (!is.null(tau) && (!is.numeric(tau) || length(tau) != 1 || tau <= 0)) {
    stop("tau must be a single positive numeric value.")
  }

  fit <- survreg(Surv(os, os_status) ~ age + sex + kps + isGTR + mgmt,
                 data = df, dist = "weibull")

  sigma <- fit$scale
  mu_vec <- predict(fit, newdata = df, type = "lp")  # linear predictor (log-scale)
  shape <- 1 / sigma                                 # Weibull shape k
  scale_vec <- exp(mu_vec)                            # Weibull scale lambda

  n <- nrow(df)
  if (!is.null(seed)) set.seed(seed)

  # helper: E[min(T, tau) | T > y] for Weibull(shape=k, scale=lambda)
  emin_tau_cond_gt <- function(y, tau, k, lambda) {
    if (tau <= y) return(tau)

    ua <- (y / lambda)^k
    ub <- (tau / lambda)^k

    Sy <- exp(-ua)
    Stau <- exp(-ub)

    # integral of t f(t) dt from y to tau:
    # lambda * Gamma(1 + 1/k) * (pgamma(ub, s) - pgamma(ua, s)), s = 1 + 1/k
    s <- 1 + 1 / k
    int_tft <- lambda * gamma(s) * (pgamma(ub, shape = s, rate = 1) -
                                    pgamma(ua, shape = s, rate = 1))

    (int_tft + tau * Stau) / Sy
  }

  if (method %in% c("simulate", "median")) {
    event_time <- numeric(n)

    for (i in seq_len(n)) {
      if (df$os_status[i] == 1) {
        event_time[i] <- df$os[i]
      } else {
        F_t0 <- pweibull(df$os[i], shape = shape, scale = scale_vec[i])

        if (method == "median") {
          u_med_cond <- F_t0 + 0.5 * (1 - F_t0)
          event_time[i] <- qweibull(u_med_cond, shape = shape, scale = scale_vec[i])
        } else {
          u <- runif(1, min = F_t0, max = 1)
          event_time[i] <- qweibull(u, shape = shape, scale = scale_vec[i])
        }
      }
    }

    df$event_time <- event_time
    return(df)
  }

  # method == "emin_tau"
  emin_tau <- numeric(n)
  for (i in seq_len(n)) {
    y <- df$os[i]
    if (df$os_status[i] == 1) {
      emin_tau[i] <- min(y, tau)
    } else {
      emin_tau[i] <- emin_tau_cond_gt(y = y, tau = tau, k = shape, lambda = scale_vec[i])
    }
  }

  df$emin_tau <- emin_tau
  df
}

