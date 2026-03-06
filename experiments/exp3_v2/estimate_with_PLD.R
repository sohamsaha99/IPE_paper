library(survival)
library(adjustedCurves)
source("./experiments/exp3_v2/weighted_rmst.R")

#' Given two datasets, find difference in RMST
#' where the treated group is weighted and
#' control group is unweighted
rmst_diff_atc_adjcurves <- function(df_trt, df_con, tau, n_boot = 500) {
  stopifnot(requireNamespace("adjustedCurves", quietly = TRUE))

  preds <- c("age", "sex", "kps", "mgmt", "isGTR")

  trt <- df_trt; con <- df_con
  # We use control indicator because we are doing ATC
  # (Adjusted treatment effect for Controls)
  trt$isCon <- 0L; con$isCon <- 1L
  dat <- rbind(trt, con)

  dat <- dat[stats::complete.cases(dat[, c("isCon", preds, "os", "os_status")]), ]
  dat$isCon <- as.factor(dat$isCon)   # adjustedCurves expects factor groups
  ps_fit <- stats::glm(isCon ~ age + sex + kps + mgmt + isGTR,
                       data = dat, family = stats::binomial())
  pC <- predict(ps_fit, type = "response")
  pC <- pmin(pmax(pC, 1e-6), 1 - 1e-6)

  # ATC weights: controls unweighted (=1), treated weighted to match controls
  w <- rep(1, nrow(dat))
  idx_t <- dat$isCon == 0
  w[idx_t] <- pC[idx_t] / (1 - pC[idx_t])

  # Scaling so total weight = number of treated
  w[idx_t] <- w[idx_t] * (sum(idx_t) / sum(w[idx_t]))
  # print(w)

  adjs <- adjustedCurves::adjustedsurv(
                                       data = dat,
                                       variable = "isCon",
                                       ev_time = "os",
                                       event = "os_status",
                                       method = "iptw_km",
                                       treatment_model = ps_fit, # <-- Weights model
                                       bootstrap = TRUE,
                                       n_boot = n_boot,
                                       conf_int = FALSE
  )

  # RMST contrast (treated - control) up to tau, with bootstrap SE/CI
  out <- adjustedCurves::adjusted_rmst(adjs,
                                       to = tau,
                                       conf_int = TRUE,
                                       contrast = "diff",
                                       group_1 = "0",
                                       group_2 = "1"
  )

  list(
       point_estimate = out$diff[1],
       sd = out$se[1],
       ci_95 = c(out$ci_lower[1], out$ci_upper[1]),
       adjs = adjs
  )
}

#' Rewrite of weighted RMST difference, where CI is based on bootstrap
rmst_diff_atc_manual <- function(df_trt, df_con, tau,
                                 n_boot = 500,
                                 conf_level = 0.95,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)

  # ---- internal: compute point estimate on given data ----
  .est <- function(trt, con) {
    trt$isCon <- 0L
    con$isCon <- 1L
    dat <- rbind(trt, con)

    ps_fit <- stats::glm(isCon ~ age + sex + kps + mgmt + isGTR,
                         data = dat, family = stats::binomial())

    pC <- stats::predict(ps_fit, type = "response")
    pC <- pmin(pmax(pC, 1e-6), 1 - 1e-6)

    # ATC weights: controls = 1, treated = pC/(1-pC), then scale treated weights
    w <- rep(1, nrow(dat))
    idx_t <- dat$isCon == 0   # treated rows in combined data
    w[idx_t] <- pC[idx_t] / (1 - pC[idx_t])
    w[idx_t] <- w[idx_t] * (sum(idx_t) / sum(w[idx_t]))  # scale treated to N_treated

    # IMPORTANT: weights must align with the rows in trt and con
    w_trt <- w[idx_t]
    w_con <- w[!idx_t]

    theta1_hat <- weighted_rmst(trt$os, trt$os_status, tau, weights = w_trt)
    theta0_hat <- weighted_rmst(con$os, con$os_status, tau, weights = w_con)

    theta1_hat$rmst[1] - theta0_hat$rmst[1]
  }

  point <- .est(df_trt, df_con)

  if (is.na(conf_level)) {
    return(list(
                estimate = point,
                conf_level = NA,
                ci = rep(NA, 2),
                boot = rep(NA, n_boot)
                ))
  }
  # ---- bootstrap ----
  n_t <- nrow(df_trt)
  n_c <- nrow(df_con)

  boot_stats <- replicate(n_boot, {
                            trt_b <- df_trt[sample.int(n_t, n_t, replace = TRUE), , drop = FALSE]
                            con_b <- df_con[sample.int(n_c, n_c, replace = TRUE), , drop = FALSE]
                            .est(trt_b, con_b)
                         })

  alpha <- 1 - conf_level
  ci <- stats::quantile(boot_stats, probs = c(alpha/2, 1 - alpha/2), na.rm = TRUE)

  list(
       estimate = point,
       conf_level = conf_level,
       ci = unname(ci),
       boot = boot_stats
  )
}
