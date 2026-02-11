library(survival)
library(adjustedCurves)

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
  print(w)

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
