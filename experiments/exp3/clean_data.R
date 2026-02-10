rm(list = ls())
# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(survival)
# library(SurvRegCensCov)

# Read GBM data
df_GBM <- read.csv("./experiments/exp3/GBM_files/AVAGLIO_CleanRR.csv", header = TRUE, stringsAsFactors = FALSE)
# Only keep necessary columns
df_GBM <- df_GBM %>%
  select(Patient, age, sex, kps90_100, eor, mgmt, os_status, os, pfs) %>%
  rename(kps = kps90_100) %>% # Binary: 0: 0-80, 1: 90-100
  filter(eor %in% c("GTR", "STR")) %>%
  mutate(eor = factor(eor)) %>% # Levels: BX, GTR, STR
  mutate(os = os / 30, pfs = pfs / 30) # Measure `os` in months instead of days
# Additionally, we note that os_status = 0 implies alive at last follow-up,
# sex = 0 implies female
# Delete rows with NA (missing) entry
df_GBM <- df_GBM %>% drop_na(age, sex, kps, eor, mgmt, os_status, os)
# After above operations, there should be 351 rows and 9 columns
print(dim(df_GBM))


# Function for generating event time for censored patients
add_event_time_conditional <- function(df, method = c("simulate","median"), seed = NULL) {
  method <- match.arg(method)
  fit <- survreg(Surv(os, os_status) ~ age + sex + kps + eor + mgmt,
                 data = df, dist = "weibull")
  sigma <- fit$scale               # survreg scale
  mu_vec <- predict(fit, newdata = df, type = "lp")  # linear predictor (log-scale)
  shape_vec <- 1 / sigma
  scale_vec <- exp(mu_vec)
  n <- nrow(df)
  if (!is.null(seed)) set.seed(seed)

  event_time <- numeric(n)

  for (i in seq_len(n)) {
    if (df$os_status[i] == 1) {
      # observed event: keep observed time
      event_time[i] <- df$os[i]
    } else {
      # censored: sample from distribution conditional on T > os
      F_t0 <- pweibull(df$os[i], shape = shape_vec, scale = scale_vec[i])
      if (method == "median") {
        # get median of the conditional distribution: invert at u = (F_t0 + 1)/2
        u_med_cond <- F_t0 + 0.5 * (1 - F_t0)
        event_time[i] <- qweibull(u_med_cond, shape = shape_vec, scale = scale_vec[i])
      } else {
        # simulate from Uniform(F_t0, 1)
        u <- runif(1, min = F_t0, max = 1)
        event_time[i] <- qweibull(u, shape = shape_vec, scale = scale_vec[i])
      }
    }
  }

  df$event_time <- event_time
  df
}

df_complete <- add_event_time_conditional(df_GBM, method = "median", seed = 20260128)

df_complete %>% write.csv(file = "./experiments/exp3/GBM_files/GBM_complete.csv", row.names = FALSE)



