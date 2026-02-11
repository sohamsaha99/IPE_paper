rm(list = ls())
# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(survival)
# library(SurvRegCensCov)

# Read GBM data
df_GBM <- read.csv("AVAGLIO_CleanRR.csv", header = TRUE, stringsAsFactors = FALSE)
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
df_GBM <- df_GBM %>% drop_na()
# After above operations, there should be 351 rows and 9 columns
print(dim(df_GBM))

# # Now fit a survival model to impute the os time for censored patients
# cox_model <- coxph(Surv(os, os_status) ~ age + sex + kps + mgmt + eor + kps:eor,
#                    data = df_GBM)
# Filter censored patients
df_censored <- df_GBM %>% filter(os_status == 0) %>%
  select(Patient, age, sex, kps, eor, mgmt, os, pfs) %>%
  rename(last_follow = os)
# # Print first three rows
newdat <- df_censored %>% head(3)
# print(newdat)
# Compute survival probabilities for these patients
# surv_probs <- survfit() # Difficult, since this model works with relative risk

# Fit Weibull regression
weibull_model <- survreg(Surv(os, os_status) ~ age + sex + kps + mgmt + eor + kps:eor,
                         data = df_GBM, dist = "weibull")
summary(weibull_model)
weibull_model_coef <- summary(weibull_model)$coef
# # Adequacy of Weibull model
# WeibullDiag(Surv(os, os_status) ~ sex + kps + mgmt + eor, data = df_GBM)

# Predict Weibull model on newdat
a1 = predict(weibull_model, newdata = newdat, type = "linear")
# Manually compute the linear predictor to check correctness
a2 = model.matrix(~ age + sex + kps + mgmt + eor + kps:eor, data = newdat) %*% weibull_model_coef
identical(as.numeric(a1), as.numeric(a2))

# Now compute the survival probability
rweib_shape <- 1 / summary(weibull_model)$scale
rweib_scale <- exp(predict(weibull_model, newdata = newdat, type = "linear"))

# Find quantiles of the Weibull distribution from base R functions
for (i in seq_along(rweib_scale)) {
  print(qweibull(p = c(0.1, 0.5, 0.9), shape = rweib_shape, scale = rweib_scale[i]))
}
# Find quantiles from the built-in predict.survreg method
print(predict(weibull_model, newdata = newdat, type = "quantile", p = c(0.1, 0.5, 0.9)))
# I checked that the last two outputs are identical, therefore the parameterization is correct

# For the censored patients, I want to simulate a value from appropriate Weibull
# distribution, conditional on os >= last_follow.
df_censored$rweib_shape <- 1 / summary(weibull_model)$scale
df_censored$rweib_scale <- as.numeric(exp(predict(weibull_model,
                                                  newdata = df_censored,
                                                  type = "linear")))
# A function to simulate T from Weibull conditional T >= t_obs
r_conditional_weibull <- function(shape, scale, t_obs) {
  # First compute the CDF at the observed value
  p_obs <- pweibull(q = t_obs, shape = shape, scale = scale)
  # print(p_obs)
  # Then sample a value uniformly from [p_obs, 1]
  p_rand <- runif(1, min = p_obs, max = 1)
  # Then compute the quantile at probability p_rand
  c(p_obs, p_rand, qweibull(p_rand, shape = shape, scale = scale))
}
# Sample from the conditional distribution for each patient
set.seed(20250319)
df_censored$os_simulated <- NA
df_censored$p_obs <- NA
df_censored$p_rand <- NA
for (i in 1:nrow(df_censored)) {
  temp <- r_conditional_weibull(shape = df_censored$rweib_shape[i],
                                scale = df_censored$rweib_scale[i],
                                t_obs = df_censored$last_follow[i])
  df_censored$p_obs[i] <- temp[1]
  df_censored$p_rand[i] <- temp[2]
  df_censored$os_simulated[i] <- temp[3]
}

# Now make a complete dataframe combining df_censored and df_GBM
df_complete <- df_GBM %>%
  left_join(df_censored %>% select(Patient, os_simulated), by = "Patient") %>%
  mutate(event_time = ifelse(os_status == 0, os_simulated, os))

df_complete %>% write.csv(file = "GBM_complete.csv", row.names = FALSE)
