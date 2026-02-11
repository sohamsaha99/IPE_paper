# rm(list = ls())
rm(list = ls()[ls() != "Lvalue_global"])
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(latex2exp)
############ Now apply MAIC on the complete data #########
df_GBM <- read.csv("GBM_complete.csv", header = TRUE, stringsAsFactors = FALSE)
ggplot(df_GBM %>%
         mutate(sex = factor(sex),
                kps = factor(kps),
                mgmt = factor(mgmt),
                eor = factor(eor))) +
  geom_point(aes(x = age, y = os, color = sex)) +
  facet_grid(eor ~ kps,
             labeller = labeller(.rows = label_both, .cols = label_both)) +
  theme_pubr(border = TRUE)

df_GBM %>% mutate(healthy = (kps == 1 & eor == "GTR" & age < 55)) %>%
  ggplot(aes(x = healthy, y = event_time, fill = healthy)) +
  geom_boxplot() +
  ylim(c(0, 40)) +
  theme_pubr() +
  scale_fill_discrete() +
  labs(y = "Overall Survival (months)",
       x = "Healthy: kps > 90, eor = \"GTR\", age < 55") +
  theme(legend.position = "none")

# Read MAIC code
source("../methods/MAIC.R", chdir = TRUE, local = TRUE)
source("../methods/find_range_var_of_thetahat.R", chdir = TRUE, local = TRUE)
source("weighted_rmst.R", chdir = TRUE, local = TRUE)
# Define a sampling method

#### STEP 1.1: Function for fitting a Weibull model and ####
## report E(min{T, tau}|Y, delta, X) for each observation ####
predict_weibull <- function(dat, tau) {
  # dat is the data of the experimental arm
  # Fit weibull regression
  weibull_model <- survreg(Surv(os, os_status) ~ age + sex + kps + mgmt + eor + kps:eor,
                           data = dat, dist = "weibull")
  # Now compute the parameters of Weibull distribution
  rweib_shape <- 1 / summary(weibull_model)$scale
  rweib_scale <- as.numeric(exp(predict(weibull_model,
                                            newdata = dat,
                                            type = "linear")))
  # print(rweib_scale)
  # print(rweib_shape)
  conditional_expectation <- rep(NA, nrow(dat))
  for (i in 1:nrow(dat)) {
    if (dat$os_status[i] == 1) {
      conditional_expectation[i] <- min(dat$os[i], tau)
    } else {
      # Find E(min{T, tau}|T > os, X)
      sh <- rweib_shape
      sc <- rweib_scale[i]
      # Find CDF at os
      p_obs <- pweibull(dat$os[i], shape = sh, scale = sc)
      # print(p_obs)
      # Simulate 100 uniform in [p_obs, 1]
      u <- runif(100, p_obs, 1)
      # Find F^{-1}(u) where F is the Weibull CDF
      w <- qweibull(u, shape = sh, scale = sc)
      # Find average of min{w, tau}
      conditional_expectation[i] <- mean(pmin(w, tau))
    }
  }
  conditional_expectation
}

#### STEP 2: Function for generating one in-silico trial ####
generate_data <- function(dat, n_trt, n_ext, tau) {
  # First generate the data for treated group
  # First choose indices depending on the age of individual patients
  ind_trt <- sample(nrow(dat), size = n_trt, replace = TRUE,
                    prob = ifelse(dat$kps == 1 & dat$eor == "GTR" & dat$age < 55, 0.25, 0.75))
  # prob = ifelse(dat$age < 55, 0.25, 0.75))
  # Now choose indices depending on the age of individual patients
  ind_ext <- sample(nrow(dat), size = n_ext, replace = TRUE,
                    prob = ifelse(dat$kps == 1 & dat$eor == "GTR" & dat$age < 55, 0.75, 0.25))
  # prob = ifelse(dat$age < 55, 0.75, 0.25))
  # Construct predictor columns of synthetic dataset
  dat_synth <- dat[c(ind_trt, ind_ext), c("age", "sex", "kps", "mgmt", "eor", "os", "os_status")]
  # Add column of treatment indicator (equivalent to: source of the patient.)
  dat_synth$isTrial <- c(rep(1, n_trt), rep(0, n_ext))
  # Find E(min{T_i, tau} | Y_i, delta_i, X_i) in the treated patients
  dat_synth$Y <- NA
  dat_synth[1:n_trt, "Y"] <- predict_weibull(dat_synth[1:n_trt, ], tau)
  # dat_synth$Y <- dat_synth$event_time + runif(nrow(dat_synth), min = 0, max = 0.5)
  # Scale age to [0, 1]
  dat_synth$age <- (dat_synth$age - min(dat$age)) / (max(dat$age) - min(dat$age))
  return(dat_synth[, c("age", "sex", "kps", "mgmt", "eor", "isTrial", "os", "os_status", "Y")])
}

#### STEP 3: Function for generating treatment arm and summary based on full trial ####
generate_summary <- function(dat, tau) {
  X <- model.matrix(~ age + sex + kps + mgmt + eor, data = dat)
  X <- X[, -1]
  X_trt <- X[dat$isTrial == 1, ]
  X_trt <- data.frame(X_trt)
  X_trt$Y <- dat[dat$isTrial == 1, "Y"]
  X_con <- X[dat$isTrial == 0, ]
  # Prepare summary info of control group
  target <- data.frame(matrix(nrow = 6, ncol = 3))
  colnames(target) <- c("xname", "value", "type")
  target$xname <- c("age", "age", "sex", "kps", "mgmt", "eorSTR")
  target$value <- c(mean(X_con[, "age"]), mean(X_con[, "age"]^2),
                    mean(X_con[, "sex"]),
                    mean(X_con[, "kps"]),
                    mean(X_con[, "mgmt"]),
                    mean(X_con[, "eorSTR"]))
  target$type <- c("first", "second", "first", "first", "first", "first")
  # The cnotrol mean will be the RMST for the control group
  control_mean <- weighted_rmst_recursion(time_var = dat[dat$isTrial == 0, "os"],
                                          event_var = dat[dat$isTrial == 0, "os_status"],
                                          time_cutoff = tau)$rmst
  # control_mean <- mean(dat[dat$isTrial == 0, "Y"])
  list(
    trt = X_trt, target = target, control_mean = control_mean
  )
}

#### Step 4: Filter synthetic data and summary based on columns to use ####
filter_data_summary <- function(dat_summary, column_names = c("age", "sex", "kps", "mgmt", "eorSTR")) {
  # Filter treatment dataframe
  dat_summary$trt <- dat_summary$trt[, c(column_names, "Y")]
  # Assign numerical values to binary covariates
  scaling_factor <- 1
  binary_columns <- NULL
  for (col in column_names) {
    # Check if the covariate is binary
    if (length(unique(dat_summary$trt[, col])) == 2) {
      binary_columns <- c(binary_columns, col)
      dat_summary$trt[, col] <- dat_summary$trt[, col] * scaling_factor
    }
  }
  # Filter summary target dataframe
  dat_summary$target <- dat_summary$target[dat_summary$target$xname %in% column_names, ]
  # Apply same scaling factor to the binary covariates in the target dataframe
  for (i in 1:nrow(dat_summary$target)) {
    col <- dat_summary$target$xname[i]
    if (col %in% binary_columns) {
      moment <- ifelse(dat_summary$target$type[i] == "first", 1, 2)
      dat_summary$target$value[i] <- dat_summary$target$value[i] * (scaling_factor)^moment
    }
  }
  dat_summary
}

# Make a scenario where MAIC is biased
B <- 1000
maic_estimate <- rep(NA, B)
matching_variables <- c("age", "sex", "kps", "mgmt", "eorSTR")
set.seed(20240728)
n_trt <- 200
n_ext <- 200
library(doParallel)
cl <- makeCluster(12)
registerDoParallel(cl)
tols_optimization <- list(
  eps = 2,
  tau0 = 0.00001,
  tau1p = 0.0001,
  # L = 2,
  L = Lvalue_global,
  distance_metric = "euclidean"
)
result <- foreach (i = 1:B, .packages = c("lpSolveAPI", "dplyr", "survival")) %dopar% {
  dat_synth <- generate_data(df_GBM, n_trt, n_ext, tau = 24)
  dat_synth_summary <- generate_summary(dat_synth, tau = 24)
  dat_synth_summary <- filter_data_summary(dat_synth_summary, column_names = matching_variables)
  maic <- MAIC(dat_synth_summary$trt, dat_synth_summary$target, dat_synth_summary$control_mean)
  maic_estimate <- maic$estimate
  sample_optimal <- find_range_of_thetahat(dat_synth_summary$trt,
                                           dat_synth_summary$target,
                                           dat_synth_summary$control_mean,
                                           tols_optimization)
  thetahat_min <- sample_optimal$thetahat_min
  thetahat_max <- sample_optimal$thetahat_max
  # weights_min[, i] <- sample_optimal$wts_min
  # weights_max[, i] <- sample_optimal$wts_max
  c(maic_estimate, thetahat_min, thetahat_max)
}
stopCluster(cl)
# saveRDS(result, file = "GBM_1000_20240319.RDS")
# Make a dataframe of results
df_result <- data.frame(matrix(NA, nrow = B, ncol = 3))
colnames(df_result) <- c("maic", "thetahat_min", "thetahat_max")
for (i in 1:B) {
  df_result[i, ] <- result[[i]]
}
df_result <- df_result %>% pivot_longer(cols = everything(),
                                        names_to = "method",
                                        values_to = "estimate")

p1 <- ggplot(df_result) +
  # geom_histogram(aes(x = estimate, fill = method), alpha = 0.6, position = "identity") +
  geom_density(aes(x = estimate, fill = method), alpha = 0.6, position = "identity") +
  scale_fill_discrete(labels = c(
    TeX(r"($\hat{\theta}_{MAIC}$)"),
    TeX(r"($\hat{\theta}_{max}$)"),
    TeX(r"($\hat{\theta}_{min}$)")
  )) +
  # xlim(c(-20, 30)) +
  theme_bw() +
  labs(x = "Estimated Difference in Restricted Mean Survival Time (months)",
       y = "Density")

p2 <- ggplot(df_result) +
  geom_boxplot(aes(y = method, x = estimate, fill = method)) +
  scale_fill_discrete(labels = c(
    TeX(r"($\hat{\theta}_{MAIC}$)"),
    TeX(r"($\hat{\theta}_{max}$)"),
    TeX(r"($\hat{\theta}_{min}$)")
  )) +
  labs(x = "",
       y = "Method") +
  # ylim(c(-20, 30)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  scale_x_continuous(position = "top")

print(ggarrange(p1, p2, nrow = 2, common.legend = TRUE, legend = "top", align = "v"))

# Save results in a dataframe
saveRDS(
  list(tols_optimization = tols_optimization,
       result = df_result),
  file = paste0("results/L_value_", tols_optimization$L, ".RDS")
)
# print(summary(maic_estimate))
# maic_estimate <- data.frame(x = maic_estimate)
# par(mfrow = c(2, 1), mar = c(4,4,2,1))
# p <- ggplot(maic_estimate, aes(x = x)) +
#   geom_density(fill = "steelblue", alpha = 0.5) +
#   geom_boxplot(aes(y = -0.05), width = 0.05, fill = "lightgray", outlier.color = "red") +
#   theme_minimal() +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   labs(x = "Estimated Survival Gain (months)",
#        y = "Density",
#        title = "MAIC") +
#   theme(axis.text.y = element_blank(),
#         axis.ticks.y = element_blank())
# 
# plot(p)

# ggplot(maic_estimate, aes(x = x, y = -0.5)) +
#   geom_boxplot(fill = "lightgray") +
#   geom_density(aes(x = x), inherit.aes = FALSE) +
#   stat_boxplot(geom = "vline", aes(xintercept = ..xlower..)) +
#   stat_boxplot(geom = "vline", aes(xintercept = ..xmiddle..)) +
#   stat_boxplot(geom = "vline", aes(xintercept = ..xupper..)) +
#   scale_fill_discrete()
# boxplot(maic_estimate); abline(h = 0, lty = "dashed")
# plot(density(maic_estimate, na.rm = T))
