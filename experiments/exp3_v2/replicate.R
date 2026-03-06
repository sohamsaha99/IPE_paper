rm(list = ls())

source("./experiments/exp3_v2/estimate_with_PLD.R")
source("./experiments/exp3_v2/estimate_with_summary.R")
source("./experiments/exp3_v2/simulation_model.R")
source("./experiments/exp3_v2/prepare_summary.R")
source("./myggplot_theme.R")

# Read data
df_trt <- read.csv("./experiments/exp3_v2/GBM_files/GBM_trt.csv", header = TRUE, stringsAsFactors = FALSE)
df_con <- read.csv("./experiments/exp3_v2/GBM_files/GBM_con.csv", header = TRUE, stringsAsFactors = FALSE)

# Keep required columns only
print(dim(df_trt))
print(dim(df_con))
df_trt <- df_trt %>%
  mutate(isGTR = as.numeric(eor == "GTR")) %>%
  select(age, sex, kps, mgmt, isGTR, os, os_status)
df_con <- df_con %>%
  mutate(isGTR = as.numeric(eor == "GTR")) %>%
  select(age, sex, kps, mgmt, isGTR, os, os_status)


#####
# # For now we only use DFCI data to make trt and control
# set.seed(20260211)
# df <- read.csv("./experiments/exp3_v2/GBM_files/GBM_con.csv", header = TRUE, stringsAsFactors = FALSE)
# df <- df %>%
#   mutate(isGTR = as.numeric(eor == "GTR")) %>%
#   select(age, sex, kps, mgmt, isGTR, os, os_status)
# idx_t <- sample(nrow(df), round(nrow(df) / 2), replace = FALSE)
# df_trt <- df[idx_t, ]
# df_con <- df[-idx_t, ]
### --- Individual Data prepared --- ###

# Prepare IPE parameters
grid_points <- expand.grid(age = seq(10, 100, by = 1),
                           sex = 0:1,
                           kps = 0:1,
                           mgmt = 0:1,
                           isGTR = 0:1
)
# grid_points <- df_trt[, c("age", "sex", "kps", "mgmt", "isGTR")]
grid_points <- data.frame(grid_points)
time_cutoff <- 24
# Resample data to replicate study
library(parallel)

# B <- 100
# results <- vector("list", B)
#
# for (b in seq_len(B)) {
#
#   trt_synth <- resample_data(df_trt, nrow(df_trt))
#   con_synth <- resample_data(df_con, nrow(df_con))
#
#   summary_synth <- prepare_summary(trt_synth, con_synth, time_cutoff)
#   target_synth <- summary_synth$target_con
#   control_rmst_synth <- summary_synth$control_rmst_value
#
#   theta_hat_result <- rmst_diff_atc_manual(
#                                            trt_synth,
#                                            con_synth,
#                                            tau = time_cutoff,
#                                            n_boot = 100,
#                                            conf_level = NA # Do not run bootstrap
#   )
#
#   theta_IPE <- estimate_IPE(
#                             df = trt_synth,
#                             target_con = target_synth,
#                             control_rmst_value = control_rmst_synth,
#                             tau = time_cutoff,
#                             grid_points = grid_points,
#                             L = 0.8
#   )
#
#   results[[b]] <- list(
#                        theta_hat_est = theta_hat_result$estimate,
#                        theta_hat_lo  = theta_hat_result$ci[1],
#                        theta_hat_hi  = theta_hat_result$ci[2],
#                        theta_min     = theta_IPE$theta_min,
#                        theta_max     = theta_IPE$theta_max
#   )
#   print(b)
# }

library(future)
library(future.apply)

B <- 96

# Use 20 workers (leave 1-2 free if you want your machine responsive)
plan(multisession, workers = 16)

# Ensure reproducible RNG across parallel workers
# (important since you're resampling/bootstrapping)
results <- future_lapply(seq_len(B), function(b) {

                           trt_synth <- resample_data(df_trt, nrow(df_trt))
                           con_synth <- resample_data(df_con, nrow(df_con))

                           summary_synth <- prepare_summary(trt_synth, con_synth, time_cutoff)
                           target_synth <- summary_synth$target_con
                           control_rmst_synth <- summary_synth$control_rmst_value

                           theta_hat_result <- rmst_diff_atc_manual(
                                                                    trt_synth,
                                                                    con_synth,
                                                                    tau = time_cutoff,
                                                                    n_boot = 400,
                                                                    conf_level = 0.95
                           )

                           theta_IPE <- estimate_IPE(
                                                     df = trt_synth,
                                                     target_con = target_synth,
                                                     control_rmst_value = control_rmst_synth,
                                                     tau = time_cutoff,
                                                     grid_points = grid_points,
                                                     L = 0.8
                           )

                           # Optional: progress message (can be noisy due to parallel output)
                           message(b)

                           list(
                                theta_hat_est = theta_hat_result$estimate,
                                theta_hat_lo  = theta_hat_result$ci[1],
                                theta_hat_hi  = theta_hat_result$ci[2],
                                theta_min     = theta_IPE$theta_min,
                                theta_max     = theta_IPE$theta_max
                           )

  }, future.seed = TRUE)

plan(sequential)
# n_cores <- detectCores() - 4   # leave four cores free

# results <- mclapply(seq_len(B),
#                     function(b) {
#
#                       trt_synth <- resample_data(df_trt, nrow(df_trt))
#                       con_synth <- resample_data(df_con, nrow(df_con))
#
#                       summary_synth <- prepare_summary(trt_synth, con_synth, time_cutoff)
#                       target_synth <- summary_synth$target_con
#                       control_rmst_synth <- summary_synth$control_rmst_value
#
#                       theta_hat_result <- rmst_diff_atc_manual(trt_synth, con_synth, tau = time_cutoff, n_boot = 100
#                       )
#
#                       theta_IPE <- estimate_IPE(
#                                                 df = trt_synth,
#                                                 target_con = target_synth,
#                                                 control_rmst_value = control_rmst_synth,
#                                                 tau = time_cutoff,
#                                                 grid_points = grid_points,
#                                                 L = 0.8
#                       )
#
#                       list(
#                            theta_hat_est = theta_hat_result$estimate,
#                            theta_hat_lo  = theta_hat_result$ci[1],
#                            theta_hat_hi  = theta_hat_result$ci[2],
#                            theta_min     = theta_IPE$theta_min,
#                            theta_max     = theta_IPE$theta_max
#                       )
#
#                     }, mc.cores = n_cores)


library(ggplot2)
library(dplyr)
library(tidyr)

# Extract values from list
df_plot <- tibble(
                  theta_hat = sapply(results, `[[`, "theta_hat_est"),
                  theta_min = sapply(results, `[[`, "theta_min"),
                  theta_max = sapply(results, `[[`, "theta_max")
)

# Convert to long format
df_long <- df_plot %>%
  pivot_longer(cols = everything(),
               names_to = "parameter",
               values_to = "value")

# Single combined boxplot
ggplot(df_long, aes(x = parameter, y = value)) +
  geom_boxplot() +
  labs(
       x = NULL,
       y = "Estimate",
       title = "Empirical Distributions of Theta Estimates"
       ) +
  theme_minimal(base_size = 30)

# Compute summary statistics
df_summary <- df_long %>%
  group_by(parameter) %>%
  summarise(
    q2.5  = quantile(value, 0.025, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    q97.5 = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
ggplot(df_summary, aes(x = parameter, y = median)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.2) +
  labs(
    x = NULL,
    y = "Estimate",
    title = "Median and 95% Quantile Intervals"
  ) +
  theme_minimal(base_size = 30)

ggplot(df_long, aes(x = value, fill = parameter)) +
  geom_histogram(alpha = 0.5, position = "identity", bins = 30) +
  labs(
    x = "Value",
    y = "Count",
    title = "Distributions"
  ) +
  theme_minimal(base_size = 30)
