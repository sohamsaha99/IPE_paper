rm(list = ls())

source("./experiments/exp3_v2/estimate_with_PLD.R")
source("./experiments/exp3_v2/estimate_with_summary.R")

# # Read data
# df_trt <- read.csv("./experiments/exp3_v2/GBM_files/GBM_trt.csv", header = TRUE, stringsAsFactors = FALSE)
# df_con <- read.csv("./experiments/exp3_v2/GBM_files/GBM_con.csv", header = TRUE, stringsAsFactors = FALSE)
#
# # Keep required columns only
# print(dim(df_trt))
# print(dim(df_con))
# df_trt <- df_trt %>%
#   mutate(isGTR = as.numeric(eor == "GTR")) %>%
#   select(age, sex, kps, mgmt, isGTR, os, os_status)
# df_con <- df_con %>%
#   mutate(isGTR = as.numeric(eor == "GTR")) %>%
#   select(age, sex, kps, mgmt, isGTR, os, os_status)
#

#####
# For now we only use DFCI data to make trt and control
df <- read.csv("./experiments/exp3_v2/GBM_files/GBM_con.csv", header = TRUE, stringsAsFactors = FALSE)
df <- df %>%
  mutate(isGTR = as.numeric(eor == "GTR")) %>%
  select(age, sex, kps, mgmt, isGTR, os, os_status)
idx_t <- sample(nrow(df), round(nrow(df)) / 2, replace = FALSE)
df_trt <- df[idx_t, ]
df_con <- df[-idx_t, ]
### --- Individual Data prepared --- ###

# Create Summary statistics
target <- data.frame(xname = c("age", "age", "sex", "kps", "mgmt", "isGTR"),
                     type  = c(1, 2, 1, 1, 1, 1)
)
target$value <- c(
                  mean(df_con$age),
                  mean(df_con$age^2),
                  mean(df_con$sex),
                  mean(df_con$kps),
                  mean(df_con$mgmt),
                  mean(df_con$isGTR)
)

time_cutoff <- 24
control_rmst_value <- weighted_rmst(df_con$os, df_con$os_status, time_cutoff)$rmst[1]
# Estimate treatment effect based on PLD control data
theta_hat <- rmst_diff_atc_adjcurves(df_trt, df_con, tau = 24, n_boot = 500)

# Estimate IPE based on summary data
theta_IPE <- estimate_IPE(df = df_trt,
                          target_con = target_con,
                          control_rmst_value = control_rmst_value,
)
