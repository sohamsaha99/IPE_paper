rm(list = ls())

source("./experiments/exp3_v2/estimate_with_PLD.R")
source("./experiments/exp3_v2/estimate_with_summary.R")
source("./myggplot_theme.R")

# # Read data
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
# For now we only use DFCI data to make trt and control
# set.seed(20260211)
# df <- read.csv("./experiments/exp3_v2/GBM_files/GBM_con.csv", header = TRUE, stringsAsFactors = FALSE)
# df <- df %>%
#   mutate(isGTR = as.numeric(eor == "GTR")) %>%
#   select(age, sex, kps, mgmt, isGTR, os, os_status)
# idx_t <- sample(nrow(df), round(nrow(df)) / 2, replace = FALSE)
# df_trt <- df[idx_t, ]
# df_con <- df[-idx_t, ]
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
### --- Summary data prepared ---

# Estimate treatment effect based on PLD control data
# theta_hat <- rmst_diff_atc_adjcurves(df_trt, df_con, tau = 24, n_boot = 500)
theta_hat_2 <- rmst_diff_atc_manual(df_trt, df_con, tau = 24, n_boot = 500)
# print(theta_hat)
print(theta_hat_2[1:3])

# Estimate IPE based on summary data
grid_points <- expand.grid(age = seq(10, 100, by = 1),
                           sex = 0:1,
                           kps = 0:1,
                           mgmt = 0:1,
                           isGTR = 0:1
)
# grid_points <- df_trt[, c("age", "sex", "kps", "mgmt", "isGTR")]
grid_points <- data.frame(grid_points)
theta_IPE <- estimate_IPE(df = df_trt,
                          target_con = target,
                          control_rmst_value = control_rmst_value,
                          tau = time_cutoff,
                          grid_points = grid_points,
                          L = 0.8
)
# Print the IPE
print(c(theta_IPE$theta_min, theta_IPE$theta_max))

# Visualize the weights
wts_df <- grid_points
wts_df$wts_min <- NA
for (i in seq_len(nrow(wts_df))) {
  wts_df$wts_min[i] <- density_ratio_wf(x = as.numeric(grid_points[i, ]),
                                        alpha = theta_IPE$alpha_min,
                                        grid_points = grid_points,
                                        L = theta_IPE$L,
                                        Sigma = theta_IPE$Sigma
  )
}
library(ggplot2)
library(ggpubr)
prepare_for_plot <- function(df) {
  df %>% mutate(sex = factor(sex, levels = 0:1, labels = c("Female", "Male")),
                kps = factor(kps, levels = 0:1, labels = c("KPS low", "KPS high")),
                mgmt = factor(mgmt, levels = 0:1, labels = c("Unmethylated", "Methylated")),
                eor = factor(isGTR, levels = 0:1, labels = c("STR", "GTR"))
  )
}

ggplot(wts_df %>% prepare_for_plot()) +
  geom_line(aes(x = age, y = wts_min, group = sex, color = sex)) +
  facet_grid(
             kps ~ mgmt + eor,
             labeller = label_both
             ) +
  theme_bw(base_size = 30)

ggplot(df_trt %>% prepare_for_plot() %>%
       mutate(wts_min = theta_IPE$wts_min)) +
  geom_point(aes(x = age, y = wts_min, group = sex, color = sex)) +
  facet_grid(
             kps ~ mgmt + eor,
             labeller = label_both
             ) +
  theme_bw(base_size = 30)

# Compare distributions f1 and f0
ggplot(df_con %>% prepare_for_plot()) +
  geom_histogram(aes(x = age, fill = "f0"),
                 bins = 7, alpha = 0.6, binwidth = 10, boundary = 15
  ) +
  geom_histogram(
                 data = df_trt %>% prepare_for_plot(),
                 aes(x = age, fill = "f1"),
                 bins = 7, alpha = 0.6, binwidth = 10, boundary = 15
  ) +
  facet_grid(sex + kps ~ mgmt + eor, labeller = label_both) +
  scale_fill_manual(
                    name = "",
                    values = c("f0" = "forestgreen", "f1" = "orange")
  ) +
  theme_minimal(base_size = 30)

# Compare distributions f1 and f0,min
ggplot(df_trt %>% prepare_for_plot()) +
  geom_histogram(aes(x = age, fill = "f1"),
                 bins = 7, alpha = 0.6, binwidth = 10, boundary = 15
                 ) +
  geom_histogram(data = df_trt %>% prepare_for_plot() %>% mutate(wts_min = theta_IPE$wts_min),
                 aes(x = age, weight = wts_min, fill = "f0min"),
                 bins = 7, alpha = 0.6, binwidth = 10, boundary = 15
  ) +
  facet_grid(
             sex + kps ~ mgmt + eor,
             labeller = label_both
             ) +
  scale_fill_manual(
                    name = "",
                    values = c("f1" = "orange", "f0min" = "steelblue")
  ) +
  theme_minimal(base_size = 30)


# Compare distributions f0 and f0,min
ggplot(df_con %>% prepare_for_plot()) +
  geom_histogram(aes(x = age, fill = "f0"),
                 bins = 7, alpha = 0.6, binwidth = 10, boundary = 15
                 ) +
  geom_histogram(data = df_trt %>% prepare_for_plot() %>% mutate(wts_min = theta_IPE$wts_min),
                 aes(x = age, weight = wts_min, fill = "f0min"),
                 bins = 7, alpha = 0.6, binwidth = 10, boundary = 15
  ) +
  facet_grid(
             sex + kps ~ mgmt + eor,
             labeller = label_both
             ) +
  scale_fill_manual(
                    name = "",
                    values = c("f0" = "forestgreen", "f0min" = "steelblue")
  ) +
  theme_minimal(base_size = 30)


