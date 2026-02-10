rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(purrr)
library(mvtnorm)
source("myggplot_theme.R")
source("./experiments/exp3/weighted_rmst.R")
source("./experiments/exp3/impute_censored.R")
source("./experiments/exp3/estimate_with_PLD.R")

FILE_PATH <- "./experiments/exp3/GBM_files/GBM_complete.csv"
df_GBM <- read.csv(FILE_PATH, header = TRUE, stringsAsFactors = FALSE)

ggplot(df_GBM %>%
         mutate(sex = factor(sex),
                kps = factor(kps),
                mgmt = factor(mgmt),
                eor = factor(eor))) +
  geom_point(aes(x = age, y = os, color = sex)) +
  facet_grid(eor ~ kps,
             labeller = labeller(.rows = label_both, .cols = label_both)) +
  theme_pubr(border = TRUE) +
  theme_mydefault(base_size = 20)

# Rescale all predictor columns so that they have mean 0 variance 1
df_GBM_processed <- df_GBM %>% mutate(isGTR = as.numeric(eor == "GTR")) %>%
  select(age, sex, kps, mgmt, isGTR, os, os_status) %>%
  mutate(across(c(age, sex, kps, mgmt, isGTR),
                ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)))

# Check normalization
df_GBM_processed %>%
  summarise(across(
    everything(),
    list(
      mean = ~ mean(., na.rm = TRUE),
      var  = ~ var(., na.rm = TRUE)
    )
  ))

# Split in treatment and control groups
set.seed(20260201)
idx_trt <- sample(x = nrow(df_GBM_processed),
                  size = floor(0.5 * nrow(df_GBM_processed)),
                  replace = FALSE)
df_trt <- df_GBM_processed[idx_trt, ]
df_con <- df_GBM_processed[-idx_trt, ]

# Cutoff
tau <- 24

# Create control population summary
target_con <- data.frame(matrix(NA, nrow = 6, ncol = 3))
colnames(target_con) <- c("xname", "value", "type")
target_con$xname     <- c("age", "age", "sex", "kps", "mgmt", "isGTR")
target_con$type      <- c("first", "second", rep("first", 4))
target_con$value <- c(mean(df_con$age),
                      mean(df_con$age^2),
                      mean(df_con$sex),
                      mean(df_con$kps),
                      mean(df_con$mgmt),
                      mean(df_con$isGTR)
)
control_rmst <- weighted_rmst(time_var = df_con$os,
                              event_var = df_con$os_status,
                              time_cutoff = tau
)
control_rmst_alt <- simtrial:::rmst_single_arm(time_var = df_con$os,
                                               event_var = df_con$os_status,
                                               tau = tau
)

# In treatnent group, replace missing event time with E(min{T, tau}|X, Y, delta)
df_trt_imputed <- impute_censored(df = df_trt, method = "emin_tau", tau = tau)

grid_points <- expand.grid(age   = seq(min(df_GBM_processed$age), max(df_GBM_processed$age), length.out = 50),
                           sex   = unique(df_GBM_processed$sex),
                           kps   = unique(df_GBM_processed$kps),
                           mgmt  = unique(df_GBM_processed$mgmt),
                           isGTR = unique(df_GBM_processed$isGTR)
)
grid_points <- data.frame(grid_points)

# Function for estimating true theta_min
estimate_theta_min <- function(df_trt_imputed,
                               target_con,
                               control_rmst_value,
                               grid_points,
                               L) {
  # Compute dmvnorm
  normal_den <- function(x, v) {
    stopifnot(ncol(x) == length(v))
    dmvnorm(x, mean = v, sigma = diag(length(v)) / (L^2))
  }
  X_cols <- c("age", "sex", "kps", "mgmt", "isGTR")
  # Y column contains E(min{T, tau}|X, Y, delta)
  df_trt_imputed$Y <- df_trt_imputed$emin_tau
  df_trt_imputed$prob <- 1 / nrow(df_trt_imputed)
  df_trt <- df_trt_imputed
  # Set linear programming:
  # min c^\top alpha
  # subject to A alpha = target_con$value
  K <- nrow(grid_points)
  cvec <- rep(NA, K)
  for (k in seq_len(K)) {
    cvec[k] <- sum(df_trt$Y * df_trt$prob *
                   normal_den(df_trt[, X_cols], as.numeric(grid_points[k, ])))
  }
  rhsvec <- c(1, target_con$value)
  dirvec <- rep("==", length(rhsvec))
  Amat <- matrix(NA, nrow = length(rhsvec), ncol = K)
  # Set up each row
  for (k in seq_len(K)) {
    den_k <- normal_den(df_trt[, X_cols], as.numeric(grid_points[k, ]))
    stopifnot(length(den_k) == nrow(df_trt))
    Amat[1, k] <- sum(df_trt$prob  * den_k)
    Amat[2, k] <- sum(df_trt$age   * df_trt$prob * den_k)
    Amat[3, k] <- sum(df_trt$age^2 * df_trt$prob * den_k)
    Amat[4, k] <- sum(df_trt$sex   * df_trt$prob * den_k)
    Amat[5, k] <- sum(df_trt$kps   * df_trt$prob * den_k)
    Amat[6, k] <- sum(df_trt$mgmt  * df_trt$prob * den_k)
    Amat[7, k] <- sum(df_trt$isGTR * df_trt$prob * den_k)
  }

  # Solve linear program
  # Solve LP: min c^T alpha  s.t. Amat %*% alpha = rhsvec, alpha >= 0
  lp_out <- lpSolve::lp(direction    = "min",
                        objective.in = cvec,
                        const.mat    = Amat,
                        const.dir    = dirvec,
                        const.rhs    = rhsvec
  )

  if (lp_out$status != 0) {
    stop("LP did not solve to optimality. lpSolve status code: ", lp_out$status)
  }

  alpha_hat <- lp_out$solution
  EY1_min <- as.numeric(lp_out$objval)
  theta_min <- EY1_min - control_rmst_value

  list(theta_min = theta_min,
       alpha = alpha_hat,
       L = L
  )
}

estimate_theta_max <- function(df_trt_imputed,
                               target_con,
                               control_rmst_value,
                               grid_points,
                               L) {

  df2 <- df_trt_imputed
  df2$emin_tau <- -1.0 * df2$emin_tau

  out <- estimate_theta_min(df_trt_imputed     = df2,
                            target_con         = target_con,
                            control_rmst_value = -control_rmst_value,
                            grid_points        = grid_points,
                            L                  = L)

  list(theta_max = -out$theta_min,
       alpha     = out$alpha,
       L = L
  )
}



estimate_theta_min(df_trt_imputed, target_con, control_rmst$rmst[1],
                   grid_points = grid_points,
                   L = 1.0
)

# Plot for range of values of L
library(ggplot2)

# Grid of L values
L_seq <- seq(0.6, 2, length.out = 20)

# Compute theta_min and theta_max
theta_min <- sapply(L_seq, function(L)
  estimate_theta_min(df_trt_imputed, target_con, control_rmst$rmst[1],
                     grid_points = grid_points, L = L
                     )$theta_min
)

theta_max <- sapply(L_seq, function(L)
  estimate_theta_max(df_trt_imputed, target_con, control_rmst$rmst[1],
                     grid_points = grid_points, L = L
                     )$theta_max
)

# Data frame for plotting
result <- data.frame(
  L = L_seq,
  theta_min = theta_min,
  theta_max = theta_max
)

# Plot
p1 <- ggplot(result, aes(x = L)) +
  geom_ribbon(aes(ymin = theta_min, ymax = theta_max), alpha = 0.3) +
  geom_line(aes(y = theta_min), linetype = "dashed") +
  geom_line(aes(y = theta_max), linetype = "dashed") +
  labs(
       x = "L",
       y = expression(theta),
       title = "Usual Summary"
       ) +
  ylim(c(-6, 6)) +
  theme_minimal() +
  theme_mydefault(base_size = 20)
print(p1)

# Visualize the optimizer weights
# w_f(x) = sum_k alpha_k * N(x | v_k, I_p / L^2)
density_ratio_wf <- function(x, alpha, grid_points, L) {
  # x: numeric vector in R^p
  # alpha: numeric vector length K
  # grid_points: data.frame / matrix with K rows, p cols (v_1,...,v_K)
  # L: positive scalar

  x <- as.numeric(x)
  V <- as.matrix(grid_points)
  alpha <- as.numeric(alpha)

  p <- length(x)
  K <- nrow(V)

  stopifnot(ncol(V) == p)
  stopifnot(length(alpha) == K)
  stopifnot(is.finite(L), L > 0)

  # Evaluate phi(x; v_k, I_p / L^2) for all k
  Sigma <- diag(p) / (L^2)
  dens <- mvtnorm::dmvnorm(V, mean = x, sigma = Sigma)  # length K

  sum(alpha * dens)
}

L_vals <- c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6)
age_vec <- seq(from = -4, to = 2.5, length.out = 100)

plot_df <- map_dfr(L_vals, function(L) {

  alpha_vec <- estimate_theta_min(
    df_trt_imputed,
    target_con,
    control_rmst$rmst[1],
    grid_points = grid_points,
    L = L
  )$alpha

  density_vec <- map_dbl(age_vec, function(age) {
    density_ratio_wf(c(age, 0, 0, 0, 0), alpha_vec, grid_points, L)
  })

  tibble(
    age = age_vec * sd(df_GBM$age) + mean(df_GBM$age),
    density_ratio = density_vec,
    L = L
  )
})

ggplot(plot_df, aes(x = age, y = density_ratio)) +
  geom_line() +
  facet_wrap(
    ~ L,
    ncol = 3,
    labeller = labeller(L = function(x) paste0("L=", x))
  ) +
  labs(
    x = "Age",
    y = "Density ratio",
    title = "Density ratio by age for different L"
  ) +
  theme_minimal() +
  theme_mydefault(base_size = 20)

###
# Generate data with additive treatment effect and run estimation methods
generate_trt_data <- function(df_trt, theta) {
  df_trt_spike <- df_trt %>%
    mutate(os = os + ifelse(os_status == 1, theta, 0))
  return(df_trt_spike)
}

theta_true_vec <- seq(from = 0, to = 9, by = 0.5)
theta_hat_matching <- rep(NA, length(theta_true_vec))
lower_int_matching <- rep(NA, length(theta_true_vec))
upper_int_matching <- rep(NA, length(theta_true_vec))
theta_min          <- rep(NA, length(theta_true_vec))
theta_max          <- rep(NA, length(theta_true_vec))
k <- 0
for (k in seq_along(theta_true_vec)) {
  df_trt_spike <- generate_trt_data(df_trt, theta_true_vec[k])
  # Estimate with matching
  ret <- estimate_trt_effet_matching(df_trt_spike, df_con, tau)
  theta_hat_matching[k] <- ret$theta_hat
  lower_int_matching[k] <- ret$lower
  upper_int_matching[k] <- ret$upper
  # Estimate IPE
  df_trt_imputed <- impute_censored(df = df_trt_spike, method = "emin_tau", tau = tau)
  theta_min[k] <- estimate_theta_min(df_trt_imputed, target_con, control_rmst$rmst[1],
                                     grid_points = grid_points,
                                     L = 1.2
  )$theta_min
  theta_max[k] <- estimate_theta_max(df_trt_imputed, target_con, control_rmst$rmst[1],
                                     grid_points = grid_points,
                                     L = 1.2
  )$theta_max
}

# Plot
library(ggplot2)

# Put results into a long data frame (one row per interval type per theta_true)
plot_df <- data.frame(
  theta_true = rep(theta_true_vec, times = 2),
  method     = rep(c("Matching CI", "IPE bounds"), each = length(theta_true_vec)),
  ymin       = c(lower_int_matching, theta_min),
  ymax       = c(upper_int_matching, theta_max)
)

# Optional: if you also want to plot the point estimate for matching
plot_df_pts <- data.frame(
  theta_true = theta_true_vec,
  theta_hat  = theta_hat_matching
)

ggplot(plot_df, aes(x = theta_true, ymin = ymin, ymax = ymax, color = method)) +
  geom_linerange(position = position_dodge(width = 0.25), linewidth = 0.9) +
  # Optional: add matching point estimate (uncomment if desired)
  geom_point(data = plot_df_pts, aes(x = theta_true, y = theta_hat),
             inherit.aes = FALSE, shape = 16, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  labs(
    x = "True treatment effect (months)",
    y = "Interval",
    color = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme_mydefault(base_size = 25)

