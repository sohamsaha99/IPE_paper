# This experiment looks at the behavior of theta_min under
# Scenario 1: Quadratic regression function,
# Regularization: Gaussian mixture,
# Correction: Different data-splitting approaches

rm(list = ls())
source("./generate/scenario1.R", chdir = TRUE, local = TRUE)
source("./methods/weighted_gaussian.R", chdir = TRUE, local = TRUE)
library(parallel)
library(ggplot2)
library(latex2exp)
library(dplyr)
set.seed(100)
L <- 8

# Create population
population <- get_population_scenario1()
with(population$categorized, plot(X, prop))
with(population$categorized, plot(X, mean_Y))

# Find population minimum
gaussian_grid_points <- ppoints(50)
population_minimum <- true_theta_min_gaussian(population,
                                              grid_points = gaussian_grid_points,
                                              L = L)
theta_true <- population_minimum$theta_min
summarize_theta_gaussian(population_minimum)

# Estimate on data
B <- 800
n1_values <- c(400, 800, 1600)
n0_values <- 100 * n1_values

# Helper to run one replication for a given (n0, n1)
run_one_rep <- function(n0, n1, grid_points, L) {
  dat <- get_data_scenario1(n0 = n0, n1 = n1)
  s_min <- estimate_theta_min_gaussian(dat, grid_points, L = L)
  cva <- splitting_basic_gaussian(dat, grid_points, L = L)
  cvb <- splitting_lagrange_gaussian(dat, grid_points, L = L)
  cvc <- splitting_derivative_gaussian(dat, grid_points, L = L)
  cvd <- splitting_inverse_gaussian(dat, grid_points, L = L)

  c(
    empirical = s_min$theta_min,
    split_basic = cva,
    split_lagrange = cvb,
    split_derivative = cvc,
    split_inverse = cvd
  )
}

# Run simulations across sample sizes
all_results <- list()
ncores <- max(1, detectCores() - 1)
for (k in seq_along(n1_values)) {
  n1 <- n1_values[k]
  n0 <- n0_values[k]

  est_list <- mclapply(X = seq_len(B),
                       FUN = function(b)
                         run_one_rep(n0 = n0, n1 = n1,
                                     grid_points = gaussian_grid_points,
                                     L = L),
                       mc.cores = ncores
  )
  est_mat <- do.call(cbind, est_list)

  ## Convert to long data.frame
  df_k <- as.data.frame(t(est_mat))
  df_k$rep <- seq_len(B)
  df_k$n1 <- n1
  df_k$n0 <- n0
  df_k$theta_true <- theta_true

  ## Long/tidy format
  df_k_long <- reshape(df_k,
                       varying = list(names(df_k)[names(df_k) %in% c("empirical","split_basic","split_lagrange", "split_derivative", "split_inverse")]),
                       v.names = "theta_hat",
                       timevar = "method",
                       times   = c("empirical", "split_basic", "split_lagrange", "split_derivative", "split_inverse"),
                       direction = "long"
  )
  rownames(df_k_long) <- NULL

  ## Scaled error: sqrt(n1)*(hat - true)
  df_k_long$scaled_error <- sqrt(df_k_long$n1) * (df_k_long$theta_hat - df_k_long$theta_true)

  all_results[[k]] <- df_k_long
  print(c(n0, n1))
}

results <- do.call(rbind, all_results)

## Peek
head(results)

# Boxplots of sqrt(n1)*(hat - true)
library(ggplot2)
source("./myggplot_theme.R", chdir = TRUE)

ggplot(results, aes(x = method, y = scaled_error)) +
  geom_boxplot() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ n1) +
  labs(
    title = expression(paste("Boxplots of ", sqrt(n[1]), "(", hat(theta), " - ", theta^"*", ")")),
    x = "Method", y = "Scaled error"
  ) +
  ylim(c(-3, 3)) +
  theme_minimal() +
  theme_mydefault(base_size = 20)

# Density plots of sqrt(n1)*(hat - true)
ggplot(results, aes(x = scaled_error, fill = method)) +
  geom_density(aes(y = after_stat(density)),
                 alpha = 0.5, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ n1) +
  labs(
       title = expression(paste("Densities of ", sqrt(n[1]), "(", hat(theta), " - ", theta^"*", ")")),
       x = "Scaled error", y = "Density", fill = "Method"
       ) +
  xlim(c(-4, 4)) +
  theme_minimal() +
  theme_mydefault(base_size = 20)

# Summary of sqrt(n1)*(hat - true)
summary_scaled <- results %>%
  group_by(n0, n1, method) %>%
  summarise(bias  = mean(scaled_error, na.rm = TRUE),
            sd    = sd(scaled_error, na.rm = TRUE),
            rmse  = sqrt(mean(scaled_error^2, na.rm = TRUE)),
            failed_prop = mean(is.na(scaled_error)),
            .groups = "drop"
  )
print(summary_scaled)


# Visualize relationship between methods
library(tidyr)
library(GGally)

results_wide <- results %>%
  select(rep, n1, method, scaled_error) %>%
  pivot_wider(names_from = method, values_from = scaled_error)

ggpairs(results_wide,
        columns = which(!names(results_wide) %in% c("rep", "n1")),
        aes(alpha = 0.6, color = factor(n1))
        ) +
theme_minimal() +
theme_mydefault(base_size = 20)
