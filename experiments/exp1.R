# This experiment looks at the behavior of theta_min under
# Scenario 1-3
# Regularization: Lipschitz, Gaussian mixture, RKHS ball, exp RKHS
# Correction: No correction

rm(list = ls())
scenario <- 1
source(paste0("./generate/scenario", scenario, ".R"), chdir = TRUE, local = TRUE)
source("./methods/weighted_gaussian.R", chdir = TRUE, local = TRUE)
source("./methods/lipschitz.R", chdir = TRUE, local = TRUE)
source("./methods/RKHS.R", chdir = TRUE, local = TRUE)
source("./methods/RKHS_log.R", chdir = TRUE, local = TRUE)

library(parallel)
library(ggplot2)
library(latex2exp)
library(dplyr)
set.seed(100)
L <- list(
          L_lipschitz = 8,
          L_gaussian = 8,
          L_rkhs = 16,
          L_rkhs_log = 16
)

# Create population
population <- get_population(scenario = scenario)
with(population$categorized, plot(X, prop))
with(population$categorized, plot(X, mean_Y))

# Find population minimum
theta_true <- list()
# Lipschitz
population_minimum_lipschitz <- true_theta_min_lipschitz(population, L = L$L_lipschitz)
theta_true$lipschitz <- population_minimum_lipschitz$theta_min
summarize_theta_lipschitz(population_minimum_lipschitz)
# Gaussian mixture
gaussian_grid_points <- ppoints(50)
population_minimum_gaussian <- true_theta_min_gaussian(population,
                                                       grid_points = gaussian_grid_points,
                                                       L = L$L_gaussian)
theta_true$gaussian <- population_minimum_gaussian$theta_min
summarize_theta_gaussian(population_minimum_gaussian)
# Gaussian RKHS norm
rkhs_sigma <- 1 / 8
population_minimum_rkhs <- true_theta_min_rkhs(population,
                                               L = L$L_rkhs,
                                               sigm = rkhs_sigma)
theta_true$rkhs <- population_minimum_rkhs$theta_min
summarize_theta_rkhs(population_minimum_rkhs)
# Gaussian RKHS log norm
rkhs_log_sigma <- 1 / 8
population_minimum_rkhs_log <- true_theta_min_rkhs_log(population,
                                                       L = L$L_rkhs_log,
                                                       sigm = rkhs_sigma)
theta_true$rkhs_log <- population_minimum_rkhs_log$theta_min
summarize_theta_rkhs_log(population_minimum_rkhs_log)

# stopifnot(FALSE)
# Estimate on data
B <- 800
n1_values <- c(200, 400, 800, 1600)
n0_values <- 100 * n1_values

# Helper to run one replication for a given (n0, n1)
run_one_rep <- function(n0, n1, gaussian_grid_points, L) {
  dat <- get_data(n0 = n0, n1 = n1, scenario = scenario)
  s_min_lipschitz <- estimate_theta_min_lipschitz(dat, L = L$L_lipschitz)
  s_min_gaussian <- estimate_theta_min_gaussian(dat, gaussian_grid_points, L = L$L_gaussian)
  s_min_rkhs <- estimate_theta_min_rkhs(dat, L = L$L_rkhs, sigm = rkhs_sigma)
  s_min_rkhs_log <- estimate_theta_min_rkhs_log(dat, L = L$L_rkhs_log, sigm = rkhs_log_sigma)

  c(
    lipschitz = s_min_lipschitz$theta_min,
    gaussian = s_min_gaussian$theta_min,
    rkhs = s_min_rkhs$theta_min,
    rkhs_log = s_min_rkhs_log$theta_min
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
                                     gaussian_grid_points = gaussian_grid_points,
                                     L = L),
                       mc.cores = ncores
  )
  est_mat <- do.call(cbind, est_list)

  ## Convert to long data.frame
  df_k <- as.data.frame(t(est_mat))
  df_k$rep <- seq_len(B)
  df_k$n1 <- n1
  df_k$n0 <- n0

  ## Long/tidy format
  df_k_long <- reshape(df_k,
                       varying = list(names(df_k)[names(df_k) %in% c("lipschitz", "gaussian", "rkhs", "rkhs_log")]),
                       v.names = "theta_hat",
                       timevar = "method",
                       times   = c("lipschitz", "gaussian", "rkhs", "rkhs_log"),
                       direction = "long"
  )
  rownames(df_k_long) <- NULL

  ## Attach the correct method-specific true value
  theta_map <- c(
                 lipschitz = theta_true$lipschitz,
                 gaussian  = theta_true$gaussian,
                 rkhs      = theta_true$rkhs,
                 rkhs_log  = theta_true$rkhs_log
  )

  df_k_long$theta_true <- theta_map[df_k_long$method]
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
