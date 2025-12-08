library(dplyr)
set.seed(101)
# Population
sample_space <- 1:3
probabilities <- c(1/3, 1/3, 1/3)
noise_y <- 0
regression_func <- function(x) {
  ifelse(x == 1, 0, -1)
}

# Population summary
population <- data.frame(X = sample_space,
                         proportion = probabilities,
                         mean_y = regression_func(sample_space)
)
# A function to solve LP based on summarized data
solve_lp <- function(categorized) {
  prop <- categorized$proportion # gives p1, p2, p3
  mean_y <- categorized$mean_y # gives y1, y2, y3
  # Problem: Minimize w1p1y1 + w2p2y2 + w3p3y3 s.t. w1p1 + w2p2 + w3p3 = 1, wi >= 0
  optimal_index <- which.min(mean_y)
  wts_min <- prop * 0
  wts_min[optimal_index] <- 1 / prop[optimal_index]
  theta_min <- min(mean_y)
  list(
       theta_min = theta_min,
       wts_min = wts_min,
       optimal_index = optimal_index
  )
}
solve_lp <- function(categorized1, categorized2, ratio_bound = 0.5) {
  # categorized1 gives objective, categorized2 gives constraint
  # Needs lpSolve
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required but not installed.")
  }

  prop   <- categorized1$proportion  # p_i
  mean_y <- categorized1$mean_y      # y_i
  prop2  <- categorized2$proportion
  n      <- length(prop)
  stopifnot(n == length(prop2))

  if (length(mean_y) != n) {
    stop("proportion and mean_y must have the same length.")
  }

  # Objective: minimize sum_i w_i * p_i * y_i
  obj <- prop * mean_y

  # Constraint 1: sum_i w_i * p_i = 1
  A   <- matrix(prop2, nrow = 1)
  dir <- "="
  rhs <- 1

  # for (i in seq_len(n)) {
  #   row <- numeric(n)
  #   row[i] <- 1
  #   A <- rbind(A, row)
  #   dir <- c(dir, "<=")
  #   rhs <- c(rhs, ratio_bound)
  # }
  # Constraint 2: each w_i p_i >= ratio_bound min{p_i}
  for (i in seq_len(n)) {
    row <- numeric(n)
    row[i] <- prop2[i]
    A <- rbind(A, row)
    dir <- c(dir, ">=")
    rhs <- c(rhs, min(prop2) * ratio_bound)
  }
  # # Constraint 2: ratio constraints w_i <= ratio_bound * w_j for all i, j
  # # Implemented as: w_i - ratio_bound * w_j <= 0
  # for (i in seq_len(n)) {
  #   for (j in seq_len(n)) {
  #     if (i == j) next
  #     row       <- numeric(n)
  #     # row[i]    <- prop[i]
  #     row[i]    <- 1
  #     # row[j]    <- -ratio_bound * prop[j]
  #     row[j]    <- -ratio_bound
  #     A         <- rbind(A, row)
  #     dir       <- c(dir, ">=")
  #     rhs       <- c(rhs, 0)
  #   }
  # }


  # Solve LP
  sol <- lpSolve::lp(
                     direction   = "min",
                     objective.in = obj,
                     const.mat   = A,
                     const.dir   = dir,
                     const.rhs   = rhs
  )

  if (sol$status != 0) {
    stop("LP did not find an optimal solution. Status code: ", sol$status)
  }

  wts_min <- sol$solution
  theta_min <- sum(obj * wts_min)  # optimal objective value

  # "Optimal index": here I return the index with the largest weight
  optimal_index <- which.max(wts_min)

  list(
       theta_min     = theta_min,
       wts_min       = wts_min,
       optimal_index = optimal_index
  )
}

theta_min_true <- solve_lp(population, population)$theta_min

results_list <- list()
k <- 0
for (n1 in c(200, 400, 800, 1600)) {
  B <- 800
  theta_min_hat <- rep(NA, B)
  theta_min_hat_cv <- rep(NA, B)
  theta_min_hat_cv2 <- rep(NA, B)
  optimal_index_hat <- rep(NA, B)
  wts_min_hat <- matrix(NA, nrow = B, ncol = length(sample_space))
  for (b in seq_len(B)) {
    trt <- data.frame(matrix(NA, nrow = n1, ncol = 2))
    colnames(trt) <- c("X", "Y")
    trt$X <- sample(sample_space, size = n1, replace = TRUE, prob = probabilities)
    trt$Y <- regression_func(trt$X)
    trt$Y <- trt$Y + rnorm(n1, mean = 0, sd = noise_y)
    # Generate split B
    trt_B <- data.frame(matrix(NA, nrow = n1, ncol = 2))
    colnames(trt_B) <- c("X", "Y")
    trt_B$X <- sample(sample_space, size = n1, replace = TRUE, prob = probabilities)
    trt_B$Y <- regression_func(trt_B$X)
    trt_B$Y <- trt_B$Y + rnorm(n1, mean = 0, sd = noise_y)


    # Summarize
    summ_trt <- trt %>%
      group_by(X) %>%
      summarise(mean_y = mean(Y),
                proportion = n() / n1
      )
      summ_trt_B <- trt_B %>%
        group_by(X) %>%
        summarise(mean_y = mean(Y),
                  proportion = n() / n1
        )
      # Create hybrid population
        summ_trt_h <- (1 - 0.00001) * summ_trt + 0.00001 * summ_trt_B
      # Solve linear program
      result <- solve_lp(summ_trt, summ_trt)
      theta_min_hat[b] <- result$theta_min
      optimal_index_hat[b] <- result$optimal_index
      wts_min_hat[b, ] <- result$wts_min
      # Plug-in weights on data B
      theta_min_hat_cv[b] <- sum(result$wts_min * summ_trt_B$proportion * summ_trt_B$mean_y)
      # Solve on hybrid population
      result_h <- solve_lp(summ_trt, summ_trt_h)
      tmp <- sum(result_h$wts_min * summ_trt$proportion * summ_trt$mean_y)
      theta_min_hat_cv2[b] <- theta_min_hat_cv[b] + (tmp - theta_min_hat[b]) * (1/0.00001)
  }
  # Put weights + optimal index into a data frame
  wts_df <- data.frame(theta_min = theta_min_hat,
                       optimal_index = optimal_index_hat,
                       theta_min_cv = theta_min_hat_cv,
                       theta_min_cv2 = theta_min_hat_cv2
  )
  for (i in seq_along(sample_space)) {
    wts_df[, paste0("w", i)] <- wts_min_hat[, i]
  }
  wts_df$n1 <- n1
  k <- k + 1
  results_list[[k]] <- wts_df
}
results <- bind_rows(results_list)
results <- results %>% mutate(scaled_error = sqrt(n1) * (theta_min - theta_min_true))
results <- results %>% mutate(scaled_error_cv = sqrt(n1) * (theta_min_cv - theta_min_true))
results <- results %>% mutate(scaled_error_cv2 = sqrt(n1) * (theta_min_cv2 - theta_min_true))
# Visualize results
library(ggplot2)
library(tidyr)
source("myggplot_theme.R", chdir = TRUE)
ggplot(results) +
  geom_boxplot(aes(x = factor(n1), y = scaled_error)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme_mydefault(base_size = 20)

ggplot(results) +
  geom_boxplot(aes(x = factor(n1), y = scaled_error_cv)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme_mydefault(base_size = 20)

ggplot(results) +
  geom_boxplot(aes(x = factor(n1), y = scaled_error_cv2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  theme_mydefault(base_size = 20)


results %>% group_by(n1) %>%
  summarise(mean(scaled_error), mean(scaled_error_cv), mean(scaled_error_cv2)) %>%
  print()
# See the optimal weights
ggplot(results %>% filter(optimal_index %in% 2:3) %>% mutate(weight = pmax(w1, w2, w3))) +
  geom_boxplot(aes(x = factor(n1), y = sqrt(n1) * (weight - 3))) +
  facet_wrap(~ optimal_index) +
  theme_mydefault(base_size = 20)
# See the optimal weights
ggplot(results %>% filter(optimal_index %in% 2:3)) +
  geom_point(aes(x = w2, y = w3, color = factor(optimal_index)), size = 0.5) +
  facet_wrap(~ factor(n1)) +
  theme_mydefault(base_size = 20)

