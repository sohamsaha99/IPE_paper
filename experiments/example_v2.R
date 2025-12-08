source("myggplot_theme.R", chdir = TRUE)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(numDeriv)
set.seed(101)
# Population
sample_space <- 1:3
probabilities <- c(1/3, 1/3, 1/3)
# probabilities <- c(0.4, 0.3, 0.3)
noise_y <- 0
regression_func <- function(x) {
  ifelse(x == 2, 1, 0)
  # ifelse(x == 1, 1 + probabilities[3] / (0.5 * probabilities[1]),
  #        ifelse(x == 2, 1, 0)
  # )
  # ifelse(x == 1, -3, ifelse(x == 2, 0, 3))
}

# Population summary
population <- data.frame(X = sample_space,
                         proportion = probabilities,
                         mean_y = regression_func(sample_space)
)

# Generate one data
get_data <- function(n) {
  trt             <- data.frame(matrix(NA, nrow = n, ncol = 2))
  colnames(trt)   <- c("X", "Y")
  trt$X           <- sample(sample_space, size = n, replace = TRUE, prob = probabilities)
  trt$Y           <- regression_func(trt$X)
  trt$Y           <- trt$Y + rnorm(n, mean = 0, sd = noise_y)
  trt
}

solve_lp <- function(h, categorized1, categorized2) {
  # categorized1 and categorized2 are two sources of summary data
  # h denotes the mixing proportion: (1-h) categorized1 + h categorized2
  # Needs lpSolve
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required but not installed.")
  }

  prop1   <- categorized1$proportion  # p_i
  prop2   <- categorized2$proportion
  mean_y1 <- categorized1$mean_y      # y_i
  mean_y2 <- categorized2$mean_y      # y_i
  N       <- length(prop1)

  if (length(mean_y1) != N || length(prop2) != N || length(mean_y2) != N) {
    stop("Length mismatch. Stop.")
  }
  if (any(categorized1$X != categorized2$X)) {
    stop("Different X values. Stop.")
  }
  prop   <- (1 - h) * prop1 + h * prop2
  mean_y <- (1 - h) * mean_y1 + h * mean_y2
  # Fix gaussian centers and sd
  centers <- 1:3
  sd_gau <- 1
  # Set up objective
  obj <- rep(NA, length(centers))
  for (k in seq_along(centers)) {
    obj[k] <- 0
    for (i in seq_len(N)) {
      obj[k] <- obj[k] + dnorm(categorized1$X[i], centers[k], sd_gau) * prop1[i] * mean_y1[i]
    }
  }

  # Constraint 1: sum_i w_i * p_i = 1
  a <- rep(NA, length(centers))
  for (k in seq_along(centers)) {
    a[k] <- 0
    for (i in seq_len(N)) {
      a[k] <- a[k] + dnorm(categorized1$X[i], centers[k], sd_gau) * prop[i]
    }
  }
  A   <- matrix(a, nrow = 1)
  dir <- "=="
  rhs <- 1

  # Solve LP
  sol <- lpSolve::lp(
                     direction    = "min",
                     objective.in = obj,
                     const.mat    = A,
                     const.dir    = dir,
                     const.rhs    = rhs
  )

  if (sol$status != 0) {
    stop("LP did not find an optimal solution. Status code: ", sol$status)
  }

  alpha_min <- sol$solution
  theta_min <- sum(obj * alpha_min)  # optimal objective value
  ## Compute A in a vectorized way
  A <- outer(categorized1$X, centers,
             function(x, m) dnorm(x, mean = m, sd = sd_gau))
  wts_min <- as.numeric(A %*% alpha_min)

  # "Optimal index": here I return the index with the largest weight
  optimal_index <- which.max(alpha_min)

  list(
       theta_min     = theta_min,
       alpha_min     = alpha_min,
       wts_min       = wts_min,
       optimal_index = optimal_index
  )
}
find_theta <- function(h, categorized1, categorized2) {
  solve_lp(h, categorized1, categorized2)$theta_min
}

library(ggplot2)

# Visualize feasible sets
feasible_set <- function(prop, prop_hat = NULL, res = 0.001) {
  stopifnot(length(prop) == 3, abs(sum(prop) - 1) <= 1e-8)
  if (!is.null(prop_hat)) {
    stopifnot(length(prop_hat) == 3, abs(sum(prop_hat) - 1) <= 1e-8)
  }

  sample_space <- 1:3

  # Grid over q1 and q2
  q1_vals <- seq(0, 1, by = res)
  q2_vals <- seq(0, 1, by = res)
  grid <- expand.grid(q1 = q1_vals, q2 = q2_vals)

  # Compute q3 and simplex constraint
  grid$q3 <- 1 - grid$q1 - grid$q2
  simplex <- grid$q3 >= 0

  centers <- 1:3
  sd_gau <- 1

  ## Compute A in a vectorized way
  A <- outer(sample_space, centers,
             function(x, m) dnorm(x, mean = m, sd = sd_gau))

  ## Precompute the inverse (or use solve(A, ...) with matrix RHS)
  A_inv <- solve(A)

  # Helper to compute feasibility for a given prop
  compute_feasible <- function(p) {
    # Only compute inside simplex to avoid wasted work
    idx <- which(simplex)

    q1 <- grid$q1[idx]
    q2 <- grid$q2[idx]
    q3 <- grid$q3[idx]

    w1 <- q1 / p[1]
    w2 <- q2 / p[2]
    w3 <- q3 / p[3]

    # Right-hand sides as a 3 x N matrix
    W <- rbind(w1, w2, w3)           # rows: i = 1,2,3; cols: grid points

    # Solve A * alpha = W for all columns at once: alpha is 3 x N
    alpha <- A_inv %*% W

    # Feasibility: all alphas >= 0
    feas_alpha <- alpha[1, ] >= 0 & alpha[2, ] >= 0 & alpha[3, ] >= 0

    # Put back into full-length logical vector
    feasible <- rep(FALSE, nrow(grid))
    feasible[idx] <- feas_alpha

    feasible & simplex
  }

  feasible_true <- compute_feasible(prop)

  # Base empty plot
  p <- ggplot()

  # True polytope: blue, semi-transparent
  p <- p + geom_raster(
                       data = subset(grid, feasible_true),
                       aes(x = q1, y = q2),
                       fill = "steelblue",
                       alpha = 0.4
  )

  # Estimated polytope: orange, semi-transparent
  if (!is.null(prop_hat)) {
    feasible_hat <- compute_feasible(prop_hat)
    p <- p + geom_raster(
                         data = subset(grid, feasible_hat),
                         aes(x = q1, y = q2),
                         fill = "orange",
                         alpha = 0.4
    )
  }

  p +
    coord_fixed() +
    xlim(c(0, 1)) +
    ylim(c(0, 1)) +
    labs(
         title = if (is.null(prop_hat)) {
           "Feasible Set in (q1, q2) Space"
         } else {
           "True (blue) vs Estimated (orange) Feasible Sets"
         },
         x = expression(q[1] == w[1] * p[1]),
         y = expression(q[2] == w[2] * p[2])
         ) +
    theme_mydefault(base_size = 20)
}


feasible_set(prop = probabilities,
             prop_hat = get_data(50) %>% group_by(X) %>% summarise(prop = n() / 50) %>% pull(prop)) +
# q1 m1 + q2 m2 = c
geom_abline(intercept = 1 / (3*regression_func(2)), slope = - regression_func(1) / regression_func(2), color = "red")


theta_min_true <- solve_lp(0, population, population)$theta_min

results_list <- list()
k <- 0
for (n1 in c(50, 100, 200, 400)) {
  B <- 800
  theta_min_hat     <- rep(NA, B)
  theta_min_hat_cv  <- rep(NA, B)
  theta_min_hat_cv2 <- rep(NA, B)
  optimal_index_hat <- rep(NA, B)
  wts_min_hat <- matrix(NA, nrow = B, ncol = length(sample_space))
  for (b in seq_len(B)) {
    n1A <- round(n1 / 2)
    n1B <- n1 - n1A
    trt_A <- get_data(n1A)
    # Generate split B
    trt_B <- get_data(n1B)
    # Full data
    trt <- bind_rows(trt_A, trt_B)
    # Summarize
    summ_trt <- trt %>%
      group_by(X) %>%
      summarise(mean_y = mean(Y),
                proportion = n() / n1
      )
    summ_trt_A <- trt_A %>%
      group_by(X) %>%
      summarise(mean_y = mean(Y),
                proportion = n() / n1A
      )
      summ_trt_B <- trt_B %>%
        group_by(X) %>%
        summarise(mean_y = mean(Y),
                  proportion = n() / n1B
        )
      # Solve linear program
      result               <- solve_lp(0, summ_trt, summ_trt)
      theta_min_hat[b]     <- result$theta_min
      optimal_index_hat[b] <- result$optimal_index
      wts_min_hat[b, ]     <- result$wts_min
      # Plug-in weights on data B
      result_A             <- solve_lp(0, summ_trt_A, summ_trt_A)
      theta_min_hat_cv[b]  <- sum(result_A$wts_min * summ_trt_B$proportion * summ_trt_B$mean_y)
      # Find derivative based on hybrid population
      derivative           <- grad(find_theta, 0, method = "Richardson", side = +1,
                                   categorized1 = summ_trt_A, categorized2 = summ_trt_B)
      theta_min_hat_cv2[b] <- theta_min_hat_cv[b] + derivative
  }
  # Put weights + optimal index into a data frame
  wts_df <- data.frame(theta_min     = theta_min_hat,
                       optimal_index = optimal_index_hat,
                       theta_min_cv  = theta_min_hat_cv,
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
results <- results %>% mutate(scaled_error     = sqrt(n1) * (theta_min - theta_min_true))
results <- results %>% mutate(scaled_error_cv  = sqrt(n1) * (theta_min_cv - theta_min_true))
results <- results %>% mutate(scaled_error_cv2 = sqrt(n1) * (theta_min_cv2 - theta_min_true))
# Visualize results
library(ggplot2)
library(tidyr)
ggplot(results) +
  geom_boxplot(aes(x = factor(n1), y = scaled_error)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  ylim(c(-2, 2)) +
  theme_mydefault(base_size = 20)

ggplot(results) +
  geom_boxplot(aes(x = factor(n1), y = scaled_error_cv)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  ylim(c(-2, 2)) +
  theme_mydefault(base_size = 20)

ggplot(results) +
  geom_boxplot(aes(x = factor(n1), y = scaled_error_cv2)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  theme_bw() +
  ylim(c(-2, 2)) +
  theme_mydefault(base_size = 20)


results %>% group_by(n1) %>%
  summarise(
            bias0 = mean(scaled_error),
            bias1 = mean(scaled_error_cv),
            bias2 = mean(scaled_error_cv2),
            sd0 = sd(scaled_error),
            sd1 = sd(scaled_error_cv),
            sd2 = sd(scaled_error_cv2)
            ) %>%
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

