rm(list = ls())
source("myggplot_theme.R", chdir = TRUE)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(latex2exp)
library(ggrepel)
library(patchwork)
library(svglite)
library(numDeriv)
set.seed(101)

##### PART 1 #####
# Define population and estimators
# Population
sample_space <- 1:3
probabilities <- c(1/3, 1/3, 1/3)
# probabilities <- c(0.4, 0.3, 0.3)
noise_y <- 0.1
regression_func <- function(x) {
  ifelse(x == 2, 1, 0)
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
    print(categorized1)
    print(categorized2)
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
      # obj[k] <- obj[k] + dnorm(categorized1$X[i], centers[k], sd_gau) * prop1[i] * mean_y1[i]
      obj[k] <- obj[k] + dnorm(categorized1$X[i], centers[k], sd_gau) * prop[i] * mean_y[i]
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
  theta_min <- sum(wts_min * prop1 * mean_y1)

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

find_wts <- function(h, categorized1, categorized2) {
  result <- solve_lp(h, categorized1, categorized2)
  result$wts_min
}

##### PART 2 #####
# Compare estimators

theta_min_true <- solve_lp(0, population, population)$theta_min

results_list <- list()
k <- 0

n_splits <- 10  # <-- number of random splits per bootstrap replicate

B <- 800
for (n1 in c(50, 100, 200, 400)) {

  theta_min_hat     <- rep(NA, B)
  theta_min_hat_cv  <- rep(NA, B)
  theta_min_hat_cv2 <- rep(NA, B)
  optimal_index_hat <- rep(NA, B)
  w0_min_hat  <- matrix(NA, nrow = B, ncol = length(sample_space))
  w0_min_dif  <- matrix(NA, nrow = B, ncol = length(sample_space))

  n1A <- round(n1 / 2)
  n1B <- n1 - n1A

  for (b in seq_len(B)) {
    ## 1) Generate full data of size n1
    trt <- get_data(n1)
    # Only used to compute a representative derivative
    # Not used in estimation
    summ_deriv <- get_data(n1) %>%
      group_by(X) %>%
      summarise(
                mean_y    = mean(Y),
                proportion = n() / n1,
                .groups = "drop"
      )

    ## Summarize full data
    summ_trt <- trt %>%
      group_by(X) %>%
      summarise(
                mean_y    = mean(Y),
                proportion = n() / n1,
                .groups = "drop"
      )

      ## 2) theta_min_hat from full data
      result                 <- solve_lp(0, summ_trt, summ_trt)
      theta_min_hat[b]       <- result$theta_min
      optimal_index_hat[b]   <- result$optimal_index
      w0_min_hat[b, ]        <- result$wts_min
      # Compute a representative derivative
      w0_min_dif[b, ]        <- jacobian(
                                  find_wts, 0,
                                  method = "Richardson", side = +1,
                                  categorized1 = summ_trt, categorized2 = summ_deriv
                                  # categorized1 = summ_trt, categorized2 = population
      )

      ## 3) Multiple random A/B splits for CV estimators
      theta_cv_splits  <- numeric(n_splits)
      theta_cv2_splits <- numeric(n_splits)

      for (s in seq_len(n_splits)) {
        # random split indices
        idx_A <- sample(seq_len(n1), n1A)
        idx_B <- setdiff(seq_len(n1), idx_A)

        trt_A <- trt[idx_A, ]
        trt_B <- trt[idx_B, ]

        summ_trt_A <- trt_A %>%
          group_by(X) %>%
          summarise(
                    mean_y    = mean(Y),
                    proportion = n() / n1A,
                    .groups = "drop"
                    ) %>%
          complete(X = sample_space,
                   fill = list(mean_y = 0, proportion = 0))

          summ_trt_B <- trt_B %>%
            group_by(X) %>%
            summarise(
                      mean_y    = mean(Y),
                      proportion = n() / n1B,
                      .groups = "drop"
                      ) %>%
            complete(X = sample_space,
                     fill = list(mean_y = 0, proportion = 0))
            # Plug-in weights: fit on A, evaluate on B
            result_A            <- solve_lp(0, summ_trt_A, summ_trt_A)
            theta_cv_splits[s]  <- sum(result_A$wts_min *
                                       summ_trt_B$proportion *
                                       summ_trt_B$mean_y)

            # Derivative based on hybrid population
            derivative          <- grad(
                                        find_theta, 0,
                                        method = "Richardson", side = +1,
                                        categorized1 = summ_trt_A, categorized2 = summ_trt_B
            )
            theta_cv2_splits[s] <- theta_cv_splits[s] + derivative
      }

      ## Store averages over splits
      theta_min_hat_cv[b]  <- mean(theta_cv_splits)
      theta_min_hat_cv2[b] <- mean(theta_cv2_splits)
  }

  # Put weights + optimal index into a data frame
  wts_df <- data.frame(
                       theta_min     = theta_min_hat,
                       optimal_index = optimal_index_hat,
                       theta_min_cv  = theta_min_hat_cv,
                       theta_min_cv2 = theta_min_hat_cv2
  )

  for (i in seq_along(sample_space)) {
    wts_df[, paste0("w0", i)] <- w0_min_hat[, i]
    wts_df[, paste0("w0", i, "_deriv")] <- w0_min_dif[, i]
  }

  wts_df$n1 <- n1
  k <- k + 1
  results_list[[k]] <- wts_df
}

results <- bind_rows(results_list)

results <- results %>%
  mutate(
         scaled_error     = sqrt(n1) * (theta_min     - theta_min_true),
         scaled_error_cv  = sqrt(n1) * (theta_min_cv  - theta_min_true),
         scaled_error_cv2 = sqrt(n1) * (theta_min_cv2 - theta_min_true),
         scaled_deriv     = sqrt(n1) * (theta_min_cv2 - theta_min_cv - theta_min_true)
  )

dir.create("results/example/v3/", recursive = TRUE, showWarnings = FALSE)
write.csv(results, "results/example/v3/results.csv", row.names = FALSE)
results <- read.csv("results/example/v3/results.csv", header = TRUE)
# Visualize results

results_long <- results %>%
  select(n1, scaled_error, scaled_error_cv, scaled_error_cv2) %>%
  pivot_longer(
    cols = c(scaled_error, scaled_error_cv, scaled_error_cv2),
    names_to = "method",
    values_to = "value"
  ) %>%
  mutate(method = factor(method,
                         levels = c("scaled_error",
                                    "scaled_error_cv",
                                    "scaled_error_cv2"),
                         labels = c("Empirical", "DS       ", "DS + Deriv")))



# ggplot(results) +
#   geom_boxplot(aes(x = factor(n1), y = scaled_error)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   theme_bw() +
#   ylim(c(-2, 2)) +
#   theme_mydefault(base_size = 20)

p1 <- ggplot(results_long, aes(x = factor(n1), y = value, fill = method)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "orange") +
  scale_fill_manual(
    values = c("yellow", "steelblue1", "steelblue4")  # optional
  ) +
  ylim(c(-2, 2)) +
  theme_mydefault(base_size = 20) +
  theme(legend.position = "top") +
  labs(
    # x = TeX(r"(Sample size $n_1$)"),
    x = "Sample size",
    # y = TeX(r"($\sqrt{n_1} \left( \hat{\theta}_{min} - \theta_{min} \right)$)"),
    y = "Scaled error",
    fill = "Estimator"
  )

print(p1)

# ggplot(results) +
#   geom_boxplot(aes(x = factor(n1), y = scaled_deriv)) +
#   geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
#   theme_bw() +
#   theme_mydefault(base_size = 20)
#
results %>% group_by(n1) %>%
  summarise(
            bias0 = mean(scaled_error),
            bias1 = mean(scaled_error_cv),
            bias2 = mean(scaled_error_cv2),
            sd0   = sd(scaled_error),
            sd1   = sd(scaled_error_cv),
            sd2   = sd(scaled_error_cv2),
            rmse0 = mean(scaled_error^2),
            rmse1 = mean(scaled_error_cv^2),
            rmse2 = mean(scaled_error_cv2^2)
            ) %>%
  print()

##### PART 3 #####
set.seed(111)
compute_vertices_analytic <- function(prop,
                                      sample_space = 1:3,
                                      centers      = 1:3,
                                      sd_gau       = 1) {
  stopifnot(length(prop) == length(sample_space))
  stopifnot(abs(sum(prop) - 1) < 1e-8)

  # A[i,k] = phi(x_i; center_k, sd_gau)
  A <- outer(sample_space, centers,
             function(x, m) dnorm(x, mean = m, sd = sd_gau))

  # a_k = sum_i A[i,k] * p_i
  a <- as.numeric(t(A) %*% prop)   # length K

  K <- length(centers)
  verts <- vector("list", K)

  for (k in seq_len(K)) {
    alpha_k <- 1 / a[k]
    # w_i = sum_k A[i,k] alpha_k = A[,k] * alpha_k (others 0)
    w <- A[, k] * alpha_k

    verts[[k]] <- data.frame(
                             w1     = w[1],
                             w2     = w[2],
                             w3     = w[3],
                             vertex = k
    )
  }

  do.call(rbind, verts)
}
verts_true <- compute_vertices_analytic(probabilities)
# Store vertices across replicates
B <- 800
verts_dat <- map_dfr(
                     .x = seq_len(B),
                     .f = function(b) {
                       # 1. Generate data
                       dat <- get_data(100)

                       # 2. Estimate prop_hat in the correct X order
                       prop_hat <- dat %>%
                         count(X) %>%
                         arrange(X) %>%
                         mutate(prop = n / sum(n)) %>%
                         pull(prop)

                       # 3. Compute vertices analytically and add replicate id
                       compute_vertices_analytic(prop_hat) %>%
                         mutate(rep = b)
                     }
)

# Ensure verts_true is ordered so the polygon doesn't self-cross
verts_true_plot <- verts_true %>%
  arrange(vertex)   # assuming vertex is 1,2,3

# ---- Arrow data (first 10 rows) ----
verts_opt <- results %>%
  filter(n1 == 100)

# Take first 10 rows
arrow_raw <- verts_opt %>%
  # filter(abs(w01_deriv) < 50, abs(w02_deriv) < 50, abs(w03_deriv) < 50) %>%
  slice_tail(n = 20) %>%
  mutate(deriv_norm = sqrt(w01_deriv^2 + w02_deriv^2))

# # Scale derivatives so the longest arrow has a fixed length on the plot
# max_len <- 0.10   # roughly 0.10 units in f-space; adjust if needed
# scale_fac <- max_len / max(arrow_raw$deriv_norm, na.rm = TRUE)
scale_fac <- 1

arrow_df <- arrow_raw %>%
  mutate(
         xend = w01 + scale_fac * w01_deriv,
         yend = w02 + scale_fac * w02_deriv
  )

# ---- Average arrow (over all rows with n1 == 200) ----
avg_arrow_df <- verts_opt %>%
  filter(abs(w01_deriv) < 5, abs(w02_deriv) < 5, abs(w03_deriv) < 5) %>%
  group_by(optimal_index) %>%
  summarise(
            w01       = mean(w01),
            w02       = mean(w02),
            w01_deriv = mean(w01_deriv),
            w02_deriv = mean(w02_deriv)
            ) %>%
  mutate(
         xend = w01 + scale_fac * w01_deriv,
         yend = w02 + scale_fac * w02_deriv
  )

p2 <- ggplot() +
  # True polytope (triangle) boundary
  geom_polygon(
    data  = verts_true_plot,
    aes(x = w1, y = w2, group = 1),
    fill  = "orange",
    # color = "orange",
    linewidth = 1,
    alpha = 0.5
  ) +
  # One representative \hat{W} polytope
  geom_polygon(
    data  = verts_dat %>% filter(rep == 1),
    aes(x = w1, y = w2, group = 1),
    fill  = "steelblue",
    color = "steelblue",
    linewidth = 0,
    alpha = 0.5
  ) +
  # 2D scatter-plot of estimated vertices across replicates
  geom_point(
    data = verts_dat,
    aes(x = w1, y = w2, group = vertex),
    alpha = 0.1,
    color = "steelblue2"
  ) +
  # Add mean of estimated vertices
  geom_point(data = verts_dat %>% group_by(vertex) %>% summarise(w1 = mean(w1), w2 = mean(w2)),
  aes(w1, w2), size = 3.5) +
  geom_hline(yintercept = theta_min_true / probabilities[2], color = "orange", linetype = "dashed") +
  coord_fixed() +
  labs(
    # title = TeX(r"(True polytope $W$ (orange), one estimated polytope $\hat{W}$ (blue) and distribution of estimated vertices)"),
    x     = expression(w(1)),
    y     = expression(w(2))
  ) +
  guides(fill = "none", color = "none") +
  theme_mydefault(base_size = 20) +
  xlim(c(0, 2.1)) + ylim(c(0.9, 1.5))

# Distribution of vertex after optimization
p3 <- ggplot() +
  # True polytope (triangle) boundary
  geom_polygon(
               data  = verts_true_plot,
               aes(x = w1, y = w2, group = 1),
               fill  = "orange",
               # color = "orange",
               linewidth = 1,
               alpha = 0.5
               ) +
  # Add regression function
  geom_hline(yintercept = theta_min_true / probabilities[2], color = "orange", linetype = "dashed") +
  # 2D scatter-plot of estimated vertices across replicates
  geom_point(
             data = verts_opt,
             aes(x = w01, y = w02, group = optimal_index),
             alpha = 0.1,
             color = "steelblue2"
             ) +
  # Black points = mean by optimal_index
  geom_point(
             data = verts_opt %>%
               group_by(optimal_index) %>%
               summarise(w01 = mean(w01), w02 = mean(w02), .groups = "drop"),
             aes(w01, w02),
             size = 3.5
             ) +
  # 10 representative arrows (transparent)
  # geom_segment(
  #              data = arrow_df,
  #              aes(x = w01, y = w02, xend = xend, yend = yend),
  #              arrow = arrow(length = unit(0.12, "cm")),
  #              linewidth = 0.4,
  #              alpha = 0.25,
  #              color = "black"
  #              ) +
  # Average arrow (solid), starting from mean(w01), mean(w02)
  geom_segment(
               data = avg_arrow_df,
               aes(x = w01, y = w02, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.38, "cm")),
               linewidth = 1.0,
               color = "black"
               ) +
  coord_fixed() +
  labs(
       # title = TeX(r"(True polytope $W$ (orange) and distribution of estimated optimal vertices)"),
       x     = expression(w(1)),
       y     = expression(w(2))
       ) +
  guides(fill = "none", color = "none") +
  theme_mydefault(base_size = 20) +
  xlim(c(0, 2.1)) + ylim(c(0.9, 1.5))

# ---- Combine the three plots vertically ----
p <- (p1 / (p2 | p3)) +
  plot_layout(heights = c(1, 1), widths = c(1, 1.8)) +
  theme(
        plot.margin = margin(0,0,0,0),
        axis.title.x = element_text(margin = margin(t=0)),
        axis.title.y = element_text(margin = margin(r=0))
  )
print(p)

# Save the three panels
ggsave("results/example/A.svg", plot = p1,
       width = 7.5, height = 10, units = "in")
ggsave("results/example/B.svg", plot = p2,
       width = 15, height = 5, units = "in")
ggsave("results/example/C.svg", plot = p3,
       width = 15, height = 5, units = "in")

# Save Screenshot
