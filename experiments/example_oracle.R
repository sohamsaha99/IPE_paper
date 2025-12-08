source("myggplot_theme.R", chdir = TRUE)
library(dplyr)
library(numDeriv)
set.seed(101)
# Population
sample_space <- 1:3
probabilities <- c(1/3, 1/3, 1/3)
# probabilities <- c(0.4, 0.3, 0.3)
noise_y <- 1
regression_func <- function(x) {
  ifelse(x == 3, 0, 1)
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

solve_lp <- function(h, categorized1, categorized2, L = 1.0) {
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
  K       <- length(prop1)

  if (length(mean_y1) != K || length(prop2) != K || length(mean_y2) != K) {
    stop("Length mismatch. Stop.")
  }
  prop   <- (1 - h) * prop1 + h * prop2
  mean_y <- (1 - h) * mean_y1 + h * mean_y2
  # Objective: minimize sum_i w_i * p1_i * y1_i
  obj    <- prop1 * mean_y1

  # Constraint 1: sum_i w_i * p_i = 1
  A   <- matrix(prop, nrow = 1)
  dir <- "=="
  rhs <- 1

  # for (i in seq_len(K)) {
  #   row    <- numeric(K)
  #   row[i] <- 1
  #   A      <- rbind(A, row)
  #   dir    <- c(dir, "<=")
  #   rhs    <- c(rhs, L)
  # }
  # # Constraint 2: each w_i ph_i >= L min{ph_i}
  # for (i in seq_len(K)) {
  #   row    <- numeric(K)
  #   row[i] <- prop[i]
  #   A      <- rbind(A, row)
  #   dir    <- c(dir, ">=")
  #   rhs    <- c(rhs, min(prop) * L)
  # }
  # # Constraint 2: ratio constraints w_i >= L * w_j for all i, j
  # # Implemented as: w_i - L * w_j >= 0
  # for (i in seq_len(K)) {
  #   for (j in seq_len(K)) {
  #     if (i == j) next
  #     row      <- numeric(K)
  #     # row[i] <- prop[i]
  #     row[i]   <- 1
  #     # row[j] <- -ratio_bound * prop[j]
  #     row[j]   <- -L
  #     A        <- rbind(A, row)
  #     dir      <- c(dir, ">=")
  #     rhs      <- c(rhs, 0)
  #   }
  # }
  # # Constraint 2: difference constraints w_i + w_j >= L for all i!=j
  # for (i in seq_len(K)) {
  #   for (j in seq_len(K)) {
  #     if (i == j) next
  #     row      <- numeric(K)
  #     # row[i] <- prop[i]
  #     row[i]   <- 1
  #     # row[j] <- -ratio_bound * prop[j]
  #     row[j]   <- 1
  #     A        <- rbind(A, row)
  #     dir      <- c(dir, ">=")
  #     rhs      <- c(rhs, L)
  #   }
  # }
  # Constraint 2: difference constraints p_iw_i + p_jw_j <= L for all i!=j
  for (i in seq_len(K)) {
    for (j in seq_len(K)) {
      if (i == j) next
      row      <- numeric(K)
      # row[i] <- prop[i]
      row[i]   <- 1
      # row[j] <- -ratio_bound * prop[j]
      row[j]   <- 1
      A        <- rbind(A, row)
      dir      <- c(dir, ">=")
      rhs      <- c(rhs, L)
    }
  }


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

  wts_min <- sol$solution
  theta_min <- sum(prop1 * mean_y1 * wts_min)  # optimal objective value

  # "Optimal index": here I return the index with the largest weight
  optimal_index <- which.max(wts_min)

  list(
       theta_min     = theta_min,
       wts_min       = wts_min,
       optimal_index = optimal_index
  )
}
find_theta <- function(h, categorized1, categorized2) {
  solve_lp(h, categorized1, categorized2)$theta_min
}

library(ggplot2)

# Visualize feasible sets
feasible_set <- function(prop, L, prop_hat = NULL, res = 0.001) {
  stopifnot(length(prop) == 3, abs(sum(prop) - 1) <= 1e-8)
  if (!is.null(prop_hat)) {
    stopifnot(length(prop_hat) == 3, abs(sum(prop_hat) - 1) <= 1e-8)
  }

  # Grid over q1 and q2
  q1_vals <- seq(0, 1, by = res)
  q2_vals <- seq(0, 1, by = res)
  grid <- expand.grid(q1 = q1_vals, q2 = q2_vals)

  # Compute q3 and simplex constraint
  grid$q3 <- 1 - grid$q1 - grid$q2
  simplex <- grid$q3 >= 0

  # Helper to compute feasibility for a given prop
  compute_feasible <- function(p) {
    w1 <- grid$q1 / p[1]
    w2 <- grid$q2 / p[2]
    w3 <- grid$q3 / p[3]
    wts <- list(w1, w2, w3)

    # # Constraint: w_i / w_j >= L
    # feasible <- simplex
    # for (i in 1:3) {
    #   for (j in 1:3) {
    #     feasible <- feasible & (wts[[i]] / wts[[j]] >= L)
    #   }
    # }
    # # Constraint: w_i - w_j <= L
    # feasible <- simplex
    # for (i in 1:3) {
    #   for (j in 1:3) {
    #     if (i == j) next
    #     feasible <- feasible & (wts[[i]] - wts[[j]] <= L)
    #   }
    # }
    # Constraint: (p_i + p_j)/L <= q_i + q_j <= L (p_i + p_j)
    feasible <- simplex
    for (i in 1:3) {
      for (j in 1:3) {
        if (i == j) next
        feasible <- feasible & (wts[[i]] * p[i] + wts[[j]] * p[j] <= L * (p[i] + p[j]))
        feasible <- feasible & (wts[[i]] * p[i] + wts[[j]] * p[j] >= (p[i] + p[j]) / L)
      }
    }

    feasible
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

  # Helper to add boundary lines for a given prop
  add_boundary_lines <- function(p, pvec, col) {
    # # w1 / w2 = L: q2 = (p2 / (p1 L)) q1
    # p <- p + geom_abline(
    #                      intercept = 0,
    #                      slope     = pvec[2] / (pvec[1] * L),
    #                      colour    = col,
    #                      linewidth = 0.7
    # )
    # # w2 / w1 = L: q2 = (p2 L / p1) q1
    # p <- p + geom_abline(
    #                      intercept = 0,
    #                      slope     = L * pvec[2] / pvec[1],
    #                      colour    = col,
    #                      linewidth = 0.7
    # )
    # # w1 / w3 = L: q2 = 1 - (1 + p3 / (L p1)) q1
    # p <- p + geom_abline(
    #                      intercept = 1,
    #                      slope     = -1 - pvec[3] / (L * pvec[1]),
    #                      colour    = col,
    #                      linewidth = 0.7
    # )
    # # w3 / w1 = L: q2 = 1 - (1 + p3 L / p1) q1
    # p <- p + geom_abline(
    #                      intercept = 1,
    #                      slope     = -1 - L * pvec[3] / pvec[1],
    #                      colour    = col,
    #                      linewidth = 0.7
    # )
    # w2 + w1 = L
    # => q2 = (- p3 L + 1) - (1 + p3/p1) q1
    p <- p + geom_abline(
                         intercept = - pvec[3] * L + 1,
                         slope     = -1 - pvec[3] / pvec[1],
                         colour    = col,
                         linewidth = 0.7
    )
    p
  }

  # Add boundaries: blue for true, orange for estimated
  p <- add_boundary_lines(p, prop, col = "blue")
  if (!is.null(prop_hat)) {
    p <- add_boundary_lines(p, prop_hat, col = "orange")
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
         x = expression(q[1]==w[1]*p[1]),
         y = expression(q[2]==w[2]*p[2])
         ) +
    theme_mydefault(base_size = 20)
}

feasible_set(prop = probabilities,
             L = 1.2,
             prop_hat = get_data(50) %>% group_by(X) %>% summarise(prop = n() / 50) %>% pull(prop)) +
# q1 m1 + q2 m2 = c
geom_abline(intercept = 1 / (3*regression_func(2)), slope = - regression_func(1) / regression_func(2), color = "red")

library(ggplot2)
library(ggrepel)

feasible_set <- function(prop_pop, L, prop_est = NULL) {
  stopifnot(length(prop_pop) == 3)
  if (!is.null(prop_est)) stopifnot(length(prop_est) == 3)

  lower_from_prop <- function(p) L * min(p)

  make_vertices <- function(lower) {
    if (3 * lower > 1) return(NULL)
    data.frame(
               q1 = c(lower, 1 - 2 * lower, lower),
               q2 = c(lower, lower,         1 - 2 * lower)
    )
  }

  # Population polytope
  lower_pop <- lower_from_prop(prop_pop)
  verts_pop <- make_vertices(lower_pop)
  if (!is.null(verts_pop)) verts_pop$type <- "Population"

  # Estimated polytope
  verts_est <- NULL
  if (!is.null(prop_est)) {
    lower_est <- lower_from_prop(prop_est)
    verts_est <- make_vertices(lower_est)
    if (!is.null(verts_est)) verts_est$type <- "Estimated"
  }

  vertices <- rbind(verts_pop, verts_est)
  if (is.null(vertices)) stop("No feasible region.")

  # Text labels
  labels_df <- vertices
  labels_df$label <- sprintf("(%.2f, %.2f)", labels_df$q1, labels_df$q2)

  # Colors for points + labels
  vertex_colors <- c(
                     "Population" = "navy",
                     "Estimated"  = "darkred"
  )

  ggplot(vertices, aes(x = q1, y = q2)) +

    # Polytopes (same fill scheme as before)
    geom_polygon(aes(fill = type, group = type),
                 alpha = 0.35, color = "black") +

    # Vertices (colored by type)
    geom_point(
               aes(color = type),
               size = 3,
               show.legend = FALSE
               ) +

    # Non-overlapping text labels (colored by type)
    geom_text_repel(
                    data = labels_df,
                    aes(label = label, color = type),
                    size = 10,
                    nudge_x = 0.03,
                    nudge_y = 0.03,
                    max.overlaps = Inf,
                    show.legend = FALSE,
                    segment.color = NA
                    ) +

    scale_fill_manual(values = c("Population" = "steelblue",
                                 "Estimated"  = "tomato")) +
    scale_color_manual(values = vertex_colors) +

    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
         title = "Feasible Polytopes",
         x = expression(q[1] == w[1] * p[1]),
         y = expression(q[2] == w[2] * p[2]),
         fill = "Polytope"
         ) +
    theme_mydefault(base_size = 20) +
    theme(legend.position = "bottom")
}

feasible_set(prop_pop = probabilities,
             L = 0.5,
             prop_est = as.numeric(rmultinom(1, size = 50, prob = probabilities)) / 50)

theta_min_true <- solve_lp(0, population, population)$theta_min

results_list <- list()
k <- 0
for (n1 in c(100, 200, 400, 800)) {
  B <- 800
  theta_min_hat     <- rep(NA, B)
  theta_min_hat_cv  <- rep(NA, B)
  theta_min_hat_cv2 <- rep(NA, B)
  optimal_index_hat <- rep(NA, B)
  wts_min_hat <- matrix(NA, nrow = B, ncol = length(sample_space))
  for (b in seq_len(B)) {
    trt <- get_data(n1)
    # Generate split B
    trt_B <- get_data(n1)
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
      # Solve linear program
      result               <- solve_lp(0, summ_trt, summ_trt)
      theta_min_hat[b]     <- result$theta_min
      optimal_index_hat[b] <- result$optimal_index
      wts_min_hat[b, ]     <- result$wts_min
      # Plug-in weights on data B
      theta_min_hat_cv[b]  <- sum(result$wts_min * summ_trt_B$proportion * summ_trt_B$mean_y)
      # Find derivative based on hybrid population
      derivative           <- grad(find_theta, 0, method = "Richardson", side = +1,
                                   categorized1 = summ_trt, categorized2 = summ_trt_B)
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

