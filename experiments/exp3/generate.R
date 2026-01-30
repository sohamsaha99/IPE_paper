rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(latex2exp)
library(ggpubr)
library(mvtnorm)
source("myggplot_theme.R")

FILE_PATH <- "./experiments/exp3/GBM_files/GBM_complete.csv"
df_GBM <- read.csv(FILE_PATH, header = TRUE, stringsAsFactors = FALSE)
# "healthy" used for resampling
df_GBM <- df_GBM %>% mutate(healthy = (kps == 1 & eor == "GTR" & age < 55))

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

df_GBM %>%
  ggplot(aes(x = healthy, y = event_time, fill = healthy)) +
  geom_boxplot() +
  ylim(c(0, 40)) +
  theme_pubr() +
  scale_fill_discrete() +
  labs(y = "Overall Survival (months)",
       x = "kps > 90, eor = \"GTR\", age < 55") +
  theme(legend.position = "none") +
  theme_mydefault(base_size = 20)

# Rescale all predictor columns so that they have mean 0 variance 1
df_GBM_processed <- df_GBM %>% mutate(isGTR = as.numeric(eor == "GTR")) %>%
  select(age, sex, kps, mgmt, isGTR, event_time, healthy) %>%
  rename(Y = event_time) %>%
  mutate(across(c(age, sex, kps, mgmt, isGTR),
                ~ (. - mean(., na.rm = TRUE)) / sd(., na.rm = TRUE)))

# Truncate outcome at 24 months
df_GBM_processed <- df_GBM_processed %>% mutate(Y = pmin(Y, 24))
# Check normalization
df_GBM_processed %>%
  summarise(across(
    everything(),
    list(
      mean = ~ mean(., na.rm = TRUE),
      var  = ~ var(., na.rm = TRUE)
    )
  ))

# Define sampling scheme
generate_data <- function(dat, n0, n1) {
  # Sampling WITH replacement
  # Treatment group
  idx_trt <- sample(nrow(dat), size = n1, replace = TRUE,
                    # prob = ifelse(dat$healthy, 0.75, 0.75))
                    prob = rep(1, nrow(dat)))
  trt <- dat[idx_trt, ] %>% select(age, sex, kps, mgmt, isGTR, Y)
  # Control group
  idx_con <- sample(nrow(dat), size = n0, replace = TRUE,
                    # prob = ifelse(dat$healthy, 0.25, 0.25))
                    prob = rep(1, nrow(dat)))
  con <- dat[idx_con, ] %>% select(age, sex, kps, mgmt, isGTR, Y)
  # Compute summary
  target <- data.frame(matrix(NA, nrow = 6, ncol = 3))
  colnames(target) <- c("xname", "value", "type")
  target$xname <- c("age", "age", "sex", "kps", "mgmt", "isGTR")
  target$type <- c("first", "second", rep("first", 4))
  target$value <- c(mean(con$age), mean(con$age^2), mean(con$sex), mean(con$kps), mean(con$mgmt), mean(con$isGTR))
  control_mean <- mean(con$Y)
  list(trt = trt, target = target, control_mean = control_mean)
}

# Find true IPE with Gaussian mixture regularization
# Create population
population_trt <- df_GBM_processed %>%
  # mutate(prob = ifelse(healthy, 0.75, 0.75)) %>%
  mutate(prob = 1) %>%
  mutate(prob = prob / sum(prob)) %>%
  select(age, sex, kps, mgmt, isGTR, Y, prob)
population_con <- df_GBM_processed %>%
  # mutate(prob = ifelse(healthy, 0.25, 0.25)) %>%
  mutate(prob = 1) %>%
  mutate(prob = prob / sum(prob)) %>%
  select(age, sex, kps, mgmt, isGTR, Y, prob)
# Create control population summary
target_con <- data.frame(matrix(NA, nrow = 6, ncol = 3))
colnames(target_con) <- c("xname", "value", "type")
target_con$xname     <- c("age", "age", "sex", "kps", "mgmt", "isGTR")
target_con$type      <- c("first", "second", rep("first", 4))
target_con$value <- c(sum(population_con$age   * population_con$prob),
                      sum(population_con$age^2 * population_con$prob),
                      sum(population_con$sex   * population_con$prob),
                      sum(population_con$kps   * population_con$prob),
                      sum(population_con$mgmt  * population_con$prob),
                      sum(population_con$isGTR * population_con$prob)
)
population_con_Y_mean <- sum(population_con$Y * population_con$prob)
grid_points <- df_GBM_processed %>% select(age, sex, kps, mgmt, isGTR)
grid_points <- expand.grid(age   = seq(min(df_GBM_processed$age), max(df_GBM_processed$age), length.out = 50),
                           sex   = unique(df_GBM_processed$sex),
                           kps   = unique(df_GBM_processed$kps),
                           mgmt  = unique(df_GBM_processed$mgmt),
                           isGTR = unique(df_GBM_processed$isGTR)
)
grid_points <- data.frame(grid_points)

# Function for finding true theta_min
true_theta_min <- function(population_trt, grid_points, L, target_con, population_con_Y_mean) {
  # Compute dmvnorm
  normal_den <- function(x, v) {
    stopifnot(ncol(x) == length(v))
    dmvnorm(x, mean = v, sigma = diag(length(v)) / (L^2))
  }
  X_cols <- c("age", "sex", "kps", "mgmt", "isGTR")
  # Set linear programming:
  # min c^\top alpha
  # subject to A alpha = target_con$value
  K <- nrow(grid_points)
  cvec <- rep(NA, K)
  for (k in seq_len(K)) {
    cvec[k] <- sum(population_trt$Y * population_trt$prob *
                   normal_den(population_trt[, X_cols], as.numeric(grid_points[k, ])))
  }
  rhsvec <- c(1, target_con$value)
  dirvec <- rep("==", length(rhsvec))
  Amat <- matrix(NA, nrow = length(rhsvec), ncol = K)
  # Set up each row
  for (k in seq_len(K)) {
    den_k <- normal_den(population_trt[, X_cols], as.numeric(grid_points[k, ]))
    stopifnot(length(den_k) == nrow(population_trt))
    Amat[1, k] <- sum(population_trt$prob  * den_k)
    Amat[2, k] <- sum(population_trt$age   * population_trt$prob * den_k)
    Amat[3, k] <- sum(population_trt$age^2 * population_trt$prob * den_k)
    Amat[4, k] <- sum(population_trt$sex   * population_trt$prob * den_k)
    Amat[5, k] <- sum(population_trt$kps   * population_trt$prob * den_k)
    Amat[6, k] <- sum(population_trt$mgmt  * population_trt$prob * den_k)
    Amat[7, k] <- sum(population_trt$isGTR * population_trt$prob * den_k)
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
  theta_min <- EY1_min - population_con_Y_mean

  list(theta_min = theta_min,
       alpha = alpha_hat
  )
}

true_theta_max <- function(population_trt, grid_points, L, target_con, population_con_Y_mean) {
  pop_trt_neg <- population_trt
  pop_trt_neg$Y <- -pop_trt_neg$Y

  ret <- true_theta_min(pop_trt_neg,
                        grid_points,
                        L,
                        target_con,
                        population_con_Y_mean = -population_con_Y_mean)

  list(theta_max = -ret$theta_min,
       alpha = ret$alpha)
}

true_theta_min(population_trt, grid_points, 0.6, target_con, population_con_Y_mean)
true_theta_max(population_trt, grid_points, 1, target_con, population_con_Y_mean)

# Plot for range of values of L
library(ggplot2)

# Grid of L values
L_seq <- seq(0.5, 2, length.out = 20)

# Compute theta_min and theta_max
theta_min <- sapply(L_seq, function(L)
  true_theta_min(population_trt, grid_points, L, target_con, population_con_Y_mean)$theta_min
)

theta_max <- sapply(L_seq, function(L)
  true_theta_max(population_trt, grid_points, L, target_con, population_con_Y_mean)$theta_max
)

# Data frame for plotting
df <- data.frame(
  L = L_seq,
  theta_min = theta_min,
  theta_max = theta_max
)

# Plot
p1 <- ggplot(df, aes(x = L)) +
  geom_ribbon(aes(ymin = theta_min, ymax = theta_max), alpha = 0.3) +
  geom_line(aes(y = theta_min), linetype = "dashed") +
  geom_line(aes(y = theta_max), linetype = "dashed") +
  labs(
       x = "L",
       y = expression(theta),
       title = "Usual Summary"
       ) +
  ylim(c(-4, 4)) +
  theme_minimal() +
  theme_mydefault(base_size = 20)
print(p1)

### Version 2 with more summary
# Create control population summary
age_breaks <- c(-Inf, -0.508, 0.085, 0.777, Inf)
target_con_v2 <- data.frame(matrix(NA, nrow = 9, ncol = 3))
colnames(target_con_v2) <- c("xname", "value", "type")
target_con_v2$xname     <- c("age", "age", "sex", "kps", "mgmt", "isGTR", rep("age", 3))
target_con_v2$type      <- c("first", "second", rep("first", 4 + 3))
target_con_v2$value <- c(sum(population_con$age   * population_con$prob),
                         sum(population_con$age^2 * population_con$prob),
                         sum(population_con$sex   * population_con$prob),
                         sum(population_con$kps   * population_con$prob),
                         sum(population_con$mgmt  * population_con$prob),
                         sum(population_con$isGTR * population_con$prob),
                         sum((population_con$age > age_breaks[1] & population_con$age <= age_breaks[2]) * population_con$prob),
                         sum((population_con$age > age_breaks[2] & population_con$age <= age_breaks[3]) * population_con$prob),
                         sum((population_con$age > age_breaks[3] & population_con$age <= age_breaks[4]) * population_con$prob)
)

# Function for finding true theta_min
true_theta_min_v2 <- function(population_trt, grid_points, L, target_con, population_con_Y_mean, age_breaks) {
  # Compute dmvnorm
  normal_den <- function(x, v) {
    stopifnot(ncol(x) == length(v))
    dmvnorm(x, mean = v, sigma = diag(length(v)) / (L^2))
  }
  X_cols <- c("age", "sex", "kps", "mgmt", "isGTR")
  # Set linear programming:
  # min c^\top alpha
  # subject to A alpha = target_con$value
  K <- nrow(grid_points)
  cvec <- rep(NA, K)
  for (k in seq_len(K)) {
    cvec[k] <- sum(population_trt$Y * population_trt$prob *
                   normal_den(population_trt[, X_cols], as.numeric(grid_points[k, ])))
  }
  rhsvec <- c(1, target_con$value)
  dirvec <- rep("==", length(rhsvec))
  Amat <- matrix(NA, nrow = length(rhsvec), ncol = K)
  # Set up each row
  for (k in seq_len(K)) {
    den_k <- normal_den(population_trt[, X_cols], as.numeric(grid_points[k, ]))
    stopifnot(length(den_k) == nrow(population_trt))
    Amat[1, k] <- sum(population_trt$prob  * den_k)
    Amat[2, k] <- sum(population_trt$age   * population_trt$prob * den_k)
    Amat[3, k] <- sum(population_trt$age^2 * population_trt$prob * den_k)
    Amat[4, k] <- sum(population_trt$sex   * population_trt$prob * den_k)
    Amat[5, k] <- sum(population_trt$kps   * population_trt$prob * den_k)
    Amat[6, k] <- sum(population_trt$mgmt  * population_trt$prob * den_k)
    Amat[7, k] <- sum(population_trt$isGTR * population_trt$prob * den_k)
    Amat[8, k] <- sum((population_trt$age > age_breaks[1] & population_trt$age <= age_breaks[2]) * population_trt$prob * den_k)
    Amat[9, k] <- sum((population_trt$age > age_breaks[2] & population_trt$age <= age_breaks[3]) * population_trt$prob * den_k)
    Amat[10, k] <- sum((population_trt$age > age_breaks[3] & population_trt$age <= age_breaks[4]) * population_trt$prob * den_k)
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
  theta_min <- EY1_min - population_con_Y_mean

  list(theta_min = theta_min,
       alpha = alpha_hat
  )
}

true_theta_max_v2 <- function(population_trt, grid_points, L, target_con, population_con_Y_mean, age_breaks) {
  pop_trt_neg <- population_trt
  pop_trt_neg$Y <- -pop_trt_neg$Y

  ret <- true_theta_min_v2(pop_trt_neg,
                           grid_points,
                           L,
                           target_con,
                           population_con_Y_mean = -population_con_Y_mean,
                           age_breaks = age_breaks)

  list(theta_max = -ret$theta_min,
       alpha = ret$alpha)
}

# Compute theta_min and theta_max
theta_min <- sapply(L_seq, function(L)
  true_theta_min_v2(population_trt, grid_points, L, target_con_v2, population_con_Y_mean, age_breaks)$theta_min
)

theta_max <- sapply(L_seq, function(L)
  true_theta_max_v2(population_trt, grid_points, L, target_con_v2, population_con_Y_mean, age_breaks)$theta_max
)

# Data frame for plotting
df <- data.frame(
  L = L_seq,
  theta_min = theta_min,
  theta_max = theta_max
)

# Plot
p2 <- ggplot(df, aes(x = L)) +
  geom_ribbon(aes(ymin = theta_min, ymax = theta_max), alpha = 0.3) +
  geom_line(aes(y = theta_min), linetype = "dashed") +
  geom_line(aes(y = theta_max), linetype = "dashed") +
  labs(
       x = "L",
       y = expression(theta),
       title = "Usual summary + age histogram (4 bins)"
       ) +
  ylim(c(-4, 4)) +
  theme_minimal() +
  theme_mydefault(base_size = 20)

print(p2)

### Version 3 with more summary
# Create control population summary
target_con_v3 <- data.frame(matrix(NA, nrow = 7, ncol = 3))
colnames(target_con_v3) <- c("xname", "value", "type")
target_con_v3$xname     <- c("age", "age", "sex", "kps", "mgmt", "isGTR", "kps1isGTR")
target_con_v3$type      <- c("first", "second", rep("first", 4 + 1))
target_con_v3$value <- c(sum(population_con$age   * population_con$prob),
                         sum(population_con$age^2 * population_con$prob),
                         sum(population_con$sex   * population_con$prob),
                         sum(population_con$kps   * population_con$prob),
                         sum(population_con$mgmt  * population_con$prob),
                         sum(population_con$isGTR * population_con$prob),
                         sum((population_con$isGTR > 0 & population_con$kps > 0) * population_con$prob)
)

# Function for finding true theta_min
true_theta_min_v3 <- function(population_trt, grid_points, L, target_con, population_con_Y_mean) {
  # Compute dmvnorm
  normal_den <- function(x, v) {
    stopifnot(ncol(x) == length(v))
    dmvnorm(x, mean = v, sigma = diag(length(v)) / (L^2))
  }
  X_cols <- c("age", "sex", "kps", "mgmt", "isGTR")
  # Set linear programming:
  # min c^\top alpha
  # subject to A alpha = target_con$value
  K <- nrow(grid_points)
  cvec <- rep(NA, K)
  for (k in seq_len(K)) {
    cvec[k] <- sum(population_trt$Y * population_trt$prob *
                   normal_den(population_trt[, X_cols], as.numeric(grid_points[k, ])))
  }
  rhsvec <- c(1, target_con$value)
  dirvec <- rep("==", length(rhsvec))
  Amat <- matrix(NA, nrow = length(rhsvec), ncol = K)
  # Set up each row
  for (k in seq_len(K)) {
    den_k <- normal_den(population_trt[, X_cols], as.numeric(grid_points[k, ]))
    stopifnot(length(den_k) == nrow(population_trt))
    Amat[1, k] <- sum(population_trt$prob  * den_k)
    Amat[2, k] <- sum(population_trt$age   * population_trt$prob * den_k)
    Amat[3, k] <- sum(population_trt$age^2 * population_trt$prob * den_k)
    Amat[4, k] <- sum(population_trt$sex   * population_trt$prob * den_k)
    Amat[5, k] <- sum(population_trt$kps   * population_trt$prob * den_k)
    Amat[6, k] <- sum(population_trt$mgmt  * population_trt$prob * den_k)
    Amat[7, k] <- sum(population_trt$isGTR * population_trt$prob * den_k)
    Amat[8, k] <- sum((population_trt$isGTR > 0 & population_trt$kps > 0) * population_trt$prob * den_k)
  }

  # print(cvec)
  # print(Amat)
  # print(rhsvec)
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
  theta_min <- EY1_min - population_con_Y_mean

  list(theta_min = theta_min,
       alpha = alpha_hat
  )
}

true_theta_max_v3 <- function(population_trt, grid_points, L, target_con, population_con_Y_mean) {
  pop_trt_neg <- population_trt
  pop_trt_neg$Y <- -pop_trt_neg$Y

  ret <- true_theta_min_v3(pop_trt_neg,
                           grid_points,
                           L,
                           target_con,
                           population_con_Y_mean = -population_con_Y_mean)

  list(theta_max = -ret$theta_min,
       alpha = ret$alpha)
}

# Compute theta_min and theta_max
theta_min <- sapply(L_seq, function(L)
  true_theta_min_v3(population_trt, grid_points, L, target_con_v3, population_con_Y_mean)$theta_min
)

theta_max <- sapply(L_seq, function(L)
  true_theta_max_v3(population_trt, grid_points, L, target_con_v3, population_con_Y_mean)$theta_max
)

# Data frame for plotting
df <- data.frame(
  L = L_seq,
  theta_min = theta_min,
  theta_max = theta_max
)

# Plot
p3 <- ggplot(df, aes(x = L)) +
  geom_ribbon(aes(ymin = theta_min, ymax = theta_max), alpha = 0.3) +
  geom_line(aes(y = theta_min), linetype = "dashed") +
  geom_line(aes(y = theta_max), linetype = "dashed") +
  labs(
       x = "L",
       y = expression(theta),
       title = "Usual summary + Joint table (kps, eor)"
       ) +
  ylim(c(-4, 4)) +
  theme_minimal() +
  theme_mydefault(base_size = 20)

print(p3)

### Version 4 with most summary
# Create control population summary
age_breaks <- c(-Inf, -0.508, 0.085, 0.777, Inf)
target_con_v4 <- data.frame(matrix(NA, nrow = 10, ncol = 3))
colnames(target_con_v4) <- c("xname", "value", "type")
target_con_v4$xname     <- c("age", "age", "sex", "kps", "mgmt", "isGTR", rep("age", 3), "kps1isGTR")
target_con_v4$type      <- c("first", "second", rep("first", 4 + 3 + 1))
target_con_v4$value <- c(sum(population_con$age   * population_con$prob),
                         sum(population_con$age^2 * population_con$prob),
                         sum(population_con$sex   * population_con$prob),
                         sum(population_con$kps   * population_con$prob),
                         sum(population_con$mgmt  * population_con$prob),
                         sum(population_con$isGTR * population_con$prob),
                         sum((population_con$age > age_breaks[1] & population_con$age <= age_breaks[2]) * population_con$prob),
                         sum((population_con$age > age_breaks[2] & population_con$age <= age_breaks[3]) * population_con$prob),
                         sum((population_con$age > age_breaks[3] & population_con$age <= age_breaks[4]) * population_con$prob),
                         sum((population_con$isGTR > 0 & population_con$kps > 0) * population_con$prob)
)

# Function for finding true theta_min
true_theta_min_v4 <- function(population_trt, grid_points, L, target_con, population_con_Y_mean, age_breaks) {
  # Compute dmvnorm
  normal_den <- function(x, v) {
    stopifnot(ncol(x) == length(v))
    dmvnorm(x, mean = v, sigma = diag(length(v)) / (L^2))
  }
  X_cols <- c("age", "sex", "kps", "mgmt", "isGTR")
  # Set linear programming:
  # min c^\top alpha
  # subject to A alpha = target_con$value
  K <- nrow(grid_points)
  cvec <- rep(NA, K)
  for (k in seq_len(K)) {
    cvec[k] <- sum(population_trt$Y * population_trt$prob *
                   normal_den(population_trt[, X_cols], as.numeric(grid_points[k, ])))
  }
  rhsvec <- c(1, target_con$value)
  dirvec <- rep("==", length(rhsvec))
  Amat <- matrix(NA, nrow = length(rhsvec), ncol = K)
  # Set up each row
  for (k in seq_len(K)) {
    den_k <- normal_den(population_trt[, X_cols], as.numeric(grid_points[k, ]))
    stopifnot(length(den_k) == nrow(population_trt))
    Amat[1, k] <- sum(population_trt$prob  * den_k)
    Amat[2, k] <- sum(population_trt$age   * population_trt$prob * den_k)
    Amat[3, k] <- sum(population_trt$age^2 * population_trt$prob * den_k)
    Amat[4, k] <- sum(population_trt$sex   * population_trt$prob * den_k)
    Amat[5, k] <- sum(population_trt$kps   * population_trt$prob * den_k)
    Amat[6, k] <- sum(population_trt$mgmt  * population_trt$prob * den_k)
    Amat[7, k] <- sum(population_trt$isGTR * population_trt$prob * den_k)
    Amat[8, k] <- sum((population_trt$age > age_breaks[1] & population_trt$age <= age_breaks[2]) * population_trt$prob * den_k)
    Amat[9, k] <- sum((population_trt$age > age_breaks[2] & population_trt$age <= age_breaks[3]) * population_trt$prob * den_k)
    Amat[10, k] <- sum((population_trt$age > age_breaks[3] & population_trt$age <= age_breaks[4]) * population_trt$prob * den_k)
    Amat[11, k] <- sum((population_trt$isGTR > 0 & population_trt$kps > 0) * population_trt$prob * den_k)
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
  theta_min <- EY1_min - population_con_Y_mean

  list(theta_min = theta_min,
       alpha = alpha_hat
  )
}

true_theta_max_v4 <- function(population_trt, grid_points, L, target_con, population_con_Y_mean, age_breaks) {
  pop_trt_neg <- population_trt
  pop_trt_neg$Y <- -pop_trt_neg$Y

  ret <- true_theta_min_v4(pop_trt_neg,
                           grid_points,
                           L,
                           target_con,
                           population_con_Y_mean = -population_con_Y_mean,
                           age_breaks = age_breaks)

  list(theta_max = -ret$theta_min,
       alpha = ret$alpha)
}

# Compute theta_min and theta_max
theta_min <- sapply(L_seq, function(L)
  true_theta_min_v4(population_trt, grid_points, L, target_con_v4, population_con_Y_mean, age_breaks)$theta_min
)

theta_max <- sapply(L_seq, function(L)
  true_theta_max_v4(population_trt, grid_points, L, target_con_v4, population_con_Y_mean, age_breaks)$theta_max
)

# Data frame for plotting
df <- data.frame(
  L = L_seq,
  theta_min = theta_min,
  theta_max = theta_max
)

# Plot
p4 <- ggplot(df, aes(x = L)) +
  geom_ribbon(aes(ymin = theta_min, ymax = theta_max), alpha = 0.3) +
  geom_line(aes(y = theta_min), linetype = "dashed") +
  geom_line(aes(y = theta_max), linetype = "dashed") +
  labs(
       x = "L",
       y = expression(theta),
       title = "Usual summary + Joint table + age histogram"
       ) +
  ylim(c(-4, 4)) +
  theme_minimal() +
  theme_mydefault(base_size = 20)

print(p4)

# Make one grid:
library(patchwork)
p1 + p2 + p3 + p4
