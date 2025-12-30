#' In this experiment, we look at the performance of MAIC estimator under
#' different scenarios
rm(list = ls())
set.seed(20240620)
#### STEP 1: Include all necessary source files ####
library(dplyr)
library(ggplot2)
library(viridis)
source("methods/MAIC.R", chdir = TRUE, local = TRUE)


n_trt <- 200
n_con <- n_trt

#### STEP 3: REPLICATE MAIC ESTIMATE ON SYNTHETIC DATASET ####
B <- 500
# Scenario 1
source("generate/scenario1.R",chdir = TRUE, local = TRUE)
scenario1_results <- rep(NA, B)
for (b in 1:B) {
  dat_synth <- get_data(n0 = n_con, n1 = n_trt, scenario = 1)
  estimate <- MAIC_estimate(dat_synth)
  scenario1_results[b] <- estimate
}
if (TRUE) {
  # Plot densities and conditional mean and sample scatter plot
  xpoints <- seq(from = 0, to = 1, length.out = 1000)
  df_densities <- data.frame(matrix(nrow = 2*length(xpoints), ncol = 0))
  df_densities$x <- c(xpoints, xpoints)
  df_densities$Density <- c(dbeta(xpoints, 3, 2),
                            dbeta(xpoints, 6, 4))
  df_densities$Group <- c(rep("Control", length(xpoints)),
                          rep("Treatment", length(xpoints)))
  df_densities$mu_Y <- 1 + 2 * df_densities$x + df_densities$x^2
  p <- ggplot(df_densities, aes(x = x, y = Density, fill = Group, color = Group)) +
    geom_area(alpha = 0.5, position = "identity") +
    scale_fill_viridis_d() + scale_color_viridis_d() + theme_bw(base_size = 18)
  ggsave("results/experiment0/Scenario1_a.png", plot = p,
         width = 5, height = 4)

  df_sample <- data.frame(matrix(nrow = n_trt + n_con, ncol = 0))
  df_sample$x <- c(rbeta(n_con, 3, 2), rbeta(n_trt, 6, 4))
  df_sample$Group <- c(rep("Control", n_con), rep("Treatment", n_trt))
  df_sample$Y <- 1 + 2 * df_sample$x + df_sample$x^2 + rnorm(n_trt + n_con, 0, 1/sqrt(2))
  p <- ggplot(df_densities, aes(x = x, y = mu_Y)) +
    geom_line() +
    geom_point(data = df_sample, aes(x = x, y = Y, color = Group, shape = Group),
               alpha = 0.7, size = 1) +
    scale_color_viridis_d() + theme_bw(base_size = 18) +
    labs(y = "Y")
  ggsave("results/experiment0/Scenario1_b.png", plot = p,
         width = 5, height = 4)
  # Save legend
  p <- p + theme_bw(base_size = 10)
  ggsave("results/experiment0/legend.png", plot = p,
         width = 5, height = 4)
}

# Scenario 2
source("generate/scenario2.R",chdir = TRUE, local = TRUE)
scenario2_results <- rep(NA, B)
for (b in 1:B) {
  dat_synth <- get_data(n0 = n_con, n1 = n_trt, scenario = 2)
  estimate <- MAIC_estimate(dat_synth)
  scenario2_results[b] <- estimate
}
if (TRUE) {
  # Plot densities and conditional mean and sample scatter plot
  xpoints <- seq(from = 0, to = 1.2, length.out = 1000)
  df_densities <- data.frame(matrix(nrow = 2*length(xpoints), ncol = 0))
  df_densities$x <- c(xpoints, xpoints)
  df_densities$Density <- c(dnorm(xpoints, 0.6, 0.2),
                            dnorm(xpoints, 0.6, 0.15))
  df_densities$Group <- c(rep("Control", length(xpoints)),
                          rep("Treatment", length(xpoints)))
  df_densities$mu_Y <- 1 + 2 * df_densities$x + df_densities$x^2 +
    2 * sin(2 * pi * df_densities$x)^2
  p <- ggplot(df_densities, aes(x = x, y = Density, fill = Group, color = Group)) +
    geom_area(alpha = 0.5, position = "identity") +
    scale_fill_viridis_d() + scale_color_viridis_d() + theme_bw(base_size = 18)
  ggsave("results/experiment0/Scenario2_a.png", plot = p,
         width = 5, height = 4)

  df_sample <- data.frame(matrix(nrow = n_trt + n_con, ncol = 0))
  df_sample$x <- c(rnorm(n_con, 0.6, 0.2), rnorm(n_trt, 0.6, 0.15))
  df_sample$Group <- c(rep("Control", n_con), rep("Treatment", n_trt))
  df_sample$Y <- 1 + 2 * df_sample$x + df_sample$x^2 +
    2 * sin(2 * pi * df_sample$x)^2 +
    rnorm(n_trt + n_con, 0, 1/sqrt(2))
  p <- ggplot(df_densities, aes(x = x, y = mu_Y)) +
    geom_line() +
    geom_point(data = df_sample, aes(x = x, y = Y, color = Group, shape = Group),
               alpha = 0.7, size = 1) +
    scale_color_viridis_d() + theme_bw(base_size = 18) +
    labs(y = "Y")
  ggsave("results/experiment0/Scenario2_b.png", plot = p,
         width = 5, height = 4)
  # Save legend
  p <- p + theme_bw(base_size = 10)
  ggsave("results/experiment0/legend.png", plot = p,
         width = 5, height = 4)
}

# Scenario 3
source("generate/scenario3.R",chdir = TRUE, local = TRUE)
scenario3_results <- rep(NA, B)
for (b in 1:B) {
  dat_synth <- get_data(n0 = n_con, n1 = n_trt, scenario = 3)
  estimate <- MAIC_estimate(dat_synth)
  scenario3_results[b] <- estimate
}
if (TRUE) {
  # Plot densities and conditional mean and sample scatter plot
  xpoints <- seq(from = 0, to = 1.0, length.out = 1000)
  df_densities <- data.frame(matrix(nrow = 2*length(xpoints), ncol = 0))
  df_densities$x <- c(xpoints, xpoints)
  df_densities$Density <- c(dbeta(xpoints, 3, 2),
                            dbeta(xpoints, 6, 4))
  df_densities$Group <- c(rep("Control", length(xpoints)),
                          rep("Treatment", length(xpoints)))
  df_densities$mu_Y <- 1 + 2 * df_densities$x + df_densities$x^2 +
    2 * sin(2 * pi * df_densities$x)^2
  p <- ggplot(df_densities, aes(x = x, y = Density, fill = Group, color = Group)) +
    geom_area(alpha = 0.5, position = "identity") +
    scale_fill_viridis_d() + scale_color_viridis_d() + theme_bw(base_size = 18)
  ggsave("results/experiment0/Scenario3_a.png", plot = p,
         width = 5, height = 4)

  df_sample <- data.frame(matrix(nrow = n_trt + n_con, ncol = 0))
  df_sample$x <- c(rbeta(n_con, 3, 2), rbeta(n_trt, 6, 4))
  df_sample$Group <- c(rep("Control", n_con), rep("Treatment", n_trt))
  df_sample$Y <- 1 + 2 * df_sample$x + df_sample$x^2 +
    2 * sin(2 * pi * df_sample$x)^2 +
    rnorm(n_trt + n_con, 0, 1/sqrt(2))
  p <- ggplot(df_densities, aes(x = x, y = mu_Y)) +
    geom_line() +
    geom_point(data = df_sample, aes(x = x, y = Y, color = Group, shape = Group),
               alpha = 0.7, size = 1) +
    scale_color_viridis_d() + theme_bw(base_size = 18) +
    labs(y = "Y")
  ggsave("results/experiment0/Scenario3_b.png", plot = p,
         width = 5, height = 4)
  # Save legend
  p <- p + theme_bw(base_size = 10)
  ggsave("results/experiment0/legend.png", plot = p,
         width = 5, height = 4)
}

#### STEP 4: COMBINE ALL RESULTS AND PLOT ####
results <- data.frame(matrix(nrow = 3*B, ncol = 2))
colnames(results) <- c("Estimate", "Scenario")
results$Estimate <- c(scenario1_results, scenario2_results, scenario3_results)
results$Scenario <- paste("Scenario", rep(1:3, each = B))
p <- ggplot(results, aes(x = Scenario, y = Estimate)) +
  geom_boxplot() +
  theme_bw(base_size = 18) +
  labs(y  = TeX(r'($\hat{\theta}_{MAIC}$)'), x = "") +
  ylim(c(-0.4, 0.4)) +
  geom_abline(slope = 0, intercept = 0, linetype = "dashed", color = "red")
ggsave("results/experiment0/boxplot.png", plot = p,
       width = 5, height = 4)
