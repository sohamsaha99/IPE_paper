library(dplyr)

# Discrete beta density on grid {0, 0.01, ..., 0.99, 1.00} =: {x_1, ..., x_K}
# Returns p1, ..., pK such that p_K = P(grid[k] < X <= grid[k+1])
beta_interval_probs <- function(shape1, shape2, grid = seq(0, 1, by = 0.01)) {
  if (length(grid) < 2L) stop("grid must have at least two points.")
  stopifnot(!is.unsorted(grid))
  if (!isTRUE(all.equal(min(grid), 0)) || !isTRUE(all.equal(max(grid), 1)))
    stop("grid[1] must be 0 and grid[K] must be 1.")
  if (any(diff(grid) <= 0)) stop("grid must be strictly increasing.")
  grid <- c(grid, Inf)
  cdf_vals <- pbeta(grid, shape1 = shape1, shape2 = shape2)
  p <- diff(cdf_vals)  # F(grid[k+1]) - F(grid[k])
  p
}

# Function to generate data give scenario and sample sizes
get_data <- function(n0, n1, scenario) {
  stopifnot(scenario == 3)
  # Fix grid to make a discrete problem
  grid <- seq(0, 1, by = 0.01)
  # Generate treatment arm
  trt <- data.frame(matrix(NA, nrow = n1, ncol = 2))
  colnames(trt) <- c("X", "Y")
  p_trt <- beta_interval_probs(shape1 = 6, shape2 = 4, grid = grid)
  trt$X <- sample(grid, size = n1, replace = TRUE, prob = p_trt)
  trt <- trt %>%
    mutate(Y = 1 + 2 * X + X^2 + 2 * sin(2 * pi * X)^2 + rnorm(n1, mean = 0, sd = 1/sqrt(2)))
  # Generate control arm
  con <- data.frame(matrix(NA, nrow = n0, ncol = 2))
  colnames(con) <- c("X", "Y")
  p_con <- beta_interval_probs(shape1 = 3, shape2 = 2, grid = grid)
  con$X <- sample(grid, size = n0, replace = TRUE, prob = p_con)
  con <- con %>%
    mutate(Y = 1 + 2 * X + X^2 + 2 * sin(2 * pi * X)^2 + rnorm(n0, mean = 0, sd = 1/sqrt(2)))
  # Generate summary of control arm
  target <- data.frame(matrix(nrow = 2, ncol = 3))
  colnames(target) <- c("xname", "value", "type")
  target$xname <- c("X", "X")
  target$value <- c(mean(con$X), mean(con$X^2))
  target$type <- c("first", "second")
  control_mean <- mean(con$Y)
  # Return everything in a list
  list(
    trt = trt, target = target, control_mean = control_mean
  )
}

# Function to generate necessary information of the population
get_population <- function(scenario) {
  stopifnot(scenario == 3)
  grid <- seq(0, 1, by = 0.01)
  p_trt <- beta_interval_probs(shape1 = 6, shape2 = 4, grid = grid)
  categorized <- data.frame(X = grid)
  categorized <- categorized %>%
    mutate(prop = p_trt) %>%
    mutate(mean_Y = 1 + 2 * X + X^2 + 2 * sin(2 * pi * X)^2)
  p_con <- beta_interval_probs(shape1 = 3, shape2 = 2, grid = grid)
  control_moments <- rep(NA, 2)
  control_moments[1] <- sum(grid * p_con)
  control_moments[2] <- sum(grid^2 * p_con)
  control_mean <- sum((1 + 2 * grid + grid^2 + 2 * sin(2 * pi * grid)^2) * p_con)
  list(categorized = categorized,
       control_moments = control_moments,
       control_mean = control_mean
  )
}
