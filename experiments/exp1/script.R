#!/usr/bin/env Rscript

# Driver that SOURCES replicate.R and calls run_replicate()
# with explicit method parameters for each method.
# Grids:
#   n1 in {100, 200, 400, 800, 1600}; n0 = 100 * n1
#   L_Lipschitz in L_LIP_GRID; L_Gaussian in L_GAU_GRID

source("experiments/exp1/replicate.R", chdir = TRUE)  # defines run_replicate()

SCENARIO <- 3
B        <- 800
CORES    <- max(1, parallel::detectCores() - 1)

L_LIP_GRID <- c(1, 2, 4, 8, 16)
L_GAU_GRID <- c(8, 10, 12, 14, 16)
L_RKHS_GRID <- c(4, 5, 6, 7, 8)
L_RKHSLOG_GRID <- c(4, 5, 6, 7, 8)
N1_GRID    <- c(100, 200, 400, 800, 1600)

# Optionally, define richer method params here:
# (You can add more fields as your estimators grow.)
make_lipschitz_params <- function(L) {
  list(L_Lipschitz = L)
}

make_gaussian_params <- function(L, gaussian_grid_points) {
  list(L_Gaussian = L, grid_points = gaussian_grid_points)
}

make_rkhs_params <- function(L, rkhs_sigma) {
  list(L_RKHS = L, rkhs_sigma = rkhs_sigma)
}

make_rkhs_log_params <- function(L, rkhs_log_sigma) {
  list(L_RKHS_log = L, rkhs_log_sigma = rkhs_log_sigma)
}

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# ---- run Lipschitz across the grid ----
cat("=== Running Lipschitz ===\n")
for (n1 in N1_GRID) {
  n0 <- 100 * n1
  for (Llip in L_LIP_GRID) {
    params <- make_lipschitz_params(Llip)
    cat(sprintf("Lipschitz: n1=%d n0=%d L=%g\n", n1, n0, Llip))
    run_replicate(method   = "lipschitz",
                  n0       = n0,
                  n1       = n1,
                  scenario = SCENARIO,
                  B        = B,
                  cores    = CORES,
                  params   = params)
  }
}

# ---- run Gaussian mixture across the grid ----
cat("=== Running Gaussian ===\n")
for (n1 in N1_GRID) {
  n0 <- 100 * n1
  for (Lg in L_GAU_GRID) {
    params <- make_gaussian_params(Lg, gaussian_grid_points = ppoints(50))
    cat(sprintf("Gaussian: n1=%d n0=%d L=%g\n", n1, n0, Lg))
    run_replicate(method   = "gaussian",
                  n0       = n0,
                  n1       = n1,
                  scenario = SCENARIO,
                  B        = B,
                  cores    = CORES,
                  params   = params)
  }
}

# ---- run RKHS across the grid ----
cat("=== Running RKHS ===\n")
for (n1 in N1_GRID) {
  n0 <- 100 * n1
  for (Lr in L_RKHS_GRID) {
    params <- make_rkhs_params(Lr, rkhs_sigma = 1 / 10)
    cat(sprintf("RKHS: n1=%d n0=%d L=%g\n", n1, n0, Lr))
    run_replicate(method   = "rkhs",
                  n0       = n0,
                  n1       = n1,
                  scenario = SCENARIO,
                  B        = B,
                  cores    = CORES,
                  params   = params)
  }
}

# ---- run RKHS log across the grid ----
cat("=== Running RKHS log ===\n")
for (n1 in N1_GRID) {
  n0 <- 100 * n1
  for (Lrl in L_RKHSLOG_GRID) {
    params <- make_rkhs_log_params(Lrl, rkhs_log_sigma = 1 / 5)
    cat(sprintf("RKHS log: n1=%d n0=%d L=%g\n", n1, n0, Lrl))
    run_replicate(method   = "rkhs_log",
                  n0       = n0,
                  n1       = n1,
                  scenario = SCENARIO,
                  B        = B,
                  cores    = CORES,
                  params   = params)
  }
}

cat("All done. Results saved to results/\n")

