#!/usr/bin/env Rscript

# A single entry-point to run *one* (n0, n1) replication set for a chosen method.
# - When sourced(): call run_replicate(method=..., n0=..., n1=..., params=list(...))
# - When executed with Rscript: minimal CLI (method, scenario, n0, n1, B, cores) and
#   default params per method (so you can test standalone).

suppressPackageStartupMessages({
  library(parallel)
  library(dplyr)
})

# ---- helpers ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

get_num <- function(kv, key, default = NA_real_) {
  if (!is.null(kv[[key]])) as.numeric(kv[[key]]) else default
}
parse_cli <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) == 0) return(list())
  kv <- strsplit(args, "=", fixed = TRUE)
  setNames(lapply(kv, `[[`, 2), sapply(kv, `[[`, 1))
}

# ---- core runner ----
run_replicate <- function(method,
                          n0, n1,
                          scenario,
                          B = 800,
                          cores = max(1, parallel::detectCores() - 1),
                          params = list()) {
  stopifnot(is.numeric(n0), is.numeric(n1), n0 > 1, n1 > 1)
  dir.create("results", showWarnings = FALSE, recursive = TRUE)

  # Common sources
  source(paste0("./generate/scenario", scenario, ".R"), chdir = TRUE, local = TRUE)

  set.seed(100)
  ncores <- max(1, as.integer(cores))

  if (identical(tolower(method), "lipschitz")) {
    # ---- Lipschitz method ----
    source("./methods/lipschitz.R", chdir = TRUE, local = TRUE)

    L_Lipschitz <- params$L_Lipschitz
    if (is.null(L_Lipschitz)) stop("params$L_Lipschitz must be provided for method='lipschitz'.")

    population <- get_population(scenario = scenario)
    truth_obj  <- true_theta_min_lipschitz(population, L = L_Lipschitz)
    if (!truth_obj$is_valid_solution) stop("Regularization is infeasible for population.")
    theta_true <- truth_obj$theta_min

    run_one_rep <- function() {
      dat <- get_data(n0 = n0, n1 = n1, scenario = scenario)
      est <- estimate_theta_min_lipschitz(dat, L = L_Lipschitz)
      est$theta_min
    }

    theta_hat <- unlist(mclapply(seq_len(B), function(b) run_one_rep(), mc.cores = ncores), use.names = FALSE)

    results <- data.frame(rep          = seq_len(B),
                          method       = "lipschitz",
                          theta_hat    = theta_hat,
                          theta_true   = theta_true,
                          n0           = n0,
                          n1           = n1,
                          scaled_error = sqrt(n1) * (theta_hat - theta_true),
                          L_Lipschitz  = L_Lipschitz,
                          scenario     = scenario,
                          stringsAsFactors = FALSE
    )

    attr(results, "params") <- params
    attr(results, "population_result") <- truth_obj
    outfile <- sprintf("results/lipschitz_scen%d_n1_%d_n0_%d_L_%g.rds",
                       scenario, n1, n0, L_Lipschitz)
    saveRDS(results, outfile)
    message("Saved: ", outfile)
    return(invisible(results))
  }

  if (identical(tolower(method), "gaussian")) {
    # ---- Gaussian mixture method ----
    source("./methods/weighted_gaussian.R", chdir = TRUE, local = TRUE)

    L_Gaussian <- params$L_Gaussian
    if (is.null(L_Gaussian)) stop("params$L_Gaussian must be provided for method='gaussian'.")

    gaussian_grid_points <- params$grid_points
    if (is.null(gaussian_grid_points)) stop("params$grid_points must be provided for method='gaussian'.")

    population <- get_population(scenario = scenario)
    truth_obj  <- true_theta_min_gaussian(population,
                                          grid_points = gaussian_grid_points,
                                          L = L_Gaussian)
    if (!truth_obj$is_valid_solution) stop("Regularization is infeasible for population.")
    theta_true <- truth_obj$theta_min

    run_one_rep <- function() {
      dat <- get_data(n0 = n0, n1 = n1, scenario = scenario)
      est <- estimate_theta_min_gaussian(dat,
                                         grid_points = gaussian_grid_points,
                                         L = L_Gaussian)
      est$theta_min
    }

    theta_hat <- unlist(mclapply(seq_len(B), function(b) run_one_rep(), mc.cores = ncores), use.names = FALSE)

    results <- data.frame(
                          rep          = seq_len(B),
                          method       = "gaussian",
                          theta_hat    = theta_hat,
                          theta_true   = theta_true,
                          n0           = n0,
                          n1           = n1,
                          scaled_error = sqrt(n1) * (theta_hat - theta_true),
                          L_Gaussian   = L_Gaussian,
                          scenario     = scenario,
                          stringsAsFactors = FALSE
    )

    attr(results, "params") <- params
    attr(results, "population_result") <- truth_obj
    outfile <- sprintf("results/gaussian_scen%d_n1_%d_n0_%d_L_%g.rds",
                       scenario, n1, n0, L_Gaussian)
    saveRDS(results, outfile)
    message("Saved: ", outfile)
    return(invisible(results))
  }

  if (identical(tolower(method), "rkhs")) {
    # ---- w \in RKHS with Gaussian kernel method ----
    source("./methods/RKHS.R", chdir = TRUE, local = TRUE)

    L_RKHS <- params$L_RKHS
    if (is.null(L_RKHS)) stop("params$L_RKHS must be provided for method='rkhs'.")

    rkhs_sigma <- params$rkhs_sigma
    if (is.null(rkhs_sigma)) stop("params$rkhs_sigma must be provided for method='rkhs'.")

    population <- get_population(scenario = scenario)
    truth_obj  <- true_theta_min_rkhs(population,
                                      L = L_RKHS,
                                      sigm = rkhs_sigma
    )
    if (!truth_obj$is_valid_solution) stop("Regularization is infeasible for population.")
    theta_true <- truth_obj$theta_min

    run_one_rep <- function() {
      dat <- get_data(n0 = n0, n1 = n1, scenario = scenario)
      est <- estimate_theta_min_rkhs(dat,
                                     L = L_RKHS,
                                     sigm = rkhs_sigma
      )
      est$theta_min
    }

    theta_hat <- unlist(mclapply(seq_len(B), function(b) run_one_rep(), mc.cores = ncores), use.names = FALSE)

    results <- data.frame(
                          rep          = seq_len(B),
                          method       = "rkhs",
                          theta_hat    = theta_hat,
                          theta_true   = theta_true,
                          n0           = n0,
                          n1           = n1,
                          scaled_error = sqrt(n1) * (theta_hat - theta_true),
                          L_RKHS   = L_RKHS,
                          scenario     = scenario,
                          stringsAsFactors = FALSE
    )

    attr(results, "params") <- params
    attr(results, "population_result") <- truth_obj
    outfile <- sprintf("results/rkhs_scen%d_n1_%d_n0_%d_L_%g.rds",
                       scenario, n1, n0, L_RKHS)
    saveRDS(results, outfile)
    message("Saved: ", outfile)
    return(invisible(results))
  }

  if (identical(tolower(method), "rkhs_log")) {
    # ---- log w \in RKHS with Gaussian kernel method ----
    source("./methods/RKHS_log.R", chdir = TRUE, local = TRUE)

    L_RKHS_log <- params$L_RKHS_log
    if (is.null(L_RKHS_log)) stop("params$L_RKHS_log must be provided for method='rkhs_log'.")

    rkhs_log_sigma <- params$rkhs_log_sigma
    if (is.null(rkhs_log_sigma)) stop("params$rkhs_log_sigma must be provided for method='rkhs_log'.")

    population <- get_population(scenario = scenario)
    truth_obj  <- true_theta_min_rkhs_log(population,
                                          L = L_RKHS_log,
                                          sigm = rkhs_log_sigma
    )
    if (!truth_obj$is_valid_solution) stop("Regularization is infeasible for population.")
    theta_true <- truth_obj$theta_min

    run_one_rep <- function() {
      dat <- get_data(n0 = n0, n1 = n1, scenario = scenario)
      est <- estimate_theta_min_rkhs_log(dat,
                                         L = L_RKHS_log,
                                         sigm = rkhs_log_sigma
      )
      est$theta_min
    }

    theta_hat <- unlist(mclapply(seq_len(B), function(b) run_one_rep(), mc.cores = ncores), use.names = FALSE)

    results <- data.frame(
                          rep          = seq_len(B),
                          method       = "rkhs_log",
                          theta_hat    = theta_hat,
                          theta_true   = theta_true,
                          n0           = n0,
                          n1           = n1,
                          scaled_error = sqrt(n1) * (theta_hat - theta_true),
                          L_RKHS_log   = L_RKHS_log,
                          scenario     = scenario,
                          stringsAsFactors = FALSE
    )

    attr(results, "params") <- params
    attr(results, "population_result") <- truth_obj
    outfile <- sprintf("results/rkhs_log_scen%d_n1_%d_n0_%d_L_%g.rds",
                       scenario, n1, n0, L_RKHS_log)
    saveRDS(results, outfile)
    message("Saved: ", outfile)
    return(invisible(results))
  }
  stop("Unknown method: '", method, "'. Use 'lipschitz' or 'gaussian' or 'rkhs' or 'rkhs_log'.")
}

# ---- if called from Rscript, parse minimal CLI and run with default params ----
if (!interactive() && sys.nframe() == 0) {
  kv <- parse_cli()
  method  <- kv[["method"]]
  n0      <- get_num(kv, "n0")
  n1      <- get_num(kv, "n1")
  scenario<- get_num(kv, "scenario", 1)
  B       <- get_num(kv, "B", 800)
  cores   <- get_num(kv, "cores", parallel::detectCores() - 1)

  if (is.null(method) || is.na(n0) || is.na(n1)) {
    stop("Usage: Rscript replicate.R method=<lipschitz|gaussian> n0=<int> n1=<int> [scenario=1 B=800 cores=...]")
  }

  # Minimal default params for standalone tests
  if (tolower(method) == "lipschitz") {
    params <- list(L_Lipschitz = 8)
  } else if (tolower(method) == "gaussian") {
    params <- list(L_Gaussian = 8, grid_K = 50)
  } else {
    stop("Unknown method: ", method)
  }

  run_replicate(method = method, n0 = n0, n1 = n1,
                scenario = scenario, B = B, cores = cores, params = params)
}

