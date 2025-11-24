#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggpubr)
})

# ---- helpers ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

must_exist <- function(path) {
  if (!dir.exists(path)) stop("Directory not found: ", path)
}

# ---- load all results ----
main <- function() {
  must_exist("results")
  files <- list.files("results", pattern = "\\.rds$", full.names = TRUE)
  if (length(files) == 0) {
    stop("No .rds results found in results/. Run script.R first to generate outputs.")
  }

  df_list <- lapply(files, function(f) {
                      obj <- readRDS(f)
                      # ensure tibble
                      as_tibble(obj)
})

  df_raw <- bind_rows(df_list) %>%
    mutate(
           method   = tolower(method),
           scenario = as.integer(scenario)
    )
  df_raw <- df_raw %>% filter(scenario == 3)

  # ---- unify L across methods (supports lipschitz, gaussian, rkhs, rkhs_log) ----
  # Accept any of the L_* columns; coalesce in a priority-free way.
  L_cols <- c("L_Lipschitz", "L_Gaussian", "L_RKHS", "L_RKHS_log")
  present_L_cols <- intersect(L_cols, names(df_raw))

  df <- df_raw %>%
    mutate(
           L = coalesce(!!!rlang::syms(present_L_cols)),
           L = suppressWarnings(as.numeric(L))
    )

  # ---- drop NA/NULL theta_hat rows and recompute scaled_error safely ----
  df <- df %>%
    filter(!is.na(theta_hat)) %>%               # ignore NA estimates
    filter(is.finite(theta_hat)) %>%            # guard against Inf
    filter(!is.na(theta_true)) %>%              # ensure truth present
    mutate(scaled_error = sqrt(n1) * (theta_hat - theta_true))

  # basic sanity checks
  needed <- c("theta_hat", "theta_true", "scaled_error", "n1", "method", "L")
  if (!all(needed %in% names(df))) {
    missing <- paste(setdiff(needed, names(df)), collapse = ", ")
    stop("Missing required columns: ", missing)
  }
  if (any(is.na(df$L))) {
    stop("Some rows have missing L. Make sure results contain one of: ",
         paste(L_cols, collapse = ", "), ".")
  }

  message("Loaded ", scales::comma(nrow(df)), " rows from ", length(files), " files.")

  dir.create("figures", showWarnings = FALSE, recursive = TRUE)

  # ---- Figure 1: Truth vs L with quantile 95% CI ----
  # Use quantiles across replicates instead of sd-based CI.
  summ1 <- df %>%
    group_by(method, scenario, L) %>%
    summarise(
              theta_true = dplyr::first(theta_true[!is.na(theta_true)]),
              mean_hat   = mean(theta_hat, na.rm = TRUE),
              median_hat = median(theta_hat, na.rm = TRUE),
              lower      = quantile(theta_hat, probs = 0.025, na.rm = TRUE, names = FALSE),
              upper      = quantile(theta_hat, probs = 0.975, na.rm = TRUE, names = FALSE),
              n_reps     = dplyr::n(),
              .groups    = "drop"
    )

    p1 <- ggplot(summ1, aes(x = L)) +
      geom_line(aes(y = theta_true), linetype = "dashed") +
      geom_point(aes(y = median_hat)) +
      geom_line(aes(y = median_hat)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      facet_wrap(~ method + paste0("scenario=", scenario), scales = "free") +
      labs(
           title = expression(paste("Estimator tracks truth as regularization ", L, " varies")),
           subtitle = "Dashed: true θ(L). Solid/points: median estimate; band: 2.5%–97.5% quantiles across replications.",
           x = "Regularization L",
           y = expression(theta),
           caption = "Each panel: method × scenario. Intervals are empirical quantiles over replicates."
           ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())

    ggsave("figures/theta_vs_L.pdf", plot = p1, width = 8.5, height = 5.0)
    ggsave("figures/theta_vs_L.png", plot = p1, width = 8.5, height = 5.0, dpi = 200)

    # ---- Figure 2: Scaled error stability across sample sizes ----
    df$n1_factor <- factor(df$n1, levels = sort(unique(df$n1)))

    p2 <- ggplot(df, aes(x = n1_factor, y = scaled_error)) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_boxplot(outlier.alpha = 0.25, width = 0.7) +
      # facet_grid(method ~ paste0("L=", L), scales = "free_y") +
      facet_wrap(method ~ paste0("L=", L), scales = "free_y") +
      labs(
           title = expression(paste("Scaled error ", sqrt(n[1]) * " ("*hat(theta) - theta[true]*") across sample sizes")),
           x = expression(n[1]),
           y = expression(sqrt(n[1]) * (hat(theta) - theta[true])),
           caption = "Each column: regularization level L. Rows: method. Boxes omit NA/NULL estimates."
           ) +
      theme_pubr() +
      theme(panel.grid.minor = element_blank())

    ggsave("figures/scaled_error_stability.pdf", plot = p2, width = 10, height = 6.5)
    ggsave("figures/scaled_error_stability.png", plot = p2, width = 10, height = 6.5, dpi = 200)

    # ---- Optional numeric check: spread of scaled_error by n1 ----
    spread <- df %>%
      group_by(method, L, n1) %>%
      summarise(
                mean_scaled = mean(scaled_error, na.rm = TRUE),
                sd_scaled   = sd(scaled_error, na.rm = TRUE),
                iqr_scaled  = IQR(scaled_error, na.rm = TRUE),
                .groups     = "drop"
      )

      write_csv(spread, "figures/scaled_error_spread_by_n1.csv")

      message(
              "Wrote figures to figures/: ",
              "\n - theta_vs_L.pdf/png",
              "\n - scaled_error_stability.pdf/png",
              "\nAlso wrote: figures/scaled_error_spread_by_n1.csv"
      )
}

if (!interactive() && sys.nframe() == 0) main()

