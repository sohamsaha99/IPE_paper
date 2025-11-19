#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(purrr)
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
                      # ensure data.frame
                      as_tibble(obj)
})

  df <- bind_rows(df_list) %>%
    mutate(method = tolower(method),
           # unify L across methods
           L = coalesce(.data$L_Lipschitz, .data$L_Gaussian),
           L = as.numeric(L),
           scenario = as.integer(scenario)
    )

  # basic sanity checks
  stopifnot(all(c("theta_hat", "theta_true", "scaled_error", "n1", "method", "L") %in% names(df)))
  if (any(is.na(df$L))) {
    stop("Some rows have missing L. Check that results contain L_Lipschitz or L_Gaussian.")
  }

  message("Loaded ", scales::comma(nrow(df)), " rows from ", length(files), " files.")

  dir.create("figures", showWarnings = FALSE, recursive = TRUE)

  # ---- Figure 1: Truth moves with L; estimator tracks it ----
  # Aggregate to show mean theta_hat with 95% CI vs L, alongside theta_true
  summ1 <- df %>%
    group_by(method, scenario, L) %>%
    summarise(
              theta_true = dplyr::first(theta_true), # constant within (method, scenario, L)
              mean_hat   = mean(theta_hat),
              # se_hat     = sd(theta_hat) / sqrt(dplyr::n()),
              se_hat     = sd(theta_hat),
              n_reps     = dplyr::n(),
              .groups = "drop"
              ) %>%
    mutate(
           lower = mean_hat - 1.96 * se_hat,
           upper = mean_hat + 1.96 * se_hat
    )

    p1 <- ggplot(summ1, aes(x = L)) +
      geom_line(aes(y = theta_true), linetype = "dashed") +
      geom_point(aes(y = mean_hat)) +
      geom_line(aes(y = mean_hat)) +
      geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
      facet_wrap(~ method + paste0("scenario=", scenario), scales = "free_y") +
      labs(
           title = expression(paste("Estimator tracks truth as regularization ", L, " varies")),
           subtitle = "Dashed: true θ(L). Solid & points: mean estimate with 95% CI across replications.",
           x = "Regularization L",
           y = expression(theta),
           caption = "Each panel: method × scenario. CIs from replication variability."
           ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())

    ggsave("figures/theta_vs_L.pdf", plot = p1, width = 8.5, height = 5.0)
    ggsave("figures/theta_vs_L.png", plot = p1, width = 8.5, height = 5.0, dpi = 200)

    # ---- Figure 2: Scaled error stability across sample sizes ----
    # Boxplots of scaled_error by n1, faceted by method and L
    df$n1_factor <- factor(df$n1, levels = sort(unique(df$n1)))

    p2 <- ggplot(df, aes(x = n1_factor, y = scaled_error)) +
      geom_hline(yintercept = 0, linewidth = 0.3) +
      geom_boxplot(outlier.alpha = 0.25, width = 0.7) +
      facet_grid(method ~ paste0("L=", L), scales = "free_y") +
      labs(
           title = expression(paste("Scaled error ", sqrt(n[1]) * " ("*hat(theta) - theta[true]*") remains stable across sample sizes")),
           x = expression(n[1]),
           y = expression(sqrt(n[1]) * (hat(theta) - theta[true])),
           caption = "Each column: regularization level L. Rows: method."
           ) +
      theme_bw() +
      theme(panel.grid.minor = element_blank())

    ggsave("figures/scaled_error_stability.pdf", plot = p2, width = 10, height = 6.5)
    ggsave("figures/scaled_error_stability.png", plot = p2, width = 10, height = 6.5, dpi = 200)

    # ---- Optional numeric check: does the spread of scaled_error change with n1? ----
    spread <- df %>%
      group_by(method, L, n1) %>%
      summarise(
                mean_scaled = mean(scaled_error),
                sd_scaled   = sd(scaled_error),
                iqr_scaled  = IQR(scaled_error),
                .groups = "drop"
      )

      write_csv(spread, "figures/scaled_error_spread_by_n1.csv")

      message("Wrote figures to figures/: ",
              "\n - theta_vs_L.pdf/png",
              "\n - scaled_error_stability.pdf/png",
              "\nAlso wrote: figures/scaled_error_spread_by_n1.csv")
}

if (!interactive() && sys.nframe() == 0) main()

