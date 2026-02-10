# Show that the scaled errors remain stable for Gaussian mixture regularization

rm(list = ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(latex2exp)
  library(readr)
  library(stringr)
  library(purrr)
  library(ggpubr)
  library(patchwork)
  library(svglite)
})

# ---- helpers ----
`%||%` <- function(a, b) if (!is.null(a)) a else b

must_exist <- function(path) {
  if (!dir.exists(path)) stop("Directory not found: ", path)
}
method_name <- "rkhs_log"
L_name <- "L_RKHS_log"
estimand_name <- "\\theta_{min}^{(ARK)}"
estimator_name <- "\\hat{\\theta}_{min}^{(ARK)}"
# ---- load all results ----
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
df_raw <- df_raw %>%
  filter(scenario %in% c(1, 2, 3)) %>%
  filter(method == method_name) %>%
  mutate(L = !!sym(L_name)) %>%
  select(rep, theta_hat, theta_true, n0, n1, scaled_error, scenario, L, method)

# ---- drop NA/NULL theta_hat rows and recompute scaled_error safely ----
df <- df_raw %>%
  filter(!is.na(theta_hat)) %>%               # ignore NA estimates
  filter(is.finite(theta_hat)) %>%            # guard against Inf
  filter(!is.na(theta_true)) %>%              # ensure truth present
  mutate(scaled_error = sqrt(n1) * (theta_hat - theta_true))

cat("Dimension of df: ", dim(df), "\n")

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

# ---- Figure 1: Truth vs L ----
summ1 <- df %>%
  group_by(n1, n0, scenario, L) %>%
  summarise(theta_true = dplyr::first(theta_true[!is.na(theta_true)]),
            mean_hat   = mean(theta_hat, na.rm = TRUE),
            median_hat = median(theta_hat, na.rm = TRUE),
            lower      = quantile(theta_hat, probs = 0.025, na.rm = TRUE, names = FALSE),
            upper      = quantile(theta_hat, probs = 0.975, na.rm = TRUE, names = FALSE),
            n_reps     = dplyr::n(),
            .groups    = "drop")# %>%
  # filter(n1 == 100)

  p1 <- ggplot(summ1, aes(x = L)) +
    geom_errorbar(aes(ymin = lower, ymax = upper, color = factor(n1), group = factor(n1)),
                  width = 0.7,
                  position = position_dodge(width = 0.7)) +
    geom_point(aes(y = median_hat, color = factor(n1), group = factor(n1)),
               size = 1, position = position_dodge(width = 0.7)) +
    geom_line(aes(y = theta_true), linetype = "dashed", color = "black", size = 0.4) +
    geom_point(aes(y = theta_true), color = "black", size = 1) +
    facet_grid(
               # rows = vars(n1),
               cols = vars(scenario),
               labeller = labeller(scenario = function(x) paste0("Scenario ", x)),
               scales = "free"
               ) +
    # labs(title = TeX(sprintf("Estimand $%s$ and estimator $%s$ as regularization $L$ varies", estimand_name, estimator_name)),
    labs(x = "Regularization L",
         y = TeX(sprintf("$%s$ and $%s$", estimand_name, estimator_name)),
         # caption = "Each panel: method × scenario. Intervals are empirical quantiles over replicates."
         ) +
    theme_bw() +
    scale_x_continuous(breaks = sort(unique(summ1$L))) +
    scale_color_discrete(name = TeX("$n_1$")) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "top",
          legend.margin = margin(t = -5, b = -5),
          plot.title = element_text(size = 10)
    )
  print(p1)
  ggsave(paste0("figures/theta_vs_L", method_name, ".svg"), plot = p1, width = 8.5, height = 2.5)
  ggsave("figures/theta_vs_L.pdf", plot = p1, width = 8.5, height = 5.0)
  ggsave("figures/theta_vs_L.png", plot = p1, width = 8.5, height = 5.0, dpi = 200)

  # ---- Figure 2: Scaled error stability across sample sizes ----
  p2 <- ggplot(df, aes(x = factor(n1), y = scaled_error)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_boxplot(outlier.alpha = 0.25, width = 0.7) +
    facet_grid(rows = vars(scenario), cols = vars(L),
               labeller = labeller(
                                   scenario = function(x) paste0("Scenario ", x),
                                   L  = function(x) paste0("L=", x)
               ),
               scales = "free_y") +
    labs(
         title = TeX(sprintf("Scaled error $\\sqrt{n_1} \\left(%s - %s\\right)$ across sample sizes", estimator_name, estimand_name)),
         x = expression(n[1]),
         y = TeX(sprintf("$ \\sqrt{n_1} \\left(%s - %s\\right) $", estimator_name, estimand_name)),
         # caption = "Columns: regularization L. Rows: Scenario."
         ) +
    ylim(c(-3, 3)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(size = 10)
    )
  print(p2)
  ggsave("figures/scaled_error_stability.pdf", plot = p2, width = 10, height = 6.5)
  ggsave("figures/scaled_error_stability.png", plot = p2, width = 10, height = 6.5, dpi = 200)

  # ---- Figure 3: Optimizer distribution for different L ----
  # Get all .rds files for the corresponding method
  get_file_name <- function(method_name, n0, n1, scenario, L) {
    sprintf("results/%s_scen%d_n1_%d_n0_%d_L_%g.rds",
                       method_name, scenario, n1, n0, L)
  }
  dens_list <- list()
  for (i in seq_len(nrow(summ1))) {
    filename <- get_file_name(method_name, summ1$n0[i], summ1$n1[i], summ1$scenario[i], summ1$L[i])
    obj <- read_rds(file = filename)
    xvalues <- attr(obj, "population_result")$X
    f0_min <- attr(obj, "population_result")$wts_min * attr(obj, "population_result")$prop
    dens_list[[i]] <- data.frame(
                                 n0 = summ1$n0[i],
                                 n1 = summ1$n1[i],
                                 scenario = summ1$scenario[i],
                                 L = summ1$L[i],
                                 X = xvalues,
                                 f0_min = f0_min
    )
  }
  dens_df <- bind_rows(dens_list)
  dens_df <- dens_df %>% filter(n1 == min(n1))
  p3 <- ggplot(dens_df, aes(x = X, y = f0_min)) +
    geom_area(alpha = 0.6) +
    facet_grid(rows = vars(scenario), cols = vars(L),
               labeller = labeller(
                                   scenario = function(x) paste0("Scenario ", x),
                                   L  = function(x) paste0("L=", x)
               ),
               scales = "free_y") +

    labs(
         title = TeX("$f_{0,min} = \\arg\\min_{f\\in \u2131} \\tilde{\\theta}(f)$ as regularization $L$ varies"),
         x = TeX("$x$"),
         y = TeX("$f_{0,min}(x)$")
         ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          plot.title = element_text(size = 10)
    )
  print(p3)
  ggsave("figures/f0_min.pdf", plot = p3, width = 10, height = 6.5)
  ggsave("figures/f0_min.png", plot = p3, width = 10, height = 6.5, dpi = 200)

  # ---- Combine the three plots vertically ----
  p <- (p1 / p2 / p3) +
    plot_layout(heights = c(1, 2, 1.5))
  # print(p)
  # Save as a single-page PDF (US Letter); adjust size if you prefer A4, etc.
  ggsave(sprintf("figures/combined_results_%s.pdf", method_name),
         plot   = p,
         width  = 8.5,  # inches
         height = 11,   # inches
         dpi    = 300)
  ggsave(sprintf("figures/combined_results_%s.png", method_name),
         plot   = p,
         width  = 8.5,  # inches
         height = 11,   # inches
         dpi    = 300)
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


