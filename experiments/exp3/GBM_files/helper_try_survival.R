rm(list = ls())
for (Lvalue_global in c(0.5, 1, 2, 4, 8)) {
  print("============================")
  print(Lvalue_global)
  cat(paste("Starting: ", Sys.time()))
  source("try_survival.R", chdir = TRUE, local = TRUE)
  cat(paste("Finished: ", Sys.time()))
}

stopifnot(FALSE)
L_vec <- c(0.5, 1, 2, 4, 8)
dat_list <- list()
k <- 0
for (L in L_vec) {
  k <- k + 1
  dat_list[[k]] <- readRDS(paste0("results/L_value_", L, ".RDS"))$result
  dat_list[[k]]$L <- L
}
# Make one dataframe
library(dplyr)
results <- bind_rows(dat_list)
# Dlete results outside [-15, 15], probably computation error
results <- results %>%
  filter(estimate < 15, estimate > -15) %>%
  # For MAIC, set L = 0
  mutate(L = ifelse(method == "maic", 0, L))


p1 <- ggplot(results, aes(x = factor(L), y = estimate, fill = method)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "#888888") +
  scale_fill_manual(name = "Estimate",
    labels = c(TeX(r'($\hat{\theta}_{MAIC}$)'), TeX(r'($\hat{\theta}_{max}$)'), TeX(r'($\hat{\theta}_{min}$)')),
    values = c("#B79F00", "#F8766D", "#00BFC4")) +
  # scale_color_discrete("Estimate", labels = c(TeX(r'($\hat{\theta}_{max}$)'), TeX(r'($\hat{\theta}_{min}$)'))) +
  # geom_line(aes(y = truth, group = type, color = type), lwd=0.5, lty = "longdash") +
  # geom_point(aes(y = truth, group = type, shape = type, color = type), size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  # scale_shape_discrete("True Values", labels = c(TeX(r'(${\theta}_{max}$)'), TeX(r'(${\theta}_{min}$)'))) +
  # scale_color_discrete("True Values", labels = c(TeX(r'(${\theta}_{max}$)'), TeX(r'(${\theta}_{min}$)'))) +
  # facet_grid(~L, labeller = labeller(L = facet_labels)) +
  labs(x = "L (Lipschitz constant)", y = TeX(r'($\hat{\theta}$: Estimated DRMST (months))')) +
  # labs(subtitle = "True and estimated values under different Lipschitz constants") +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.direction = "horizontal") +
  guides(fill = guide_legend(override.aes = list(shape = NA)),
         shape = guide_legend(override.aes = list(fill = NA))) +
  ylim(c(-6, 6))

p1
p2 <- ggplot(results %>% filter(method == "maic"), aes(y = estimate)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, color = "#888888") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(y = TeX(r'($\hat{\theta}$: Estimated DRMST)')) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(override.aes = list(shape = NA)),
         shape = guide_legend(override.aes = list(fill = NA))) +
  ylim(c(-6, 6))
p2

L_vec <- c(0.5, 1, 2, 4, 8)
list_boxplot <- vector("list", length = length(L_vec))

k <- 0
for (L in L_vec) {
  k <- k + 1
  results <- read.csv(paste0("results/experiment3/", L, "_result.csv"), header = TRUE)
  truevar_min <- var(results[, 1])
  truevar_max <- var(results[, 2])
  dat_boxplot <- data.frame(matrix(nrow = length(results$thetahat_min), ncol = 4))
  colnames(dat_boxplot) <- c("L", "bootvar", "truevar", "type")
  dat_boxplot$L <- L
  dat_boxplot$bootvar <- results$bootvar_min
  dat_boxplot$truevar <- truevar_min
  dat_boxplot$type <- "minimum"
  list_boxplot[[k]] <- dat_boxplot
  k <- k + 1
  dat_boxplot <- data.frame(matrix(nrow = length(results$thetahat_min), ncol = 4))
  colnames(dat_boxplot) <- c("L", "bootvar", "truevar", "type")
  dat_boxplot$L <- L
  dat_boxplot$bootvar <- results$bootvar_max
  dat_boxplot$truevar <- truevar_max
  dat_boxplot$type <- "maximum"
  list_boxplot[[k]] <- dat_boxplot
  k <- k + 1
}
dat_boxplot <- bind_rows(list_boxplot)
facet_labels <- c(TeX(r'(Var($\hat{\theta}_{min}$|EC Data))'),
                  TeX(r'(Var($\hat{\theta}_{max}$|EC Data))'))
facet_labels <- c("", "")
names(facet_labels) <- c("minimum", "maximum")
p1 <- ggplot(dat_boxplot, aes(x = factor(L), y = sqrt(bootvar), fill = type)) +
  geom_boxplot(alpha = 0.4, outlier.shape = NA, color = "#888888") +
  geom_line(aes(y = sqrt(truevar), group = type, color = type), lwd=0.5, lty = "longdash") +
  geom_point(aes(y = sqrt(truevar), group = type, shape = type, color = type), size = 2) +
  scale_shape_discrete("Empirical s.d.", labels = c(TeX(r'($\hat{\theta}_{max}$)'), TeX(r'(\hat${\theta}_{min}$)'))) +
  scale_color_discrete("Empirical s.d.", labels = c(TeX(r'($\hat{\theta}_{max}$)'), TeX(r'($\hat{\theta}_{min}$)'))) +
  scale_fill_discrete("Estimated s.d.", labels = c(TeX(r'($\hat{\theta}_{max}$)'), TeX(r'($\hat{\theta}_{min}$)'))) +
  ylim(c(0, 0.3)) +
  labs(x = "L (Lipschitz constant)", y = TeX(r'(sd($\hat{\theta}$|EC Summary))')) +
  facet_grid(type~., labeller = labeller(type = facet_labels)) +
  theme_minimal() +
  theme(legend.position = "bottom")

ggsave("plots/experiment3/boxplot.png", plot = p1, width = 14, height = 10, dpi = 400, units = "cm")
