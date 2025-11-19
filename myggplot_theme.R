# my_theme.R

library(ggplot2)

theme_mydefault <- function(base_size = 16, base_family = "sans") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
          plot.title   = element_text(size = 1.4 * base_size, face = "bold", hjust = 0.5),
          axis.title   = element_text(size = 1.2 * base_size),
          axis.text    = element_text(size = 1.0 * base_size),
          legend.title = element_text(size = 1.1 * base_size),
          legend.text  = element_text(size = 0.9 * base_size)
    )
}

