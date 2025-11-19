library(dplyr)
#' Given data with X, Y columns, categorize based on X
make_table <- function(dat) {
  # Check if dat has two columns X and Y
  stopifnot(ncol(dat) == 2)
  stopifnot(setequal(colnames(dat), c("X", "Y")))
  # Construct categorized data
  categorized <- dat %>%
    group_by(X) %>%
    summarize(count = n(),
              mean_Y = mean(Y),
              mean_Ysq = mean(Y^2)) %>%
    # tidyr::complete(X = seq(0, 1, by = 0.01),
    #                 fill = list(count = 0,
    #                             mean_Y = 0,
    #                             mean_Ysq = 0)) %>%
    arrange(X)
  categorized$prop <- categorized$count / sum(categorized$count)
  categorized
}
