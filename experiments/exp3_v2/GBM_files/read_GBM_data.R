
df_GBM <- read.csv("gbm_TRG.csv", header = TRUE, stringsAsFactors = FALSE)
# Read GBM data file
df_GBM[, "eor"] <- as.factor(df_GBM[, "eor"])
df_GBM <- df_GBM[, -1]
# Read output of Weibull regression
weibull_model <- read.csv("weib_2.csv", header = TRUE, stringsAsFactors = FALSE)
# Pre-processing of Weibull regression parameters
weib_reg <- list()
weib_reg$weibull_beta <- weibull_model$wr.coefficients[1:6]
names(weib_reg$weibull_beta) <- c("age", "sex", "kps", "mgmt", "eorGTR", "eorSTR")
weib_reg$weibull_scale <- exp(weibull_model$wr.coefficients[7])
weib_reg$weibull_shape <- exp(weibull_model$wr.coefficients[8])
