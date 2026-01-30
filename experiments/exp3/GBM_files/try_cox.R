rm(list = ls())
# Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(survival)

dat <- read.csv("GBM_complete.csv", header = TRUE, stringsAsFactors = TRUE)
dat <- dat %>% select(age, sex, kps, eor, mgmt, os_status, os)

# Fit a Cox PH model
cox_model <- coxph(Surv(os, os_status) ~ age + sex + kps + eor + mgmt,
                   data = dat)

summary(cox_model)

# Construct a dataset for prediction
# Suppose we want to predict the survival probability for the first individual
newdat <- dat[rep(2, 110), ]
colnames(newdat) <- colnames(dat)

# Change "os" column to time at which we want to compute the survival probability
newdat$os <- seq(from = 0, to = 40, length.out = nrow(newdat))
newdat$os_status <- 1 # Placeholder. Not used.

# Now predict based on the fitted model
plot(newdat$os, predict(cox_model, newdata = newdat, type = "survival"),
     ty = "l", xlab = "t", ylab = "S(t|X) = P(T > t|X)", ylim = c(0, 1))

# A function to get the survival probabilities