library(MatchIt)
library(survival)
library(survRM2)

#####
# ps.mod <- glm(
#   treat ~ age + sex + kps + mgmt + isGTR,
#   data = df,
#   family = binomial
# )
#
# df$ps <- predict(ps.mod, type = "response")
# df$w <- with(df,
#   treat / ps + (1 - treat) / (1 - ps)
# )
# rmst.w <- rmst2(
#   time = df$os,
#   status = df$os_status,
#   arm = df$treat,
#   tau = tau,
#   weight = df$w
# )
#
# rmst.w$unadjusted.result
#

# Function for finding RMST difference after matching
estimate_trt_effet_matching <- function(df_trt, df_con, tau) {
  df_trt$treat <- 1
  df_con$treat <- 0

  df <- rbind(df_trt, df_con)

  m.out <- matchit(
                   treat ~ age + sex + kps + mgmt + isGTR,
                   data = df,
                   method = "nearest",
                   distance = "logit",
                   ratio = 1
  )

  summary(m.out)
  df.matched <- match.data(m.out)

  rmst.out <- rmst2(
                    time = df.matched$os,
                    status = df.matched$os_status,
                    arm = df.matched$treat,
                    tau = tau
  )

  theta_hat <- rmst.out$unadjusted.result[1, 1]
  lower <- rmst.out$unadjusted.result[1, 2]
  upper <- rmst.out$unadjusted.result[1, 3]
  pvalue <- rmst.out$unadjusted.result[1, 4]

  list(
       theta_hat = theta_hat,
       lower = lower,
       upper = upper,
       pvalue = pvalue
  )
}
# estimate_trt_effet_matching(df_trt, df_con, 24)
