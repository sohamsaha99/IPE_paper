prepare_summary <- function(df_trt, df_con, tau) {
  # Create Summary statistics
  target <- data.frame(xname = c("age", "age", "sex", "kps", "mgmt", "isGTR"),
                       type  = c(1, 2, 1, 1, 1, 1)
  )
  target$value <- c(
                    mean(df_con$age),
                    mean(df_con$age^2),
                    mean(df_con$sex),
                    mean(df_con$kps),
                    mean(df_con$mgmt),
                    mean(df_con$isGTR)
  )

  control_rmst_value <- weighted_rmst(df_con$os, df_con$os_status, tau)$rmst[1]
  list(df_trt = df_trt, target_con = target, control_rmst_value = control_rmst_value)
}
