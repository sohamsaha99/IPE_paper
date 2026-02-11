# Function for generating event time based on Weibull regression fitted model
ET<-function(dat, trt_eff, a_shape, b_scale, beta_coefficient, interaction_eff = 0.0){
  Z_a <- model.matrix( ~ age + sex + kps + mgmt + eor, 
                       data = dat)
  Z_a <- Z_a[, -1]
  Z_a <- Z_a - matrix(colMeans(Z_a), nrow = nrow(Z_a), ncol = ncol(Z_a), 
                      byrow = T)
  # print(colMeans(Z_a))
  # if (ncol(Z_a) != 6) {
  #   print(Z_a)
  # }
  lin_pred_a <- as.vector(Z_a %*% beta_coefficient[colnames(Z_a)])
  # Modify the linear predictor to include interaction between kps and eor
  lin_pred_a <- lin_pred_a + interaction_eff * (dat$kps == 1) * (dat$mgmt == 1)
  print(sum((dat$kps == 1) * (dat$mgmt == 1)))
  # Apply treatment effect depending on treatment assignment
  trt_eff <- trt_eff * dat$isTrial
  # print(length(lin_pred_a)); print(length(trt_eff))
  # print(lin_pred_a)
  stopifnot(length(lin_pred_a) == length(trt_eff))
  x_event_time <- rweibull(nrow(dat),
                           shape = a_shape,
                           scale = b_scale*(exp(-(lin_pred_a - trt_eff) / a_shape)))
  # print(x_event_time)
  return(x_event_time)
}
