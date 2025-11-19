#' Apply MAIC on a tabulated data coming from tabulate(trt)
#' Return the estimates \eta_1 and \eta_2 used to define the parametric density ratio
MAIC_categorized <- function(categorized, control_moments) {
  FUN <- function(par) {
    eta1 <- par[1]
    eta2 <- par[2]
    f <- function(x, prop, eta1, eta2, control_moments) {
      exp(eta1 * (x - control_moments[1]) + eta2 * (x^2 - control_moments[2])) *
        prop
    }
    sum(f(categorized$X, categorized$prop,
          eta1 = eta1, eta2 = eta2,
          control_moments = control_moments))
  }
  eta_hat <- optim(par = c(0, 0), fn = FUN, method = "BFGS")$par
  eta_hat
}

