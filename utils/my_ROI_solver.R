library(ROI)
library(ROI.plugin.highs)
library(ROI.plugin.lpsolve)
library(ROI.plugin.glpk)

#' Function to solve an optimization problem defined using OP
#' We shall try the available solvers in sequence until we get success
#' We shall return an object similar to what ROI_solve returns
my_ROI_solver <- function(op, ...) {
  avaialable_solvers <- c("glpk", "lpsolve", "highs")
  for (method in avaialable_solvers) {
    model <- try(ROI_solve(op, solver = method, ...), silent = TRUE)
    if (!inherits(model, "try-error") && model$status$code == 0) {
      model$solver <- method
      return(model)
    }
  }
  if (inherits(model, "try-error")) {
    print("[[my_ROI_solver]]- All solvers encountered runtime error or failed.")
    return(list(solver = "NONE", status = list(code = 1)))
  }
  print("[[my_ROI_solver]]- All solvers failed")
  model$solver <- "NONE"
  model
}

# Make dual of the problem min cx subject to Ax dir b, x>= 0
make_dual_lp <- function(A, b, c, dir) {
  stopifnot(is.matrix(A), is.numeric(A))
  m <- nrow(A); n <- ncol(A)
  stopifnot(length(b) == m, length(c) == n, length(dir) == m)
  if (!all(dir %in% c("<=", ">=", "==")))
    stop("dir must be in {'<=','>=','=='}")

  # Dual variables y have one per primal constraint (m of them)
  # Set bounds for y according to the primal row sense
  lb <- rep(-Inf, m)
  ub <- rep(Inf, m)
  lb[dir == ">="] <- 0 # y >= 0
  ub[dir == "<="] <- 0 # y <= 0
  # dir == "==" stays free (-Inf, +Inf)

  ROI::OP(
          objective   = ROI::L_objective(b),                 # max b^T y
          constraints = ROI::L_constraint(L = t(A),          # A^T y <= c
                                          dir = rep("<=", n),
                                          rhs = c),
          bounds      = ROI::V_bound(lb = lb, ub = ub),      # signs of y
          maximum     = TRUE
  )
}

