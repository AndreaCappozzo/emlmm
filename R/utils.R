#' @export
control_EM <- function(itermax = 1000,
                       tol = 1e-08,
                       err = .Machine$double.xmax / 2) {
  list(itermax = itermax,
       tol = tol,
       err = err)
}

# Function that checks which columns of X are also in the random effect Z, returning a penalty factor that identifies which column of
# X to penalize. Note: if the intercept is found in Z then that element is eliminated from the final vector, so it can be used in
# glmnet where the intercept is never penalized
#' @export
pen_factor_builder <- function(X,Z){
  # Note: X and Z may contain the intercept
  p <- ncol(X)
  q <- ncol(Z)
  penalty_factor_glmnet <- rep(1, ncol(X))

  for (q_col in 1:q) {
    for (p_col in 1:p) {
      if(all(Z[,q_col]==X[,p_col])){
        penalty_factor_glmnet[p_col] <- 0
        break
      }
    }
  }
  # Note: we assume that the intercept is always the first column of both X and Z
  if(all(Z[,1]==1) & all(X[,1]==1)){
    penalty_factor_glmnet <- penalty_factor_glmnet[-1]
  }
  penalty_factor_glmnet
}
