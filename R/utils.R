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

update_BETA_f <- function(penalty_type) {
  switch(
    penalty_type,
    "elastic-net" = update_BETA_elastic_net,
    "group-lasso" = update_BETA_group_lasso,
    "net-reg" = update_BETA_netreg
  )
}

penalty_value_f <-
  function(penalty_type,
           BETA,
           alpha,
           lambda,
           lambda_1=0, # lambda_1 and lambda_2 are used in net-reg only
           lambda_2=0) {
    switch(
      penalty_type,
      "elastic-net" = lambda*((1 - alpha) * norm(BETA[-1, , drop = FALSE], type = "F") ^ 2 / 2 + # convex comb of L2 and L1 norm
                               alpha*(sum(abs(BETA[-1,])))),
      "group-lasso" = lambda * ((1 - alpha) * norm(BETA[-1, , drop = FALSE], type = "F") ^
                                  2 / 2 + # [-1,] cos the intercepts are not penalized
                                  alpha * sum(apply(BETA[-1, , drop = FALSE], 1, function(b)
                                    sqrt(sum(b ^ 2))))),
      "net-reg" = "FIXME"
    )
  }


# Different updates for BETA ----------------------------------------------

update_BETA_group_lasso <- function(X,Y,alpha,lambda,...){
  penalized_regression <-
    glmnet::glmnet(
      x = X,
      y = Y ,
      family = "mgaussian",
      standardize = FALSE,
      alpha = alpha,
      lambda = lambda
    )
  as.matrix(Reduce(f = cbind,stats::coef(penalized_regression)))
}

update_BETA_netreg <- function(X,Y,alpha,lambda, lambda_1,lambda_2,...){
  penalized_regression <-
    glmnet::glmnet(
      x = X,
      y = Y ,
      family = "mgaussian",
      standardize = FALSE,
      alpha = alpha,
      lambda = lambda
    )
  as.matrix(Reduce(f = cbind,stats::coef(penalized_regression)))
}

update_BETA_elastic_net <- function(X,Y,alpha,lambda, I_r,...){
  p_no_intercept <- ncol(X)
  r <- ncol(Y)
  vec_Y <- c(Y)
  vec_X <- I_r%x%cbind(1,X)
  penalized_regression <-
    glmnet::glmnet(
      x = vec_X,
      y = vec_Y ,
      family = "gaussian",
      standardize = FALSE,
      alpha = alpha,
      # thresh = .Machine$double.eps,
      intercept = FALSE, # I do not need a single intercept as I need to estimate r intercepts
      lambda = lambda
    )

  matrix(stats::coef(penalized_regression)[-1],nrow = p_no_intercept+1,ncol = r) # I remove the intercept with [-1] (that in any case was not estimated since intercept = FALSE)
}
