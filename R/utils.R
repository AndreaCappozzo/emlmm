#' @export
control_EM <- function(itermax = 1000,
                       tol = 1e-08,
                       err = .Machine$double.xmax / 2,
                         BETA_zero_tol=1e-7, CD_threshold=1e-8) {
  list(
    itermax = itermax,
    tol = tol,
    err = err,
    BETA_zero_tol = BETA_zero_tol,
    CD_threshold = CD_threshold # Convergence threshold for coordinate descent
  )
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
    "net-reg" = update_BETA_netreg,
    "sparse-grouplasso"=update_BETA_sparse_group_lasso
  )
}

penalty_value_f <-
  function(penalty_type,
           BETA,
           alpha,
           lambda,
           lambda_X = 0, # lambda_1 and lambda_2 are used in net-reg only
           lambda_Y = 0,
           G_X = NULL,
           G_Y = NULL) {
    switch(
      penalty_type,
      "elastic-net" = lambda*((1 - alpha) * norm(BETA[-1, , drop = FALSE], type = "F") ^ 2 / 2 + # convex comb of L2 and L1 norm
                               alpha*(sum(abs(BETA[-1,])))),
      "group-lasso" = lambda * ((1 - alpha) * norm(BETA[-1, , drop = FALSE], type = "F") ^
                                  2 / 2 + # [-1,] cos the intercepts are not penalized
                                  alpha * sum(apply(BETA[-1, , drop = FALSE], 1, function(b)
                                    sqrt(sum(b ^ 2))))),
      "sparse-grouplasso" = lambda * (alpha * (sum(abs(BETA[-1,]))) + # [-1,] cos the intercepts are not penalized
                                    (1 - alpha) * sum(apply(BETA[-1, , drop = FALSE], 1, function(b)
                                    sqrt(sum(b ^ 2))))),
      "net-reg" = lambda*(sum(abs(BETA[-1,]))) +
        ifelse(lambda_X==0,0,lambda_X*sum(diag((t(BETA[-1,])%*%(diag(colSums(G_X))-G_X)%*%BETA[-1,]))))+
        ifelse(lambda_Y==0,0,lambda_Y*sum(diag((t(BETA[-1,])%*%(diag(colSums(G_Y))-G_Y)%*%BETA[-1,]))))
    )
  }


# Different updates for BETA ----------------------------------------------

update_BETA_group_lasso <- function(X,Y,alpha,lambda,CD_threshold,...){
  penalized_regression <-
    glmnet::glmnet(
      x = X,
      y = Y ,
      family = "mgaussian",
      thresh = CD_threshold,
      standardize = FALSE,
      alpha = alpha,
      lambda = lambda
    )
  as.matrix(Reduce(f = cbind,stats::coef(penalized_regression)))
}

update_BETA_sparse_group_lasso <-
  function(X, Y, alpha, lambda, CD_threshold, ...) {
    penalized_regression <-
      lsgl::fit(
        x = X,
        y = Y,
        intercept = TRUE,
        alpha = alpha,
        lambda = lambda,
        d = 1,
        algorithm.config = lsgl::lsgl.algorithm.config(tolerance_penalized_main_equation_loop = CD_threshold,
                                                       verbose = FALSE)
      )

    as.matrix(penalized_regression$beta[[1]])
  }

update_BETA_netreg <-
  function(X,
           Y,
           lambda,
           lambda_X=0,
           lambda_Y=0,
           G_X=NULL,
           G_Y=NULL,
           CD_threshold,
           ...) {

  penalized_regression <-
    netReg::edgenet(
      X = X,
      Y = Y,
      G.X=G_X,
      G.Y=G_Y,
      family = "gaussian",
      lambda = lambda,
      psigx = lambda_X,
      psigy = lambda_Y,
      thresh = CD_threshold
    )

  rbind(penalized_regression$alpha,
        penalized_regression$beta)
}

update_BETA_elastic_net <- function(X,Y,alpha,lambda, I_r,CD_threshold,...){
  p_no_intercept <- ncol(X)
  r <- ncol(Y)
  vec_Y <- c(Y)
  vec_X <- I_r%x%cbind(1,X)
  lasso_weigths <- rep(1,ncol(vec_X))
  lasso_weigths[seq(1,ncol(vec_X),by=(p_no_intercept+1))] <- 0 # in this way I do not penalize the r intercepts
  penalized_regression <-
    glmnet::glmnet(
      x = vec_X,
      y = vec_Y ,
      family = "gaussian",
      standardize = FALSE,
      alpha = alpha,
      thresh = CD_threshold,
      penalty.factor = lasso_weigths,
      intercept = FALSE, # I do not need a single intercept as I need to estimate r intercepts
      lambda = lambda
    )

  matrix(stats::coef(penalized_regression)[-1],nrow = p_no_intercept+1,ncol = r) # I remove the intercept with [-1] (that in any case was not estimated since intercept = FALSE)
}
