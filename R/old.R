# Univariate lmm lasso with the same transformation as per the multivariate version: I divide y by lambda and I
# postmultiply the resulting betas at each M iteration. Results are equivalent to ecm_lmm_lasso
# ecm_lmm_lasso_2 <-
#   function(X ,
#            y,
#            Z,
#            group,
#            lambda = NULL,
#            control_EM_algorithm = control_EM()) {
#     X <-  data.matrix(X)
#     X_no_intercept <-
#       X[, -1, drop = FALSE] # needed for penalized regression
#     # y <-  data.matrix(y)
#     Z <-  data.matrix(Z)
#     q <- ncol(Z) # # ran eff
#     p <- ncol(X) # # fix eff
#     N <- nrow(X)
#     J <- length(unique(group)) # number of groups
#
#     if (is.null(lambda)) { #FIXME
#       # if lambda is not provided, an initial guess is obtained by performing 10-CV on the fixed effect model only
#       fit_cv_glmnet <-
#         glmnet::cv.glmnet(
#           x = X_no_intercept,
#           y = y,
#           nfolds = 10,
#           alpha = 1,
#           family = "gaussian",
#           standardize = FALSE
#         )
#       lambda_range <-
#         range(fit_cv_glmnet$lambda)*N # this should be used to have a guess about the range in which lambda may vary for the given dataset
#       lambda <- fit_cv_glmnet$lambda.min*N
#     } else{
#       lambda_range = NULL
#     }
#
#     # Set initial values
#     sigma2 <- 1
#     Omega <- diag(q)
#
#     if (lambda == 0) {
#       beta <- (stats::lm.fit(y = y, x = X))$coefficients
#     } else{
#       penalized_regression_0 <-
#         glmnet::glmnet(
#           x = X_no_intercept,
#           y = y/sigma2 ,
#           family = "gaussian",
#           # standardize = FALSE,
#           alpha = 1,
#           lambda = (lambda)/N  # NOTE: the obj func in Rohart 2014 is multiplied by 2 wrt to the one I am using, and I think they should divide lambda by n as well
#         )
#       beta <-
#         as.vector(stats::coef(penalized_regression_0)*sigma2)
#     }
#
#     # EM parameters
#     itermax <- control_EM_algorithm$itermax
#     tol <- control_EM_algorithm$tol
#     err <- control_EM_algorithm$err
#     iter <- 0
#     # I need to monitor the penalized likelihood
#     loglik_pen <- loglik_pen_prev <- -.Machine$integer.max / 2
#     loglik_pen_vec <- NULL
#
#     crit <- TRUE
#     group_indicator <- as.numeric(group)
#
#     while (crit) {
#       res_fixed <- y - X %*% beta
#
#
#       # E step ------------------------------------------------------------------
#
#       e_step_lmm <- estep_lmm_cpp(
#         res_fixed = res_fixed,
#         Z = Z,
#         group_indicator = group_indicator,
#         sigma2 = sigma2,inv_Omega = solve(Omega),
#         J = J
#       )
#
#       raneff_i <- e_step_lmm$raneff_i
#       est_second_moment <- e_step_lmm$est_second_moment
#       est_second_moment_error <- e_step_lmm$est_second_moment_error
#
#       # M step ------------------------------------------------------------------
#
#       penalized_regression <-
#         glmnet::glmnet(
#           x = X_no_intercept,
#           y = (y - raneff_i)/sigma2,
#           family = "gaussian",
#           standardize = FALSE,
#           alpha = 1,
#           lambda = (lambda)/N
#         )
#
#       beta <-
#         as.vector(stats::coef(penalized_regression))*sigma2 # I include the intercept in the set of estimated parameters
#       Omega <- as.matrix(est_second_moment / J)
#       sigma2 <- c(stats::var(y - X %*% beta - raneff_i)*(N-1)/N) + est_second_moment_error/N
#
#       #### log lik evaluation-------------------------------------------------
#
#       loglik <- log_lik_lmm_cpp(
#         y = y,
#         Z = Z,
#         X = X,
#         group_indicator = group_indicator,
#         beta = beta,
#         Omega = Omega,
#         sigma2 = sigma2,
#         J = J
#       )
#
#       loglik_pen <-
#         loglik - lambda * sum(abs(beta[-1])) # objective function FIXME once penalty_factor glm is computed
#
#       # check convergence
#       err <-
#         abs(loglik_pen - loglik_pen_prev) / (1 + abs(loglik_pen))
#       loglik_pen_prev <- loglik_pen
#       loglik_pen_vec <- c(loglik_pen_vec, loglik_pen)
#       iter <- iter + 1
#       crit <- (err > tol & iter < itermax)
#     }
#
#     # Add parameters names
#     names(beta) <- rownames(stats::coef(penalized_regression))
#
#     # Collect results
#     return(
#       list(
#         beta = beta,
#         Omega = Omega,
#         sigma2 = sigma2,
#         mu_raneff = e_step_lmm$mu_raneff,
#         loglik = loglik,
#         loglik_pen = loglik_pen,
#         loglik_pen_trace = loglik_pen_vec,
#         lambda = lambda,
#         lambda_range = lambda_range,
#         iter=iter
#       )
#     )
#   }
