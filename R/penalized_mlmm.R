# Y is a matrix of responses, with r the number of response variables
#' @export
ecm_mlmm_penalized <-
  function(X,
           Y,
           Z,
           group,
           lambda,
           alpha=1, # alpha is for elastic net type of penalty for the group_lasso penalty_type only. alpha=1 means lasso
           penalty_type = c("group_lasso", "netReg"),
           control_EM_algorithm = control_EM()) {


    X <-  data.matrix(X)
    X_no_intercept <-
      X[, -1, drop = FALSE] # needed for penalized regression
    Y <-  data.matrix(Y)
    Z <-  data.matrix(Z)

    q <- ncol(Z) # # ran eff
    p <- ncol(X) # # fix eff
    r <- ncol(Y) # # of responses
    N <- nrow(X)
    J <- length(unique(group))

    # Set initial values

    SIGMA <- diag(r)
    PSI <- diag(q*r)

    if (lambda == 0) {
    BETA <- (stats::lm.fit(y = Y, x = X))$coefficients

    } else{
    #FIXME define a function that returns BETA for which the way BETA is computed depends on penalty_type
      penalized_regression_0 <-
        glmnet::glmnet(
          x = X_no_intercept,
          y = Y , # I do not premultiply by SIGMA as in the init SIGMA is a diagonal matrix
          family = "mgaussian",
          alpha = alpha,
          lambda = lambda/N
        )
      BETA <-as.matrix(Reduce(f = cbind,stats::coef(penalized_regression_0)))
    }



    # Objects to which I apply the vec operator and other useful quantities
    vec_XB <- c(X %*% BETA) # Nr x 1
    vec_Y <- c(Y) # Nr x 1
    I_r <- diag(r)

    # EM parameters
    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    err <- control_EM_algorithm$err

    iter <- 0
    loglik_pen <- loglik_pen_prev <- -.Machine$integer.max / 2
    loglik_pen_vec <- NULL

    crit <- TRUE
    group_indicator <- as.numeric(group)
    vec_group_indicator <- rep(group_indicator,r)

    while (crit) {

      # Compute useful quantities
      vec_res_fixed <- matrix((vec_Y - vec_XB),ncol = 1)
      PSI_inv <- solve(PSI)

      eig_SIGMA <- eigen(SIGMA, symmetric = TRUE)
      SIGMA_half <- eig_SIGMA$vectors%*%diag(sqrt(eig_SIGMA$values))%*%t(eig_SIGMA$vectors)
      SIGMA_inv_half <- eig_SIGMA$vectors%*%diag(sqrt(1/eig_SIGMA$values))%*%t(eig_SIGMA$vectors)
      SIGMA_inv <- SIGMA_inv_half%*%SIGMA_inv_half

      # E step ------------------------------------------------------------------

      e_step_lmm <- estep_mlmm_cpp(
        vec_res_fixed = vec_res_fixed,
        Z = Z,
        group_indicator = group_indicator,
        vec_group_indicator = vec_group_indicator,
        inv_Psi = PSI_inv,
        inv_Sigma = SIGMA_inv,
        I_r=I_r,
        r = r,
        J = J
      )

      vec_raneff_i <- e_step_lmm$vec_raneff_i
      est_second_moment <- e_step_lmm$est_second_moment

      raneff_i <- matrix(vec_raneff_i, nrow = N, ncol = r)

      # M step ------------------------------------------------------------------

      # FIXME wrap BETA calculation into a function that depends on the penalty

      Sigma_inv_half_Y_tilde <- t(tcrossprod(x = SIGMA_inv_half, y = Y - raneff_i))

      penalized_regression <-
        glmnet::glmnet(
          x = X_no_intercept,
          y = Sigma_inv_half_Y_tilde,
          family = "mgaussian",
          standardize = FALSE,
          alpha = alpha,
          lambda = (lambda)/N
        )

      BETA <-
        as.matrix(Reduce(f = cbind,stats::coef(penalized_regression))%*%SIGMA_half)

      PSI <- as.matrix(est_second_moment / J)

      SIGMA <- stats::cov(Y - X %*% BETA - raneff_i)*(N-1)/N
      vec_XB <- c(X %*% BETA) # Nr x 1

      #### log lik evaluation-------------------------------------------------

      loglik <- log_lik_mlmm_cpp(
        vec_Y=vec_Y,
        vec_XB=vec_XB,
        Z = Z,
        group_indicator = group_indicator,
        vec_group_indicator=vec_group_indicator,
        PSI = PSI,
        SIGMA=SIGMA,
        I_r=I_r,
        J = J
      )

      penalty_value <-
        (1-alpha)*norm(BETA[-1, , drop = FALSE],type = "F")^2/2+ # [-1,] cos the intercepts are not penalized
        alpha*sum(apply(BETA[-1, , drop = FALSE], 1, function(b)
          sqrt(sum(b ^ 2))))#FIXME wrap it in a function that changes according to the penalty type

      loglik_pen <-
        loglik - lambda * penalty_value

      # check convergence
      err <-
        abs(loglik_pen - loglik_pen_prev) / (1 + abs(loglik_pen))
      loglik_pen_prev <- loglik_pen
      loglik_pen_vec <- c(loglik_pen_vec, loglik_pen)
      iter <- iter + 1
      crit <- (err > tol & iter < itermax)
    }

    mu_raneff <- array(c(e_step_lmm$mat_mu_raneff),dim=c(q,r,J)) # FIXME check when q>1

    OUT <-       list(
      BETA = BETA,
      PSI = PSI,
      SIGMA = SIGMA,
      mu_raneff = mu_raneff,
      loglik = loglik,
      loglik_pen = loglik_pen,
      loglik_pen_trace = loglik_pen_vec,
      lambda = lambda,
      penalty_value = penalty_value,
      iter = iter
    )

    class(OUT) <- "mlmm"
    OUT
  }
