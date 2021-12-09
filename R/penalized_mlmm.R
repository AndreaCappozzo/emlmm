# Y is a matrix of responses, with r the number of response variables
#' @export
ecm_mlmm_penalized <-
  function(X,
           Y,
           Z,
           group,
           lambda,
           alpha=1, # alpha is for elastic net type of penalty aka group-lasso and elastic-net penalty_type only. alpha=1 means lasso
           penalty_type = c("group-lasso", "elastic-net","net-reg"),
           control_EM_algorithm = control_EM()) {

    # Depending on the chosen penalty, I define the corresponding function to be used in the update in the M-step
    update_BETA <- update_BETA_f(penalty_type = penalty_type)

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

    # Objects to which I apply the vec operator and other useful quantities
    vec_Y <- c(Y) # Nr x 1
    I_r <- diag(r)

    if (lambda == 0) {
    BETA <- (stats::lm.fit(y = Y, x = X))$coefficients

    } else{
      BETA <-
        update_BETA(
          X = X_no_intercept,
          Y = Y, # I do not premultiply by SIGMA as in the init SIGMA is a diagonal matrix
          alpha = alpha,
          lambda = lambda / N,
          I_r = I_r
        )
    }


    vec_XB <- c(X %*% BETA) # Nr x 1


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

      Sigma_inv_half_Y_tilde <- t(tcrossprod(x = SIGMA_inv_half, y = Y - raneff_i))

      BETA <-
        update_BETA(
          X = X_no_intercept,
          Y = Sigma_inv_half_Y_tilde,
          alpha = alpha,
          lambda = lambda / N,
          I_r = I_r
        )%*%SIGMA_half # I postmultiply by SIGMA_half as per derivation in the paper

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
        penalty_value_f(
          penalty_type = penalty_type,
          BETA = BETA,
          alpha = alpha,
          lambda = lambda
        )

      loglik_pen <-
        loglik - penalty_value

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
