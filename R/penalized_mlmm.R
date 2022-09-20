# Y is a matrix of responses, with r the number of response variables
#' @export
ecm_mlmm_penalized <-
  function(X,
           Y,
           Z,
           group,
           lambda,
           lambda_X=0, # lambda_X, lambda_Y, G_X and G_Y are required only when penalty_type== "net-reg"
           lambda_Y=0,
           G_X=NULL,
           G_Y=NULL,
           alpha=1, # alpha is for elastic net type of penalty aka group-lasso and elastic-net penalty_type only. alpha=1 means lasso
           penalty_type = c("group-lasso", "elastic-net","net-reg"),
           control_EM_algorithm = control_EM()) {
    start_time <- Sys.time()
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
    CD_threshold <- control_EM_algorithm$CD_threshold

    if (lambda == 0) {
    BETA <- (stats::lm.fit(y = Y, x = X))$coefficients

    } else{
      BETA <-
        update_BETA(
          X = X_no_intercept,
          Y = Y,
          # I do not premultiply by SIGMA as in the init SIGMA is a diagonal matrix
          alpha = alpha,
          lambda = lambda / N,
          CD_threshold = CD_threshold,
          I_r = I_r,
          lambda_X = lambda_X,
          lambda_Y = lambda_Y,
          G_X = G_X,
          G_Y = G_Y
        )
    }


    vec_XB <- c(X %*% BETA) # Nr x 1


    # EM parameters
    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    err <- control_EM_algorithm$err
    BETA_zero_tol <- control_EM_algorithm$BETA_zero_tol


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
      #
      vec_raneff_i <- e_step_lmm$vec_raneff_i
      est_second_moment <- e_step_lmm$est_second_moment
      est_second_moment_error <- e_step_lmm$est_second_moment_error

      raneff_i <- matrix(vec_raneff_i, nrow = N, ncol = r)

      # # Original R code # FIXME implement everything in cpp
      # est_second_moment_error <- matrix(0,nrow = r, ncol = r)
      # est_second_moment <- matrix(0,nrow = q*r, ncol = q*r)
      # mu_raneff <- array(dim=c(q,r,J))
      # vec_raneff_i <- vector(mode = "numeric", length = length(vec_Y))
      #
      #
      # for (j in 1:J) {
      #   # iterate over different groups
      #   rows_j <- which(group_indicator == j)
      #   vec_rows_j <- which(vec_group_indicator == j)
      #   n_j <- length(rows_j)
      #   I_nj <- diag(n_j)
      #
      #   Z_j <- Z[rows_j, , drop = FALSE]
      #
      #   vec_res_fixed_j <- vec_res_fixed[vec_rows_j, , drop = FALSE]
      #
      #   common_component_j <-
      #     (I_r %x% t(Z_j)) %*% (SIGMA_inv %x% I_nj)
      #   Gamma_j <-
      #     solve(common_component_j %*% (I_r %x% Z_j) + PSI_inv)
      #   vec_DELTA_j <-
      #     Gamma_j %*% common_component_j %*% vec_res_fixed_j
      #   mu_raneff[, , j] <- matrix(vec_DELTA_j, nrow = q, ncol = r)
      #   vec_raneff_i[vec_rows_j] <- (I_r %x% Z_j) %*% vec_DELTA_j
      #   est_second_moment <-
      #     est_second_moment + Gamma_j + vec_DELTA_j %*% t(vec_DELTA_j)
      #
      #   var_ei <- (I_r%x%Z_j) %*% Gamma_j %*% (I_r%x%t(Z_j))
      #   slice_rows <- slice_cols <- split(1:(n_j*r), ceiling(seq_along(1:(n_j*r))/n_j))
      #
      #   for (row in 1:r) {
      #     for (col in 1:r) {
      #       est_second_moment_error[row, col] <-
      #         est_second_moment_error[row, col] + sum(diag(var_ei[slice_rows[[row]], slice_cols[[col]]]))
      #     }
      #   }
      # }
      #
      # raneff_i <- matrix(vec_raneff_i, nrow = N, ncol = r)

      # M step ------------------------------------------------------------------

      Sigma_inv_half_Y_tilde <- t(tcrossprod(x = SIGMA_inv_half, y = Y - raneff_i))

      BETA_tmp <-
        update_BETA(
          X = X_no_intercept,
          Y = Sigma_inv_half_Y_tilde,
          alpha = alpha,
          lambda = lambda / N,
          I_r = I_r,
          CD_threshold=CD_threshold
        )

      BETA <-
        (tcrossprod(BETA_tmp, SIGMA_half))

      BETA <- # should I have values of BETA very close to 0, I can set them to be exactly 0 if they are smaller than BETA_zero_tol
        ifelse(abs(BETA) >= BETA_zero_tol |
                 lambda == 0 | penalty_type == "group-lasso",
               BETA,
               0)

      # BETA_fixed <-
      #   update_BETA(
      #     X = X_no_intercept,
      #     Y = (Y - raneff_i),
      #     alpha = alpha,
      #     lambda = lambda / N,
      #     I_r = I_r,
      #     CD_threshold=CD_threshold
      #   )
      #
      # BETA <- tcrossprod(BETA_tmp, SIGMA_half)*(abs(BETA_fixed)>0)

      # BETA <- # Correct but numerically unstable for elastic-net and net-reg
      #   tcrossprod(update_BETA(
      #     X = X_no_intercept,
      #     Y = Sigma_inv_half_Y_tilde,
      #     alpha = alpha,
      #     lambda = lambda / N,
      #     I_r = I_r,
      #     CD_threshold=CD_threshold
      #   ), SIGMA_half) # I postmultiply by SIGMA_half as per derivation in the paper

      PSI <- as.matrix(est_second_moment / J)

      # SIGMA <- stats::cov(Y - X %*% BETA - raneff_i)*(N-1)/N
      SIGMA <- (crossprod(Y - X %*% BETA - raneff_i) +est_second_moment_error)/N
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
          lambda = lambda,
          lambda_X = lambda_X,
          lambda_Y = lambda_Y,
          G_X = G_X,
          G_Y = G_Y
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

    # BETA <- # should I have values of BETA very close to 0, I can set them to be exactly 0 if they are smaller than BETA_zero_tol
    #   ifelse(abs(BETA) >= BETA_zero_tol |
    #            lambda == 0 | penalty_type == "group-lasso",
    #          BETA,
    #          0)

    mu_raneff <- array(c(e_step_lmm$mat_mu_raneff),dim=c(q,r,J)) # FIXME check when q>1
    end_time <- Sys.time()
    elapsed_time <- difftime(end_time,start_time,units = "secs")
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
      iter = iter,
      elapsed_time=elapsed_time
    )

    class(OUT) <- "mlmm"
    OUT
  }
