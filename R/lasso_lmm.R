

# LMM model with lasso penalty on the fixed effects ----------------------

# Standard EM algorithm for LMM with lasso penalty
#' @export
em_lmm_lasso <-
  function(X ,
           y,
           Z,
           group,
           lambda = NULL,
           control_EM_algorithm = control_EM()) {
    X <-  data.matrix(X)
    X_no_intercept <-
      X[, -1, drop = FALSE] # needed for penalized regression
    y <-  data.matrix(y)
    Z <-  data.matrix(Z)
    q <- ncol(Z) # # ran eff
    p <- ncol(X) # # fix eff
    N <- nrow(X)
    J <- length(unique(group)) # number of groups

    if (is.null(lambda)) {
      # if lambda is not provided, an initial guess is obtained by performing 10-CV on the fixed effect model only
      fit_cv_glmnet <-
        glmnet::cv.glmnet(
          x = X_no_intercept,
          y = y,
          nfolds = 10,
          alpha = 1,
          family = "gaussian",
          standardize = FALSE
        )
      lambda_range <-
        range(fit_cv_glmnet$lambda) # this should be used to have a guess about the range in which lambda may vary for the given dataset
      lambda <- fit_cv_glmnet$lambda.min
    } else{
      lambda_range = NULL
    }

    # Set initial values
    if (lambda == 0) {
      beta <- (lm.fit(y = y, x = X))$coefficients
    } else{
      beta <-
        c(mean(y), rep(0, (p - 1))) # I set all betas equal to 0 but the intercept
    }

    sigma2 <- 1
    Omega <- diag(q)

    mu_raneff <- matrix(nrow = q, ncol = J)
    raneff_i <- vector(mode = "numeric", length = N)

    # EM parameters
    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    err <- control_EM_algorithm$err
    iter <- 0
    # I need to monitor the penalized likelihood
    loglik_pen <- loglik_pen_prev <- -.Machine$integer.max / 2
    loglik_pen_vec <- NULL

    crit <- TRUE
    group_indicator <- as.numeric(group)

    while (crit) {
      res_fixed <- y - X %*% beta
      est_second_moment <- 0

      # E step ------------------------------------------------------------------

      for (j in 1:J) {
        # iterate over different groups
        rows_j <- which(group_indicator == j)
        Z_j <- Z[rows_j, , drop = FALSE]
        res_fixed_j <- res_fixed[rows_j, drop = FALSE]
        Gamma_j <- solve(t(Z_j) %*% Z_j / sigma2 + solve(Omega))
        mu_j <- (Gamma_j %*% t(Z_j) %*% res_fixed_j) / sigma2
        mu_raneff[, j] <- mu_j
        raneff_i[rows_j] <- Z_j %*% mu_j
        est_second_moment <-
          est_second_moment + Gamma_j + mu_j %*% t(mu_j)
      }

      # M step ------------------------------------------------------------------

      penalized_regression <-
        glmnet::glmnet(
          x = X_no_intercept,
          y = y - raneff_i,
          family = "gaussian",
          standardize = FALSE,
          alpha = 1,
          lambda = lambda
        )
      beta <-
        as.vector(stats::coef(penalized_regression)) # I include the intercept in the set of estimated parameters
      Omega <- as.matrix(est_second_moment / J)
      sigma2 <- mean(y * (y - X %*% beta - raneff_i))

      #### log lik evaluation-------------------------------------------------

      loglik <- 0

      for (j in 1:J) {
        rows_j <- which(group_indicator == j)
        Z_j <- Z[rows_j, , drop = FALSE]
        y_j <- y[rows_j, drop = FALSE]
        X_j <- X[rows_j, , drop = FALSE]
        G_j <-
          Z_j %*% Omega %*% t(Z_j) + diag(sigma2, nrow = length(rows_j))
        loglik <-
          loglik + mvtnorm::dmvnorm(
            x = y_j,
            mean = c(X_j %*% beta),
            sigma = G_j,
            log = TRUE
          )
      }

      loglik_pen <-
        loglik - lambda * sum(abs(beta)) # objective function

      # check convergence
      err <-
        abs(loglik_pen - loglik_pen_prev) / (1 + abs(loglik_pen))
      loglik_pen_prev <- loglik_pen
      loglik_pen_vec <- c(loglik_pen_vec, loglik_pen)
      iter <- iter + 1
      crit <- (err > tol & iter < itermax)
    }

    return(
      list(
        beta = beta,
        Omega = Omega,
        sigma2 = sigma2,
        mu_raneff = mu_raneff,
        loglik = loglik,
        loglik_pen = loglik_pen,
        loglik_pen_trace = loglik_pen_vec,
        lambda = lambda,
        lambda_range = lambda_range
      )
    )
  }

# ECM algorithm for LMM with lasso penalty as described in Rohart 2014 (http://dx.doi.org/10.1016/j.csda.2014.06.022)
#' @export
ecm_lmm_lasso <-
  function(X ,
           y,
           Z,
           group,
           lambda = NULL,
           control_EM_algorithm = control_EM()) {

    X <-  data.matrix(X)
    X_no_intercept <-
      X[, -1, drop = FALSE] # needed for penalized regression
    y <-  data.matrix(y)
    Z <-  data.matrix(Z)
    q <- ncol(Z) # # ran eff
    p <- ncol(X) # # fix eff
    N <- nrow(X)
    J <- length(unique(group)) # number of groups

    if (is.null(lambda)) {
      # if lambda is not provided, an initial guess is obtained by performing 10-CV on the fixed effect model only
      fit_cv_glmnet <-
        glmnet::cv.glmnet(
          x = X_no_intercept,
          y = y,
          nfolds = 10,
          alpha = 1,
          family = "gaussian",
          standardize = FALSE
        )
      lambda_range <-
        range(fit_cv_glmnet$lambda) # this should be used to have a guess about the range in which lambda may vary for the given dataset
      lambda <- fit_cv_glmnet$lambda.min
    } else{
      lambda_range = NULL
    }

    # Set initial values
    if (lambda == 0) {
      beta <- (lm.fit(y = y, x = X))$coefficients
    } else{
      beta <-
        c(mean(y), rep(0, (p - 1))) # I set all betas equal to 0 but the intercept FIXME
    }

    sigma2 <- 1
    Omega <- diag(q)

    mu_raneff <- matrix(nrow = q, ncol = J)
    raneff_i <- vector(mode = "numeric", length = N)

    # EM parameters
    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    err <- control_EM_algorithm$err
    iter <- 0

    # I need to monitor the penalized likelihood
    loglik_pen <- loglik_pen_prev <- -.Machine$integer.max / 2
    loglik_pen_vec <- NULL

    crit <- TRUE
    group_indicator <- as.numeric(group)

    while (crit) {

      res_fixed <- y - X %*% beta
      est_second_moment <- 0

      # First E step ------------------------------------------------------------------
      # Compute the BLURP
      for (j in 1:J) {
        # iterate over different groups
        rows_j <- which(group_indicator == j)
        Z_j <- Z[rows_j, , drop = FALSE]
        res_fixed_j <- res_fixed[rows_j, drop = FALSE]
        Gamma_j <- solve(t(Z_j) %*% Z_j / sigma2 + solve(Omega))
        mu_j <- (Gamma_j %*% t(Z_j) %*% res_fixed_j) / sigma2
        mu_raneff[, j] <- mu_j
        raneff_i[rows_j] <- Z_j %*% mu_j
        # est_second_moment <-
        #   est_second_moment + Gamma_j + mu_j %*% t(mu_j)
      }

      # First M step ------------------------------------------------------------------
      # Compute the penalized betas

      penalized_regression <-
        glmnet::glmnet(
          x = X_no_intercept,
          y = y - raneff_i,
          family = "gaussian",
          standardize = FALSE,
          alpha = 1,
          lambda = lambda
        )

      beta <-
        as.vector(stats::coef(penalized_regression)) # I include the intercept in the set of estimated parameters

      # Second E step ------------------------------------------------------------------
      # Compute the BLURP and the estimated second moments
      for (j in 1:J) {
        # iterate over different groups
        rows_j <- which(group_indicator == j)
        Z_j <- Z[rows_j, , drop = FALSE]
        res_fixed_j <- res_fixed[rows_j, drop = FALSE]
        Gamma_j <- solve(t(Z_j) %*% Z_j / sigma2 + solve(Omega))
        mu_j <- (Gamma_j %*% t(Z_j) %*% res_fixed_j) / sigma2
        mu_raneff[, j] <- mu_j
        raneff_i[rows_j] <- Z_j %*% mu_j
        est_second_moment <-
          est_second_moment + Gamma_j + mu_j %*% t(mu_j)
      }

      # Second M step -----------------------------------------------------------
      # Compute the variance parameters

      Omega <- as.matrix(est_second_moment / J)
      sigma2 <- mean(y * (y - X %*% beta - raneff_i))

      #### log lik evaluation-------------------------------------------------

      loglik <- 0

      for (j in 1:J) {
        rows_j <- which(group_indicator == j)
        Z_j <- Z[rows_j, , drop = FALSE]
        y_j <- y[rows_j, drop = FALSE]
        X_j <- X[rows_j, , drop = FALSE]
        G_j <-
          Z_j %*% Omega %*% t(Z_j) + diag(sigma2, nrow = length(rows_j))
        loglik <-
          loglik + mvtnorm::dmvnorm(
            x = y_j,
            mean = c(X_j %*% beta),
            sigma = G_j,
            log = TRUE
          )
      }

      loglik_pen <-
        loglik - lambda * sum(abs(beta)) # objective function

      # check convergence
      err <-
        abs(loglik_pen - loglik_pen_prev) / (1 + abs(loglik_pen))
      loglik_pen_prev <- loglik_pen
      loglik_pen_vec <- c(loglik_pen_vec, loglik_pen)
      iter <- iter + 1
      crit <- (err > tol & iter < itermax)
    }

    return(
      list(
        beta = beta,
        Omega = Omega,
        sigma2 = sigma2,
        mu_raneff = mu_raneff,
        loglik = loglik,
        loglik_pen = loglik_pen,
        loglik_pen_trace = loglik_pen_vec,
        lambda = lambda,
        lambda_range = lambda_range
      )
    )
  }