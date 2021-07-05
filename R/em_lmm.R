
# Standard linear mixed model fitted using an EM algorithm ----------------
# It provides the same results of lme4::lmer()

em_lmm <-
  function(X,
           y,
           Z,
           group) {
    X <-  data.matrix(X)
    y <-  data.matrix(y)
    Z <-  data.matrix(Z)
    q <- ncol(Z) # # ran eff
    p <- ncol(X) # # fix eff
    N <- nrow(X)
    J <- length(unique(group))

    # Set initial values
    beta <- (lm.fit(y = y, x = X))$coefficients
    sigma2 <- 1
    Omega <- diag(q)

    mu_raneff <- matrix(nrow = q, ncol = J)
    raneff_i <- vector(mode = "numeric", length = N)

    # EM parameters
    itermax <- 1000
    tol <- 1e-8

    iter <- 0
    loglik <- loglik_prev <- -.Machine$integer.max / 2
    loglik_vec <- NULL

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

      beta <-
        (lm.fit(y = y - raneff_i, x = X))$coefficients # insert penalty term here
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

      # check convergence
      err <-
        abs(loglik - loglik_prev) / (1 + abs(loglik))
      loglik_prev <- loglik
      loglik_vec <- c(loglik_vec, loglik)
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
        loglik_trace = loglik_vec
      )
    )
  }
