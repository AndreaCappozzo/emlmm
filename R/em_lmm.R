
# Standard linear mixed model fitted using an EM algorithm ----------------
# It provides the same results of lme4::lmer()
#' @export
ecm_lmm <-
  function(X,
           y,
           Z,
           group,
           control_EM_algorithm = control_EM()) {

    X <-  data.matrix(X)
    y <-  data.matrix(y)
    Z <-  data.matrix(Z)
    q <- ncol(Z) # # ran eff
    p <- ncol(X) # # fix eff
    N <- nrow(X)
    J <- length(unique(group))

    # Set initial values
    beta <- (stats::lm.fit(y = y, x = X))$coefficients
    sigma2 <- 1
    Omega <- diag(q)

    # EM parameters
    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    err <- control_EM_algorithm$err

    iter <- 0
    loglik <- loglik_prev <- -.Machine$integer.max / 2
    loglik_vec <- NULL

    crit <- TRUE
    group_indicator <- as.numeric(group)

    while (crit) {
      res_fixed <- y - X %*% beta
      est_second_moment <- 0

      # E step ------------------------------------------------------------------

      e_step_lmm <- estep_lmm_cpp(
        res_fixed = res_fixed,
        Z = Z,
        group_indicator = group_indicator,
        sigma2 = sigma2,inv_Omega = solve(Omega),
        J = J
      )

      raneff_i <- e_step_lmm$raneff_i
      est_second_moment <- e_step_lmm$est_second_moment

      # M step ------------------------------------------------------------------

      beta <-
        (stats::lm.fit(y = y - raneff_i, x = X))$coefficients
      Omega <- as.matrix(est_second_moment / J)
      sigma2 <- mean(y * (y - X %*% beta - raneff_i))

      #### log lik evaluation-------------------------------------------------

      loglik <- log_lik_lmm_cpp(
        y = y,
        Z = Z,
        X = X,
        group_indicator = group_indicator,
        beta = beta,
        Omega = Omega,
        sigma2 = sigma2,
        J = J
      )

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
        mu_raneff = e_step_lmm$mu_raneff,
        loglik = loglik,
        loglik_trace = loglik_vec
      )
    )
  }
