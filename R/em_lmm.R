
# Standard linear mixed model fitted using an EM algorithm with both univariate and multivariate responses ----------------

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

# ecm_lmm_nocpp <-
#   function(X,
#            y,
#            Z,
#            group,
#            control_EM_algorithm = control_EM()) {
#
#     X <-  data.matrix(X)
#     y <-  data.matrix(y)
#     Z <-  data.matrix(Z)
#     q <- ncol(Z) # # ran eff
#     p <- ncol(X) # # fix eff
#     N <- nrow(X)
#     J <- length(unique(group))
#
#     # Set initial values
#     beta <- (stats::lm.fit(y = y, x = X))$coefficients
#     sigma2 <- 1
#     Omega <- diag(q)
#
#     # EM parameters
#     itermax <- control_EM_algorithm$itermax
#     tol <- control_EM_algorithm$tol
#     err <- control_EM_algorithm$err
#
#     iter <- 0
#     loglik <- loglik_prev <- -.Machine$integer.max / 2
#     loglik_vec <- NULL
#
#     crit <- TRUE
#     group_indicator <- as.numeric(group)
#
#     while (crit) {
#       res_fixed <- y - X %*% beta
#       est_second_moment <- 0
#
#       # E step ------------------------------------------------------------------
#
#       mu_raneff <- matrix(nrow = q, ncol = J)
#       est_second_moment_error <- 0
#       est_second_moment <- 0
#       raneff_i <- numeric(N)
#       for (j in 1:J) {
#         # iterate over different groups
#         rows_j <- which(group_indicator == j)
#         Z_j <- Z[rows_j, , drop = FALSE]
#         res_fixed_j <- res_fixed[rows_j, drop = FALSE]
#         Gamma_j <- solve(t(Z_j) %*% Z_j / sigma2 + solve(Omega))
#         mu_j <- (Gamma_j %*% t(Z_j) %*% res_fixed_j) / sigma2
#         mu_raneff[, j] <- mu_j
#         raneff_i[rows_j] <- Z_j %*% mu_j
#         est_second_moment <-
#           est_second_moment + Gamma_j + mu_j %*% t(mu_j)
#         est_second_moment_error <- est_second_moment_error + sum(diag(Z_j%*%Gamma_j%*%t(Z_j))) # second piece A.1 Rohart 2014
#       }
#
#       # M step ------------------------------------------------------------------
#
#       beta <-
#         (stats::lm.fit(y = y - raneff_i, x = X))$coefficients
#       Omega <- as.matrix(est_second_moment / J)
#       sigma2 <- mean(y * (y - X %*% beta - raneff_i))
#
#       #### log lik evaluation-------------------------------------------------
#
#       loglik <- 0
#
#       for (j in 1:J) {
#         rows_j <- which(group_indicator == j)
#         Z_j <- Z[rows_j, , drop = FALSE]
#         y_j <- y[rows_j, drop = FALSE]
#         X_j <- X[rows_j, , drop = FALSE]
#         G_j <-
#           Z_j %*% Omega %*% t(Z_j) + diag(sigma2, nrow = length(rows_j))
#         loglik <-
#           loglik + mvtnorm::dmvnorm(
#             x = y_j,
#             mean = c(X_j %*% beta),
#             sigma = G_j,
#             log = TRUE
#           )
#       }
#
#       # check convergence
#       err <-
#         abs(loglik - loglik_prev) / (1 + abs(loglik))
#       loglik_prev <- loglik
#       loglik_vec <- c(loglik_vec, loglik)
#       iter <- iter + 1
#       crit <- (err > tol & iter < itermax)
#     }
#
#     return(
#       list(
#         beta = beta,
#         Omega = Omega,
#         sigma2 = sigma2,
#         mu_raneff = mu_raneff,
#         loglik = loglik,
#         loglik_trace = loglik_vec
#       )
#     )
#   }

# Y is a matrix of responses, with r the number of response variables
#' @export
ecm_mlmm <-
  function(X,
           Y,
           Z,
           group,
           control_EM_algorithm = control_EM()) {

    X <-  data.matrix(X)
    Y <-  data.matrix(Y)
    Z <-  data.matrix(Z)

    q <- ncol(Z) # # ran eff
    p <- ncol(X) # # fix eff
    r <- ncol(Y) # # of responses
    N <- nrow(X)
    J <- length(unique(group))

    # Set initial values
    BETA <- (stats::lm.fit(y = Y, x = X))$coefficients
    SIGMA <- diag(r)
    PSI <- diag(q*r)

    # Objects to which I apply the vec operator and other useful quantities
    vec_XB <- c(X %*% BETA) # Nr x 1
    vec_Y <- c(Y) # Nr x 1
    I_r <- diag(r)

    # EM parameters
    itermax <- control_EM_algorithm$itermax
    tol <- control_EM_algorithm$tol
    err <- control_EM_algorithm$err

    iter <- 0
    loglik <- loglik_prev <- -.Machine$integer.max / 2
    loglik_vec <- NULL

    crit <- TRUE
    group_indicator <- as.numeric(group)
    vec_group_indicator <- rep(group_indicator,r)

    while (crit) {


      vec_res_fixed <- matrix((vec_Y - vec_XB),ncol = 1)
      PSI_inv <- solve(PSI)
      SIGMA_inv <- solve(SIGMA)

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

      # # Original R code
      # # est_second_moment_error <- matrix(0,nrow = r, ncol = r)
      # est_second_moment <- matrix(0,nrow = q*r, ncol = q*r)
      # mu_raneff <- array(dim=c(q,r,J))
      # vec_raneff_i <- vector(mode = "numeric", length = length(vec_Y))
      #
      #
      # for (j in 1:J) {
      #   # iterate over different groups
      #   rows_j <- which(group_indicator == j)
      #   vec_rows_j <- which(vec_group_indicator == j)
      #   I_nj <- diag(length(rows_j))
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
      #     est_second_moment + Gamma_j + vec_DELTA_j %*% t(vec_DELTA_j) # FIXME check this part
      #   # VAR_vec_E_j <- (I_r%x%Z_j)%*%Gamma_j%*%(I_r%x%t(Z_j))
      #
      #   # for (s in 1:r) {
      #   #   est_second_moment_error[s,s] <-
      #   # }
      #   # for (s in 1:r) {
      #   #   for (t in 1:r) {
      #   #     est_second_moment_error[s, t] <-
      #   #       est_second_moment_error[s, t] + #sum(diag(Z_j %*% Gamma_j %*% t(Z_j))) # FIXME
      #   #   }
      #   }

      raneff_i <- matrix(vec_raneff_i, nrow = N, ncol = r)

      # M step ------------------------------------------------------------------

      BETA <-
        (stats::lm.fit(y = Y - raneff_i, x = X))$coefficients

      PSI <- as.matrix(est_second_moment / J)

      SIGMA <- stats::cov(Y - X %*% BETA - raneff_i)*(N-1)/N # FIXME potentially a piece is missing here
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

      # Original R code
      # loglik <- 0
      #
      #
      # for (j in 1:J) {
      #   rows_j <- which(group_indicator == j)
      #   vec_rows_j <- which(vec_group_indicator == j)
      #   I_nj <- diag(length(rows_j))
      #   vec_XB_j <- vec_XB[vec_rows_j]
      #   Z_j <- Z[rows_j, , drop = FALSE]
      #   vec_Y_j <- vec_Y[vec_rows_j]
      #   G_j <-
      #     (I_r %x% Z_j)%*%PSI%*%(I_r %x% t(Z_j)) + SIGMA%x%I_nj
      #   loglik <-
      #     loglik + mvtnorm::dmvnorm(
      #       x = vec_Y_j,
      #       mean = vec_XB_j,
      #       sigma = G_j,
      #       log = TRUE
      #     )
      # }

      # check convergence
      err <-
        abs(loglik - loglik_prev) / (1 + abs(loglik))
      loglik_prev <- loglik
      loglik_vec <- c(loglik_vec, loglik)
      iter <- iter + 1
      crit <- (err > tol & iter < itermax)
    }

    mu_raneff <- array(c(e_step_lmm$mat_mu_raneff),dim=c(q,r,J)) # FIXME check when q>1

    return(
      list(
        BETA = BETA,
        PSI = PSI,
        SIGMA = SIGMA,
        mu_raneff = mu_raneff,
        loglik = loglik,
        loglik_trace = loglik_vec
      )
    )
  }
