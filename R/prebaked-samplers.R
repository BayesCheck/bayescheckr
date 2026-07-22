#Prebaked Samplers for known models

#to do next: implementation to have default workflows you can select for run_geweke and run_sbc e.g. an optional parameter like model_type = "own" (default, add in custom samplers) or "iid_normal", "linreg", "revjump", and so on

#iid Normal -------

#' iidnorm_prior_sampler
#'
#' Draws a sample from the conjugate Normal-Inverse-Gamma prior for the
#' IID Normal model.
#'
#' @param prior_hyper_params List containing the prior hyperparameters.
#' @return Returns a named list containing the sampled mean and variance.
#' @export
iidnorm_prior_sampler <- function(prior_hyper_params) {
  a0  <- prior_hyper_params$a0
  b0  <- prior_hyper_params$b0
  mu0 <- prior_hyper_params$mu0
  v20 <- prior_hyper_params$v20
  list(
    mu      = rnorm(1, mu0, sqrt(v20)),
    sigmasq = 1 / rgamma(1, a0, rate = b0)
  )
}

#' iidnorm_likelihood_sampler
#'
#' Simulates IID Normal observations using the supplied parameter values.
#'
#' @param n_obs Number of observations to generate.
#' @param theta Named list containing the model parameters.
#' @param prior_hyper_params Prior hyperparameters (unused).
#' @return Returns a numeric vector of simulated observations.
#' @export
iidnorm_likelihood_sampler <- function(n_obs, theta, prior_hyper_params) {
  rnorm(n_obs, theta$mu, sqrt(theta$sigmasq))
}

#' iidnorm_posterior_sampler
#'
#' Runs a Gibbs sampler for the IID Normal model with conjugate priors.
#'
#' @param n_draws Number of posterior draws.
#' @param y Observed data.
#' @param prior_hyper_params List containing the prior hyperparameters.
#' @param init Initial values for the Markov chain.
#' @return Returns a matrix of posterior draws.
#' @export
iidnorm_posterior_sampler <- function(n_draws, y, prior_hyper_params,
                                      init = c(mean(y), var(y))) {
  n <- length(y); sum.y <- sum(y)
  a0  <- prior_hyper_params$a0
  b0  <- prior_hyper_params$b0
  mu0 <- prior_hyper_params$mu0
  v20 <- prior_hyper_params$v20

  THETA <- matrix(0, nrow = n_draws, ncol = 2)

  theta <- c(mean(y), var(y))

  for (s in 1:n_draws) {
    v2n <- 1 / (1/v20 + n/theta[2])
    mun <- v2n * (mu0/v20 + sum.y/theta[2])
    theta[1] <- rnorm(1, mun, sqrt(v2n))

    an <- a0 + n/2
    bn <- b0 + sum((y - theta[1])^2) / 2
    theta[2] <- 1 / rgamma(1, shape = an, rate = bn)

    THETA[s, ] <- theta
  }

  return(THETA)
  # need to convert matrix -> list of named theta-lists (new contract)

}

#Simple linear regression model -------

#' linreg_prior_sampler
#'
#' Draws regression coefficients and the error variance from the Bayesian
#' linear regression prior.
#'
#' @param prior_hyper_params List containing the prior hyperparameters.
#' @return Returns a named list containing the sampled coefficients and variance.
#' @export
linreg_prior_sampler <- function(prior_hyper_params){

  beta_0 <- prior_hyper_params$beta_0
  v_0    <- prior_hyper_params$v_0
  a_0    <- prior_hyper_params$a_0
  b_0    <- prior_hyper_params$b_0

  #simulate betas
  beta <- as.vector(mvtnorm::rmvnorm(1, beta_0, v_0))

  #simulate sigma squared
  sig2 <- 1 / rgamma(1, shape = a_0, rate = b_0)

  beta_list <- setNames(as.list(beta), paste0("beta", seq_along(beta)))
  return(c(beta_list, list(sigmasq = sig2)))
}

#' linreg_likelihood_sampler
#'
#' Simulates responses from a Bayesian linear regression model.
#'
#' @param n_obs Number of observations to generate.
#' @param theta Named list containing the regression coefficients and variance.
#' @param prior_hyper_params List containing the design matrix and prior hyperparameters.
#' @return Returns a numeric vector of simulated responses.
#' @export
linreg_likelihood_sampler <- function(n_obs, theta, prior_hyper_params){

  X <- prior_hyper_params$X

  p <- length(theta) - 1

  beta <- unlist(theta[1:p])
  sig2 <- theta[[p + 1]]

  y <- as.numeric(X %*% beta) + rnorm(n_obs, mean = 0, sd = sqrt(sig2))

  return(y)
}

#' linreg_posterior_sampler
#'
#' Runs a Gibbs sampler for Bayesian linear regression with conjugate priors.
#'
#' @param n_draws Number of posterior draws.
#' @param y Observed response vector.
#' @param prior_hyper_params List containing the prior hyperparameters.
#' @param init Initial values for the Markov chain.
#' @return Returns a matrix of posterior draws.
#' @export
linreg_posterior_sampler <- function(n_draws, y, prior_hyper_params,
                              init = NULL)
{

  # unpack everything
  X      <- prior_hyper_params$X
  beta_0 <- prior_hyper_params$beta_0
  v_0    <- prior_hyper_params$v_0
  a_0    <- prior_hyper_params$a_0
  b_0    <- prior_hyper_params$b_0

  #pre-allocate storage and initialize
  p <- ncol(X)
  n_obs <- length(y)
  THETA <- matrix(0, nrow = n_draws, ncol = p+1)
  #colnames(THETA) <- c("beta1", "beta2", "sigmasq")
  colnames(THETA) <- c(paste0("beta", 1:p), "sigmasq")

  #initialize values

  if (is.null(init)) {
    theta <- c(rep(0, p), 1)
  } else {
    theta <- unlist(init)
  }

  v_0_inv <- solve(v_0)
  X_trans <- t(X)

  for(s in 1:n_draws){

    beta_n <- solve( v_0_inv + (X_trans%*%X)/theta[p+1] ) %*%
      ( v_0_inv%*%beta_0 + (X_trans%*%y)/theta[p+1] )

    v_n <- solve( v_0_inv + (X_trans%*%X)/theta[p+1] )

    theta[1:p] <- mvtnorm::rmvnorm(1, beta_n, v_n)

    #draw from updated sigma
    ssr <- t(y - X%*%theta[1:p])%*%(y - X%*%theta[1:p])

    a_n <- a_0 + n_obs/2
    b_n <- ssr/2 + b_0

    theta[p+1] <- 1 / rgamma(1, shape = a_n, rate = b_n)

    THETA[s, ] <- theta

  }

  return(THETA)
}

#Reversible jump MCMC, Andrieu & Doucet (1999) paper -----

#' doucet_build_D
#'
#' Constructs the Fourier design matrix used by the Andrieu and Doucet
#' reversible-jump MCMC example.
#'
#' @return Returns the Fourier design matrix.
#' @export
doucet_build_D <- function(omega_k, n) {
  k <- length(omega_k)

  i_idx <- 0:(n-1)

  if (k == 0) {
    return(matrix(0, nrow = n, ncol = 0))
  }

  D <- matrix(0, nrow = length(i_idx), ncol = 2*k)
  for (j in seq_len(k)) {
    D[, 2*j - 1] <- cos(omega_k[j]*i_idx)
    D[, 2*j]     <- sin(omega_k[j]*i_idx)
  }
  return(D)
}

#' doucet_yPy
#'
#' Computes the quadratic form used in the Andrieu and Doucet reversible-jump
#' MCMC algorithm.
#'
#' @return Returns the computed quadratic form.
#' @export
doucet_yPy <- function(omega_k, y, hyper) {
  k <- length(omega_k)
  if (k == 0) return(sum(y^2))
  D <- doucet_build_D(omega_k, length(y))
  DtD <- crossprod(D)
  Mk_inv <- ((hyper$delta2 + 1) / hyper$delta2) * DtD
  Dty <- crossprod(D, y)
  sol <- tryCatch(solve(Mk_inv, Dty), error = function(e) NULL)
  if (is.null(sol)) return(sum(y^2))
  sum(y^2) - as.numeric(crossprod(Dty, sol))
}

#' doucet_inner_update
#'
#' Performs the within-model Gibbs and Metropolis-Hastings updates for the
#' reversible-jump MCMC sampler.
#'
#' @return Returns the updated model parameters.
#' @export
doucet_inner_update <- function(omega, a, sigma2, y, prior_hyper_params) {
  hyper <- prior_hyper_params
  k <- length(omega)
  n_obs <- length(y)
  mh_lambda <- 0.3
  sigma_rw <- 1 / (5 * n_obs)

  # FFT proposal setup (same as before)
  fft_y <- fft(y)
  periodogram <- Mod(fft_y)^2
  Np <- n_obs
  periodogram <- periodogram[1:Np]
  freq_grid <- (0:(Np - 1)) * pi / Np
  period_probs <- periodogram / sum(periodogram)

  propose_fft_freq <- function() {
    bin <- sample.int(Np, size = 1, prob = period_probs)
    lo <- freq_grid[bin]
    hi <- lo + pi / Np
    runif(1, lo, min(hi, pi))
  }

  density_fft_freq <- function(w) {
    bin <- min(Np, max(1, floor(w / (pi / Np)) + 1))
    period_probs[bin] / (pi / Np)
  }

  # MH on frequencies (same as before)
  for (j in seq_len(k)) {
    omega_minus_j <- omega[-j]
    u <- runif(1)

    if (u < mh_lambda) {
      omega_prop <- propose_fft_freq()
      q_fwd <- density_fft_freq(omega_prop)
      q_bwd <- density_fft_freq(omega[j])
    } else {
      omega_prop <- rnorm(1, mean = omega[j], sd = sigma_rw)
      if (omega_prop <= 0 || omega_prop >= pi) next
      q_fwd <- dnorm(omega_prop, mean = omega[j], sd = sigma_rw)
      q_bwd <- dnorm(omega[j], mean = omega_prop, sd = sigma_rw)
    }

    omega_cand <- sort(c(omega_minus_j, omega_prop))
    omega_old <- sort(c(omega_minus_j, omega[j]))

    s_cand <- doucet_yPy(omega_cand, y, hyper)
    s_old <- doucet_yPy(omega_old, y, hyper)

    log_target_ratio <- ((n_obs + hyper$v0) / 2) *
      (log(hyper$gamma0 + s_old) - log(hyper$gamma0 + s_cand))
    log_q_ratio <- log(q_bwd) - log(q_fwd)
    log_alpha <- log_target_ratio + log_q_ratio

    if (log(runif(1)) < log_alpha) {
      omega <- omega_cand
    }
  }

  # Gibbs (same as before)
  if (k == 0) {
    Pky <- sum(y^2)
    shape_sigma2 <- (hyper$v0 + n_obs) / 2
    rate_sigma2 <- (hyper$gamma0 + Pky) / 2
    sigma2 <- 1 / rgamma(1, shape = shape_sigma2, rate = rate_sigma2)
    a <- numeric(0)
  } else {
    D <- doucet_build_D(omega, n_obs)
    DtD <- crossprod(D)
    Mk_inv <- ((hyper$delta2 + 1) / hyper$delta2) * DtD
    Mk <- solve(Mk_inv)
    m_k <- Mk %*% crossprod(D, y)
    P_k <- diag(n_obs) - D %*% Mk %*% t(D)
    Pky <- as.numeric(t(y) %*% P_k %*% y)

    shape_sigma2 <- (hyper$v0 + n_obs) / 2
    rate_sigma2 <- (hyper$gamma0 + Pky) / 2
    sigma2 <- 1 / rgamma(1, shape = shape_sigma2, rate = rate_sigma2)
    a <- as.vector(MASS::mvrnorm(1, mu = as.numeric(m_k),
                                 Sigma = sigma2 * Mk))
  }

  list(omega = omega, a = a, sigma2 = sigma2)
}

#end of Doucet helper functions

#' doucet_prior_sampler
#'
#' Draws parameters from the prior distribution for the Andrieu and Doucet
#' reversible-jump MCMC example.
#'
#' @param prior_hyper_params List containing the prior hyperparameters.
#' @return Returns a named vector containing the sampled model parameters.
#' @export
doucet_prior_sampler <- function(prior_hyper_params) {

  n <- prior_hyper_params$n
  k_max <- prior_hyper_params$k_max
  Lambda <- prior_hyper_params$Lambda
  delta2 <- prior_hyper_params$delta2
  v0 <- prior_hyper_params$v0
  gamma0 <- prior_hyper_params$gamma0

  # 1. sample k ~ truncated Poisson on {0, ..., k_max}
  repeat {
    k <- rpois(1, Lambda)
    if (k <= k_max) break
  }

  # 2. sample omega_k ~ iid Uniform(0, pi), k-length vector
  if (k > 0) {
    omega <- runif(k, 0, pi)
  } else {
    omega <- numeric(0)
  }

  # 3. sample sigma2 ~ IG(v0/2, gamma0/2)
  sigma2 <- 1 / rgamma(1, shape = (v0/2), rate = (gamma0/2))

  # 4. sample a_k | omega_k, sigma2, k ~ N_k(0, delta2 * sigma2 * (D^T D)^{-1})
  if (k > 0) {
    D <- doucet_build_D(omega, n = n)
    DtD_inv <- solve(t(D) %*% D)
    Sigma_a <- sigma2 * delta2 * DtD_inv #sigma2 value x Sigma matrix described
    a <- as.vector(mvtnorm::rmvnorm(1,
                                    mean = rep(0, 2*k), #mu = 0 vector
                                    sigma = Sigma_a)) #covariance matrix
  } else {
    a <- numeric(0) #empty (0) vector
  }

  # 5. pad to k_max and assemble named vector
  omega_padded <- rep(NA_real_, k_max)
  if (k > 0) omega_padded[1:k] <- omega

  a_padded <- rep(NA_real_, 2 * k_max)
  if (k > 0) a_padded[1:(2 * k)] <- a

  theta <- c(k = k,
             setNames(omega_padded, paste0("omega", 1:k_max)),
             setNames(a_padded, paste0("a", 1:(2 * k_max))),
             sigma2 = sigma2)

  return(theta)
}

#now to call SBC and parallelize, only need to wrap local({doucet_prior_sampler})

#' doucet_likelihood_sampler
#'
#' Simulates observations from the Fourier regression model used in the
#' Andrieu and Doucet example.
#'
#' @param n_obs Number of observations to generate.
#' @param theta Model parameters.
#' @param prior_hyper_params List containing the prior hyperparameters.
#' @return Returns a numeric vector of simulated observations.
#' @export
doucet_likelihood_sampler <- function(n_obs, theta, prior_hyper_params) {

  # infer k_max from vector length: 1(k) + k_max(omega) + 2*k_max(a) + 1(sigma2)
  k_max <- (length(theta) - 2) / 3

  k <- theta["k"]
  sigma2 <- theta["sigma2"]

  if (k > 0) {
    omega <- theta[paste0("omega", 1:k)]   # only first k slots, rest are NA
    a <- theta[paste0("a", 1:(2 * k))]  # only first 2k slots
    D <- doucet_build_D(omega, n_obs)
    mu <- as.vector(D %*% a)
  } else {
    mu <- rep(0, n_obs)
  }

  y <- rnorm(n_obs, mean = mu, sd = sqrt(sigma2))
  return(y)
}

#' doucet_posterior_sampler
#'
#' Runs the reversible-jump MCMC sampler described by Andrieu and Doucet (1999).
#'
#' @param n_draws Number of posterior draws.
#' @param y Observed data.
#' @param prior_hyper_params List containing the prior hyperparameters.
#' @param init Initial values for the Markov chain.
#' @param buggy Logical indicating whether to run the intentionally flawed sampler.
#' @return Returns a matrix of posterior draws.
#' @export
doucet_posterior_sampler <- function(n_draws,
                                     y,
                                     prior_hyper_params,
                                     init = NULL,
                                     buggy = FALSE) {
  hyper <- prior_hyper_params
  n_obs <- length(y)
  k_max <- hyper$k_max
  Lambda <- hyper$Lambda

  # Initialize from prior or use provided init
  if (is.null(init)) {
    init_state <- doucet_prior_sampler(prior_hyper_params)
  } else {
    init_state <- init
  }
  k <- init_state["k"]
  omega <- sort(init_state[grepl("^omega", names(init_state))])
  a <- init_state[grepl("^a", names(init_state))]
  sigma2 <- init_state["sigma2"]

  c_prob <- 0.5
  b <- function(k) {
    if (k == k_max) return(0)
    c_prob * min(1, dpois(k + 1, Lambda) / dpois(k, Lambda))
  }
  d <- function(k) {
    if (k == 0) return(0)
    c_prob * min(1, dpois(k - 1, Lambda) / dpois(k, Lambda))
  }

  col_names <- c("k", paste0("omega", 1:k_max),
                 paste0("a", 1:(2*k_max)),
                 "sigma2")
  out <- matrix(NA_real_, nrow = n_draws, ncol = length(col_names),
                dimnames = list(NULL, col_names))

  for (i in seq_len(n_draws)) {
    bk <- b(k)
    dk <- d(k)
    u <- runif(1)

    if (u < bk) {
      omega_new <- runif(1, 0, pi)
      omega_prop <- sort(c(omega, omega_new))

      s_old <- doucet_yPy(omega, y, hyper)
      s_prop <- doucet_yPy(omega_prop, y, hyper)

      r_birth <- ((hyper$gamma0 + s_old) / (hyper$gamma0 + s_prop))^((n_obs + hyper$v0) / 2) *
        1 / (1 + hyper$delta2)

      #This is the mistake! Acceptance ratio for birth moves was too small,
      #overly restricting the movement of the MCMC chain
      if (buggy) r_birth <- r_birth / (k + 1)

      if (runif(1) < min(1, r_birth)) {
        k <- k + 1
        omega <- omega_prop
        updated <- doucet_inner_update(omega, rep(0, 2*k), sigma2, y, hyper)
        a <- updated$a
        sigma2 <- updated$sigma2
      }

    } else if (u < bk + dk) {
      j_kill <- sample.int(k, 1)
      omega_prop <- omega[-j_kill]

      s_old <- doucet_yPy(omega, y, hyper)
      s_prop <- doucet_yPy(omega_prop, y, hyper)

      r_birth <- ((hyper$gamma0 + s_prop) / (hyper$gamma0 + s_old))^((n_obs + hyper$v0) / 2) *
        1 / (1 + hyper$delta2)

      if (runif(1) < min(1, 1 / r_birth)) {
        k <- k - 1
        omega <- omega_prop
        a <- if (k > 0) a[-c(2*j_kill - 1, 2*j_kill)] else numeric(0)
        updated <- doucet_inner_update(omega, a, sigma2, y, hyper)
        a <- updated$a
        sigma2 <- updated$sigma2
      }

    } else {
      updated <- doucet_inner_update(omega, a, sigma2, y, hyper)
      omega <- updated$omega
      a <- updated$a
      sigma2 <- updated$sigma2
    }

    omega_pad <- rep(NA_real_, k_max)
    if (k > 0) omega_pad[1:k] <- omega

    a_pad <- rep(NA_real_, 2 * k_max)
    if (k > 0) a_pad[1:(2*k)] <- a

    out[i, ] <- c(k, omega_pad, a_pad, sigma2)
  }

  return(out)
}

#Assembly error: Primiceri (2005) ------

#' primiceri_prior_sampler
#'
#' Draws parameters from the prior distribution for the Primiceri (2005)
#' stochastic volatility example.
#'
#' @param prior_hyper_params List containing the prior hyperparameters.
#' @return Returns the sampled model parameters.
#' @export
primiceri_prior_sampler <- function(prior_hyper_params){

  theta_bar_0 <- prior_hyper_params["theta_bar_0"]
  v_0         <- prior_hyper_params["v_0"]
  a_0         <- prior_hyper_params["a_0"]
  b_0         <- prior_hyper_params["b_0"]
  n           <- prior_hyper_params["n"]

  theta  <- rnorm(1, mean = theta_bar_0, sd = sqrt(v_0))
  u      <- rnorm(1, mean = a_0, sd = sqrt(b_0))
  sigmasq <- exp(u)
  s      <- sample(1:7, size = n, replace = TRUE, prob = ksc_pi)

  c(theta = theta, sigmasq = sigmasq)

}

#' primiceri_likelihood_sampler
#'
#' Simulates observations from the Primiceri stochastic volatility example.
#'
#' @param n_obs Number of observations to generate.
#' @param theta Model parameters.
#' @param prior_hyper_params List containing the prior hyperparameters.
#' @return Returns a numeric vector of simulated observations.
#' @export
primiceri_likelihood_sampler <- function(n_obs, theta, prior_hyper_params){
  rnorm(n_obs, mean = theta[["theta"]], sd = sqrt(theta[["sigmasq"]]))
}

#' primiceri_posterior_sampler
#'
#' Runs the posterior sampler for the Primiceri (2005) stochastic volatility
#' example.
#'
#' @param n_draws Number of posterior draws.
#' @param y Observed data.
#' @param prior_hyper_params List containing the prior hyperparameters.
#' @param init Initial values for the Markov chain.
#' @param buggy Logical indicating whether to run the intentionally flawed sampler.
#' @return Returns a matrix of posterior draws.
#' @export
primiceri_posterior_sampler <- function(n_draws,
                              y,
                              prior_hyper_params,
                              init = NULL,
                              buggy = FALSE){

  # extract hyperparameters
  theta_bar_0 <- prior_hyper_params["theta_bar_0"]
  v_0         <- prior_hyper_params["v_0"]
  a_0         <- prior_hyper_params["a_0"]
  b_0         <- prior_hyper_params["b_0"]
  n           <- length(y)

  # preallocate storage
  THETA <- matrix(0, nrow = n_draws, ncol = 2 + n)
  colnames(THETA) <- c("theta", "sigmasq", paste0("s_", 1:n))

  # initialize from init if provided, otherwise use data-based defaults
  if (!is.null(init)) {
    theta_cur <- as.numeric(init["theta"])
    u_cur     <- log(as.numeric(init["sigmasq"]))
    s_names   <- grep("^s_", names(init), value = TRUE)
    s_cur     <- if (length(s_names) == n) as.numeric(init[s_names]) else rep(5, n)
  } else {
    theta_cur <- mean(y)
    u_cur     <- log(var(y))
    s_cur     <- rep(5, n)
  }

  for (i in 1:n_draws) {

    # step 1: draw u | theta, s, y
    z     <- log((y - theta_cur)^2 + 1e-20)
    z_tilde <- z - (ksc_mu[s_cur] - 1.2704)   # mean-corrected observations
    b_n   <- 1 / (1/b_0 + sum(1/ksc_tau[s_cur]))
    a_n   <- b_n * (a_0/b_0 + sum(z_tilde / ksc_tau[s_cur]))
    u_prop <- rnorm(1, mean = a_n, sd = sqrt(b_n))

    # true log-likelihood at proposed and current u
    loglik_true <- function(u) {
      sum(dnorm(y, mean = theta_cur, sd = sqrt(exp(u)), log = TRUE))
    }

    # mixture approximation log-likelihood at proposed and current u
    loglik_mix <- function(u) {
      # for each observation i, compute log(sum_j pi_j * N(z_i; u + m_j - 1.2704, v_j))
      log_mix_i <- sapply(1:n, function(ii) {
        log(sum(ksc_pi * dnorm(z[ii], mean = u + ksc_mu - 1.2704, sd = sqrt(ksc_tau))))
      })
      sum(log_mix_i)
    }

    # log acceptance ratio
    log_alpha <- (loglik_true(u_prop) - loglik_true(u_cur)) -
      (loglik_mix(u_prop)  - loglik_mix(u_cur))

    # MH step
    if (log(runif(1)) < log_alpha) {
      u_cur <- u_prop
    }

    if (buggy) {
      # step 2: draw s | u, theta, y
      z <- log((y - theta_cur)^2 + 1e-20)  # recompute with OUTDATED theta_cur
      s_cur <- sapply(1:n, function(ii) {
        log_w <- log(ksc_pi) + dnorm(z[ii],
                                     mean = u_cur + ksc_mu - 1.2704,
                                     sd   = sqrt(ksc_tau),
                                     log  = TRUE)
        log_w <- log_w - max(log_w)
        w <- exp(log_w) / sum(exp(log_w))
        sample(1:7, size = 1, prob = w)
      })

      # step 3: draw theta | u, y
      # not conditioned on s -> mistake
      sigmasq_cur <- exp(u_cur)
      v_n         <- 1 / (1/v_0 + n/sigmasq_cur)
      theta_n     <- v_n * (theta_bar_0/v_0 + sum(y)/sigmasq_cur)
      theta_cur   <- rnorm(1, mean = theta_n, sd = sqrt(v_n))
    }

    else {
      # step 2: draw theta | u, y
      sigmasq_cur <- exp(u_cur)
      v_n         <- 1 / (1/v_0 + n/sigmasq_cur)
      theta_n     <- v_n * (theta_bar_0/v_0 + sum(y)/sigmasq_cur)
      theta_cur   <- rnorm(1, mean = theta_n, sd = sqrt(v_n))

      # step 3: draw s | u, y
      z <- log((y - theta_cur)^2 + 1e-20)  # recompute with updated theta_cur
      s_cur <- sapply(1:n, function(ii) {
        log_w <- log(ksc_pi) + dnorm(z[ii],
                                     mean = u_cur + ksc_mu - 1.2704,
                                     sd   = sqrt(ksc_tau),
                                     log  = TRUE)
        log_w <- log_w - max(log_w)
        w <- exp(log_w) / sum(exp(log_w))
        sample(1:7, size = 1, prob = w)
      })
    }

    # store
    THETA[i, "theta"]            <- theta_cur
    THETA[i, "sigmasq"]          <- exp(u_cur)
    THETA[i, paste0("s_", 1:n)]  <- s_cur
  }

  return(THETA)
}

.make_sampler_registry <- function() {
  list(
    prior = list(
      iidnorm = iidnorm_prior_sampler,
      linreg  = linreg_prior_sampler,
      doucet  = doucet_prior_sampler,
      primiceri = primiceri_prior_sampler
    ),
    likelihood = list(
      iidnorm = iidnorm_likelihood_sampler,
      linreg  = linreg_likelihood_sampler,
      doucet  = doucet_likelihood_sampler,
      primiceri = primiceri_likelihood_sampler
    ),
    posterior = list(
      iidnorm = iidnorm_posterior_sampler,
      linreg  = linreg_posterior_sampler,
      doucet  = doucet_posterior_sampler
    ),
    default_hyper = list(
      iidnorm = list(
        a0  = 1, b0  = 1, mu0 = 0, v20 = 1
      ),
      linreg = list(
        beta_0 = c(0, 0),
        v_0    = diag(2),
        a_0    = 1,
        b_0    = 1,
        X      = matrix(rnorm(10 * 2), nrow = 10, ncol = 2)
      ),
      doucet = list(
        n            = 64,
        k_max        = 3,
        Lambda       = 3,
        delta2       = 100,
        v0           = 2,
        gamma0       = 2,
        lambda       = 0.3,
        sigma_rw     = 1 / (5 * 64),
        n_mcmc_sweeps = 50
      )
    )
  )
}
