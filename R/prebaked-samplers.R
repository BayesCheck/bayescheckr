#Prebaked Samplers for known models

#to do next: implementation to have default workflows you can select for run_geweke and run_sbc e.g. an optional parameter like model_type = "own" (default, add in custom samplers) or "iid_normal", "linreg", "revjump", and so on

#iid Normal -------

#' @export
iidnorm_prior_sampler <- function(prior_hyper_params) {
  a0 <- prior_hyper_params$a0
  b0 <- prior_hyper_params$b0
  mu0 <- prior_hyper_params$mu0
  v20 <- prior_hyper_params$v20

  # simulate mu ~ N(mu0, v20)
  mu <- rnorm(1, mu0, sqrt(v20))

  # independently simulate sigma2 ~ IG(a0, b0)
  sigma2 <- 1 / rgamma(1, a0, rate = b0)

  return(c(mu = mu, sigma2 = sigma2))
}

#' @export
iidnorm_likelihood_sampler <- function(n_obs, theta) {
  mu <- theta[1]
  sigma2 <- theta[2]

  y <- rnorm(n_obs, mu, sqrt(sigma2))

  return(y)
}

#' @export
iidnorm_posterior_sampler <- function(ndraws,
                                      y,
                                      prior_hyper_params,
                                      init = c(mean(y), var(y))
                                      ) {
  a0 <- prior_hyper_params$a0
  b0 <- prior_hyper_params$b0
  mu0 <- prior_hyper_params$mu0
  v20 <- prior_hyper_params$v20
  n <- length(y)

  # preallocate storage
  THETA <- matrix(0, nrow = ndraws, ncol = length(init))

  # initialize
  theta <- init

  for(s in 1:ndraws){

    # draw mean given variance and data

    v2n = 1 / (1/v20 + n/theta[2])
    mun = v2n * (mu0/v20 + sum(y)/theta[2])

    theta[1] <- rnorm(1, mun, sqrt(v2n))

    # draw variance given mean and data

    an = a0 + n/2
    bn = b0 + sum((y - theta[1])^2) / 2

    theta[2] <- 1 / rgamma(1, shape = an, rate = bn)

    # store the latest draw in the sth row of the THETA matrix

    THETA[s, ] <- theta

  }

  colnames(THETA) <- c("mu", "sigma2")
  return(THETA)
}

#Simple linear regression model -------

#' @export
linreg_prior_sampler <- function(prior_hyper_params) {
  beta_0 <- prior_hyper_params$beta_0
  v_0    <- prior_hyper_params$v_0
  a_0    <- prior_hyper_params$a_0
  b_0    <- prior_hyper_params$b_0

  #simulate betas
  beta <- as.vector(mvtnorm::rmvnorm(1, beta_0, v_0))

  #simulate sigma squared
  sig2 <- 1 / rgamma(1, shape = a_0, rate = b_0)

  return(c(beta1 = beta[1], beta2 = beta[2], sigmasq = sig2))
}

#' @export
linreg_likelihood_sampler <- function(n_obs, theta) {
  X <- prior_hyper_params$X
  #p <- ncol(X)
  p <- length(theta) - 1
  beta <- theta[1:p]
  sig2 <- theta[p+1]

  y <- as.numeric(X %*% beta) + rnorm(n_obs, mean = 0, sd = sqrt(sig2))

  return(y)
}

#' @export
linreg_posterior_sampler <- function(n_draws,
                                     y,
                                     prior_hyper_params,
                                     init = c(rep(0, p), 1)
                                     ) {
  X      <- prior_hyper_params$X
  beta_0 <- prior_hyper_params$beta_0
  v_0    <- prior_hyper_params$v_0
  a_0    <- prior_hyper_params$a_0
  b_0    <- prior_hyper_params$b_0

  #pre-allocate storage and initialize
  p <- ncol(X)
  n_obs <- length(y)
  THETA <- matrix(0, nrow = ndraws, ncol = p+1)
  colnames(THETA) <- c("beta1", "beta2", "sigmasq")

  #initialize values
  theta <- init

  v_0_inv <- solve(v_0)
  X_trans <- t(X)

  for(s in 1:ndraws){

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

# PENDING/TO DO NEXT: Doucet & Primiceri correct implementations

#Reversible jump MCMC, Andrieu & Doucet (1999) paper -----

#first, define all Doucet helper functions so they're exported
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
    D <- build_D(omega, n = n)
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

#' @export
doucet_likelihood_sampler <- function(n_obs, theta) {

  # infer k_max from vector length: 1(k) + k_max(omega) + 2*k_max(a) + 1(sigma2)
  k_max <- (length(theta) - 2) / 3

  k <- theta["k"]
  sigma2 <- theta["sigma2"]

  if (k > 0) {
    omega <- theta[paste0("omega", 1:k)]   # only first k slots, rest are NA
    a <- theta[paste0("a", 1:(2 * k))]  # only first 2k slots
    D <- build_D(omega, n_obs)
    mu <- as.vector(D %*% a)
  } else {
    mu <- rep(0, n_obs)
  }

  y <- rnorm(n_obs, mean = mu, sd = sqrt(sigma2))
  return(y)
}

#' @export
doucet_posterior_sampler <- function(n_draws,
                                     y,
                                     prior_hyper_params,
                                     init = NULL) {
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

  col_names <- c("k", paste0("omega_", 1:k_max),
                 paste0("a_c", 1:k_max), paste0("a_s", 1:k_max),
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

#' @export
primiceri_likelihood_sampler <- function(n_obs, theta){
  rnorm(n_obs, mean = theta[["theta"]], sd = sqrt(theta[["sigmasq"]]))
}

.make_sampler_registry <- function() {
  list(
    prior = list(
      iidnorm = iidnorm_prior_sampler,
      linreg  = linreg_prior_sampler,
      doucet  = doucet_prior_sampler
    ),
    likelihood = list(
      iidnorm = iidnorm_likelihood_sampler,
      linreg  = linreg_likelihood_sampler,
      doucet  = doucet_likelihood_sampler
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

# pending: Primiceri posterior (will add on Wednesday)

# want to add modifications/params to run_geweke & run_sbc

# IMPORTANT TODO: generic function that takes outputs from any classifier
# and computes the test statistic and p-value per the paper in Lopez-Paz & Oquab (2017)

# PSEUDOCODE TIME
#
# My vision is that run_geweke and run_sbc have fully automated workflows
# for some pre-fitted menu of models:
#   - iidnorm (for iid Normal)
#   - linreg (for linear regression)
#   - doucet (for the reversible-jump MCMC Andrieu & Doucet (1999) paper)
#   - primiceri (for the assembly error in Primiceri (2005))
#
# Then if you input one of these string arguments into the OPTIONAL param
# of model into run_geweke or run_sbc, you don't have to feed in any of the
# typical sampler arguments like posterior_sampler, prior_sampler, or anything else
# because it's already done for you. All that's needed are the other initializing
# arguments: n_obs, n_draws = 1e4, n_burn = 1000, n_thin = 5, prior_hyper_params,
# test_fns = NULL, parallelize = TRUE (note that most already have defaults
# as written - maybe I put defaults for all of them, even hyperparams, so that
# it's possible to run one of these workflows with JUST the model type selected
# and get the "vanilla" baseline version of results back).
#
# e.g., run_geweke(model = "doucet", prior_hyper_params = prior_hp_list [predefined by user], n_obs = 10)
# or even just run_geweke(model = "doucet") and the same goes for run_sbc.
#
#
# so like again taking run_geweke as my guinea pig here but this all applies to
# run_sbc as well:
#
#   prior_list <- list("iidnorm" = iidnorm_prior_sampler,
#                      "linreg" = linreg_prior_sampler)
#                      ... #and so on
#   likelihood_list <- list("iidnorm" = iidnorm_likelihood_sampler,
#                      "linreg" = linreg_likelihood_sampler)
#                      ... #and so on
#   posterior_list <- list("iidnorm" =  iidnorm_prior_sampler,
#                      "linreg" = linreg_posterior_sampler)
#                      ... #and so on
#
#   run_geweke <- function(model = NULL,
#                          prior_hyper_params = [],
#                          n_obs = [],
#                          n_sims = [], #for sbc
#                          ...
#                          ) { #maybe to start/keep things clean we still have user define prior_hyper_params and n_obs every time
#     if (is.null(model)) run_geweke(#take all the normal arguments from the user)
#
#     else:  #where model name string = index key thru which we find  relevant function in list of samplers
#       prior_sampler <- prior_list[[model]]
#       likelihood_sampler <- likelihood_list[[model]]
#       posterior_sampler <- posterior_list[[model]]
#
#       #after pulling all of our defaults run the same workflow
#
#       run_geweke(prior_sampler = prior_sampler,
#                  likelihood_sampler = likelihood_sampler,
#                  posterior_sampler = posterior_sampler,
#                  n_obs = n_obs,
#                  n_draws = n_draws, #if included in the params above in the ... / else leave blank for defaults
#                  #and so on for other parameters as specified, else use defaults like n_thin = 5
#                  )
#
#   }
#
#   One thing to flag: the SBC versions of samplers don't need an init argument because they're not doing
#   successive-conditional draws like Geweke needs to, but Geweke def needs that argument
#   to do all the init = theta stuff. So the samplers need to be robust to both possibilities
#   and thus include init as an argument. If this would cause a bug (i.e., calling init as a param = SBC bug)
#   then you'll have to call different versions of the prebaked samplers based on whether
#   the function called is Geweke or SBC. In which case there need to be 2 versions
#   of each list of samplers, one with the geweke draws, one with the sbc draws.
#   Also, if there are problems in calling the list of samplers as globally defined
#   objects, write the code to create these lists (which will call other
#   package functions like iidnorm_prior_sampler) inside the function being called
#   if absolutely needed or if there's a computationally less heavy/redundant
#   alternative implement that instead to get accessible sampler lists in the way
#   I've outlined above.
#
