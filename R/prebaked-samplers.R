#Prebaked Samplers for known models

#final step (delete when done): make sure to export everything

#iid Normal:

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

#simple linear regression model:

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
