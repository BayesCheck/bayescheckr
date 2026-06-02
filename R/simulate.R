# simulate.R

run_simulations <- function(prior_sampler,
                            likelihood_sampler,
                            posterior_sampler,
                            n_sims        = 100,
                            n_obs         = 50,
                            n_draws       = 1000,
                            prior_hyper_params) {

  sims <- vector("list", n_sims)

  for (i in seq_len(n_sims)) {

    # 1. draw from joint: theta_tilde ~ prior, y ~ likelihood | theta_tilde
    theta_tilde <- prior_sampler(prior_hyper_params)          # named numeric vector
    y           <- likelihood_sampler(n_obs, theta_tilde)     # numeric vector of data

    # 2. draw from posterior | y  ->  [n_draws x n_params] matrix, named columns
    draws <- posterior_sampler(ndraws            = n_draws,
                               y                 = y,
                               prior_hyper_params = prior_hyper_params)

    # store everything both frameworks will need
    sims[[i]] <- list(
      theta_tilde = theta_tilde,   # named numeric vector, e.g. c(mu=..., sigmasq=...)
      y           = y,             # the generated data
      draws       = draws          # matrix, colnames = param names
    )
  }

  structure(
    list(
      sims               = sims,
      n_sims             = n_sims,
      n_obs              = n_obs,
      n_draws            = n_draws,
      prior_hyper_params = prior_hyper_params
    ),
    class = "bayescheckr_sims"
  )
}

