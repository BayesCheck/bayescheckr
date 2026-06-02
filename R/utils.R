# Utils

#' prior_sampler
#'
#' Draws 1 from the prior. relies on having no other inputs than the function and drawing from it once
#' @param prior_input the function that details the actual distribution function
#' @param prior_hyperparams the vector of hyperparameters
prior_sampler <- function(prior_input, prior_hyperparams){
  prior_sample <- prior_input(prior_hyperparams)
  return(prior_sample)
}


#' likelihood_sampler
#'
#'Draws N times from the likelihood using the draw from the prior
#'NOTE: format the output later
#' @param n_obs number of observations
#' @param prior_samp the value drawn from the prior
#' @param likelihood_input function that tells the code what distribution to draw from and how many times
likelihood_sampler <- function(n_obs, prior_samp, likelihood_input) {
  n_obs <- as.integer(n_obs)
  likelihood_draws <- as.vector(likelihood_input(n_obs, prior_samp))
  return(likelihood_draws)
}


#' posterior_sampler
#'
#' The posterior given to us to validate, as opposed to the one we approximate with gibbs.
#' NOTE: format the output later
#' @param y the vector of likelihood draws
#' @param n_obs number of observations
#' @param posterior_input_function whatever the user thinks the posterior draws are supposed to be.
#' @param prior_hyperparams the vector of hyperparameters
posterior_sampler <- function(y, n_obs, posterior_input_function, prior_hyperparams){
  n_obs <- as.integer(n_obs)
  posterior_matrix <- posterior_input_function(y, n_obs, prior_hyperparams)
  return(posterior_matrix)
}



#' Rank Mean Test
#'
#' Test function for mean
#' @param results the res object
#' @param L number of posterior draws
rank_mean_test <- function(results, L) {

  stopifnot(inherits(results, "SBC_results"))

  ranks <- results$stats$rank

  observed_mean <- mean(ranks)
  expected_mean <- L / 2

  list(
    observed_mean = observed_mean,
    expected_mean = expected_mean,
    difference = observed_mean - expected_mean
  )
}

#' Rank Variance Test
#'
#' Test function for varience
#' @param results the res object
#' @param L number of posterior draws
rank_variance_test <- function(results, L) {

  stopifnot(inherits(results, "SBC_results"))

  ranks <- results$stats$rank

  observed_variance <- var(ranks)
  expected_variance <- L * (L + 2) / 12

  list(
    observed_variance = observed_variance,
    expected_variance = expected_variance,
    ratio = observed_variance / expected_variance
  )
}

#' Rank KS Test
#'
#' Test function for KS test
#' @param results the res object
#' @param L number of posterior draws
rank_ks_test <- function(results, L) {

  stopifnot(inherits(results, "SBC_results"))

  ranks <- results$stats$rank

  ks.test(
    ranks / L,
    "punif"
  )
}


#' Rank KL Divergence
#'
#' Test function for KL Divergence
#' @param results the res object
#' @param L number of posterior draws
rank_kl_divergence <- function(results, L) {

  stopifnot(inherits(results, "SBC_results"))

  ranks <- results$stats$rank

  counts <- tabulate(ranks + 1, nbins = L + 1)

  p_obs <- counts / sum(counts)
  p_exp <- rep(1 / (L + 1), L + 1)

  # Remove zero bins to avoid 0*log(0)
  nonzero <- p_obs > 0

  kl <- sum(
    p_obs[nonzero] *
      log(p_obs[nonzero] / p_exp[nonzero])
  )

  list(
    KL_divergence = kl,
    observed_probs = p_obs,
    expected_probs = p_exp
  )
}


#' compute_ranks
#'
#' Raphael pls explain
#' @param sims Raphael pls explain
#compute_ranks <- function(sims) {
#  param_names <- names(sims$sims[[1]]$theta_tilde)
#  n_params    <- length(param_names)

#  ranks <- matrix(NA,
#                  nrow     = sims$n_sims,
#                  ncol     = n_params,
#                 dimnames = list(NULL, param_names))
#
#  for (i in seq_len(sims$n_sims)) {
#    theta_tilde <- sims$sims[[i]]$theta_tilde
#    draws       <- sims$sims[[i]]$draws
#    for (p in param_names) {
#      ranks[i, p] <- sum(draws[, p] < theta_tilde[p])
#    }
#  }
#
#  structure(
#    list(ranks = ranks, n_draws = sims$n_draws, n_sims = sims$n_sims),
#    class = "bayescheckr_ranks"
#  )
#}


