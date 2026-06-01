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
