#testfunctions

#quants <- derived_quantities(
#  product = mu * sigmasq,
#  log_like = sum(dnorm(y, mean = mu, sd = sqrt(sigmasq), log = TRUE)),
#)
#

#' Run Tests
#'
#' Test function from users input, ran all at once.
#' @param results the res object
#' @param L number of posterior draws
#' @param test
run_tests <- function(results, L, tests) {
  lapply(
    tests,
    function(f) f(results, L)
  )
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
