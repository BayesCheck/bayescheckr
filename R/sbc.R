# sbc.R

#'compute_ranks
#'Computes ranks for parameters (no test function yet).
#'@param sims Object created by run_simulations (a list of simulations, where
#'each 'entry' contains a draws matrix)
compute_ranks <- function(sims) {

  #1. Extract parameter names from the first simulation

  param_names <- names(sims$sims[[1]]$theta_tilde)
  n_params <- length(param_names)

  #2. Preallocate the ranks matrix: n_sims x n_params

  ranks <- matrix(NA,
                  nrow     = sims$n_sims,
                  ncol     = n_params,
                  dimnames = list(NULL, param_names))

  #3. Loop over simulations and compute rank for each parameter

  for (i in seq_len(sims$n_sims)) {

    theta_tilde <- sims$sims[[i]]$theta_tilde   # named vector, length n_params
    draws <- sims$sims[[i]]$draws               # matrix [n_draws x n_params]

    for (p in param_names) {
      ranks[i, p] <- sum(draws[, p] < theta_tilde[p])
    }
  }

  #4. Return a classed object

  structure(
    list(
      ranks = ranks,
      n_draws = sims$n_draws,
      n_sims = sims$n_sims
    ),
    class = "bayescheckr_ranks"
  )
}
