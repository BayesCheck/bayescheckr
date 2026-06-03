#' @export
make_backend <- function(posterior_sampler, n_draws,
                          prior_hyper_params, param_names) {
  structure(
    list(sampling_func = posterior_sampler,
         ndraws = n_draws,
         prior_hyper_params = prior_hyper_params,
         param_names = param_names),
    class = "bayescheckr_backend"
  )
}
#' @export
SBC_fit_bayescheckr_backend <- function(backend, generated, cores){
  res_raw <- backend$sampling_func(
    ndraws = backend$ndraws,
    y = generated$y,
    prior_hyper_params = backend$prior_hyper_params
  )
  colnames(res_raw) <- backend$param_names
  posterior::as_draws_matrix(res_raw)
}
#' @export
make_generator_single <- function(prior_sampler, likelihood_sampler, n_obs,
                                   prior_hyper_params) {
  function() {
    theta <- prior_sampler(prior_hyper_params)
    y <- likelihood_sampler(n_obs, theta)
    list(
      variables = as.list(theta),
      generated = list(y = y,
                       n = n_obs)
    )
  }
}

#5. The main entry point.
#' @export
run_sbc <- function(prior_sampler,
                    likelihood_sampler,
                    posterior_sampler,
                    n_sims,
                    n_obs,
                    n_draws,
                    prior_hyper_params,
                    parallelize = FALSE){

  #1. Get param names from one trial draw of the prior
  param_names <- names(prior_sampler(prior_hyper_params))

  #2. Build the generator and dataset
  generator <- SBC::SBC_generator_function(
    make_generator_single(prior_sampler, likelihood_sampler, n_obs,
                           prior_hyper_params)
  )
  dataset <- SBC::generate_datasets(generator, n_sims)

  #3. Build the backend
  backend <- make_backend(posterior_sampler, n_draws, prior_hyper_params,
                           param_names)

  #4. Run SBC, maybe use parallelization
  if (parallelize == TRUE) {
    future::plan(future::multisession)
    on.exit(future::plan(future::sequential))
    globals <- c("SBC_fit_bayescheckr_backend")
    res <- SBC::compute_SBC(dataset, backend, globals = globals)
  } else {
    res <- SBC::compute_SBC(dataset, backend)
  }

  structure(
    list(
      sbc_result = res,      # the SBC_results object
      dataset    = dataset   # the SBC_datasets object
    ),
    class = "bayescheckr_sbc"
  )
}

#' A full picture of what run_sbc() returns:
#' result <- run_sbc(...)
#'
#' result$sbc_result                    # native SBC_results, use with SBC:: functions
#' result$dataset$variables             # n_sims x n_params matrix of theta_tilde draws
#' result$dataset$generated[[i]]$y      # y vector for simulation i
#' result$sbc_result$fits[[i]]          # n_draws x n_params posterior draws for simulation i

#===============================================================================
#TEST FUNCTIONS STUFF
#===============================================================================
#' @export
ranks_to_sbc_results <- function(ranked) {

  # build the minimal stats dataframe the SBC plotting functions need
  n_sims   <- ranked$n_sims
  n_tests  <- ncol(ranked$ranks)

  stats <- data.frame(
    sim_id   = rep(seq_len(n_sims), each = n_tests),
    variable = rep(colnames(ranked$ranks), times = n_sims),
    rank     = as.vector(t(ranked$ranks)),
    max_rank = ranked$n_draws
  )

  # bare minimum SBC_results shell — just enough for plotting
  structure(
    list(
      stats               = stats,
      fits                = vector("list", n_sims),
      backend_diagnostics = NULL,
      default_diagnostics = data.frame(sim_id = seq_len(n_sims)),
      outputs             = vector("list", n_sims),
      messages            = vector("list", n_sims),
      warnings            = vector("list", n_sims),
      errors              = vector("list", n_sims)
    ),
    class = "SBC_results"
  )
}


