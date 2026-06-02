.make_backend <- function(posterior_sampler, n_draws,
                          prior_hyper_params, param_names) {
  structure(
    list(sampling_func = posterior_sampler,
         ndraws = n_draws,
         prior_hyper_params = prior_hyper_params,
         param_names = param_names),
    class = "bayescheckr_backend"
  )
}

SBC_fit.bayescheckr_backend <- function(backend, generated, cores){
  res_raw <- backend$sampling_func(
    ndraws = backend$ndraws,
    y = generated$y,
    prior_hyper_params = backend$prior_hyper_params
  )
  colnames(res_raw) <- backend$param_names
  posterior::as_draws_matrix(res_raw)
}

.make_generator_single <- function(prior_sampler, likelihood_sampler, n_obs,
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
    .make_generator_single(prior_sampler, likelihood_sampler, n_obs,
                           prior_hyper_params)
  )
  dataset <- SBC::generate_datasets(generator, n_sims)

  #3. Build the backend
  backend <- .make_backend(posterior_sampler, n_draws, prior_hyper_params,
                           param_names)

  #4. Run SBC, maybe use parallelization
  if (parallelize == TRUE) {
    future::plan(future::multisession)
    on.exit(future::plan(future::sequential))
    globals <- c("SBC_fit.bayescheckr_backend")
    res <- SBC::compute_SBC(dataset, backend, globals = globals)
  } else {
    res <- SBC::compute_SBC(dataset, backend)
  }

  res
}



