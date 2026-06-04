#' @import SBC
NULL

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

#' @method SBC_fit bayescheckr_backend
#' @export
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

# The main entry point.
#' @export
run_sbc <- function(prior_sampler,
                    likelihood_sampler,
                    posterior_sampler,
                    n_sims,
                    n_obs,
                    n_draws,
                    prior_hyper_params,
                    parallelize = FALSE){

  # register S3 method at call time, not load time
  registerS3method("SBC_fit", "bayescheckr_backend",
                   SBC_fit.bayescheckr_backend,
                   envir = asNamespace("SBC"))

  #0. Validate stuff:
  theta_test <- prior_sampler(prior_hyper_params)
  if (!is.numeric(theta_test))
    stop("prior_sampler must return a numeric vector.")
  if (is.null(names(theta_test)) || any(names(theta_test) == ""))
    stop("prior_sampler must return a named numeric vector. ",
         "e.g. return(c(mu = mu, sigmasq = sigma2))")

  #1. Get param names from one trial draw of the prior
  param_names <- names(theta_test)

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
    globals <- c("SBC_fit.bayescheckr_backend", ".make_backend",
                 ".make_generator_single")
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

#' Convert a bayescheckr_ranks object to an SBC_results object for plotting
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

