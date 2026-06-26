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

.make_generator_single_testfns <- function(prior_sampler, likelihood_sampler,
                                           n_obs, prior_hyper_params, test_fns) {
  function() {
    theta <- prior_sampler(prior_hyper_params)
    y     <- likelihood_sampler(n_obs, theta)
    # variables: apply each test function to the prior draw
    # this is what SBC will rank the posterior values against
    vars <- lapply(test_fns, function(f) f(theta, y))
    list(
      variables = vars,
      generated = list(y = y, n = n_obs)
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
                    n_draws = 1e4,
                    n_burnin = 0,
                    n_thin = 10,
                    prior_hyper_params,
                    test_fns = NULL,          # <-- NEW argument
                    parallelize = FALSE,
                    globals = NULL) {

  theta_test <- prior_sampler(prior_hyper_params)

  # --- CHANGED: validation branches on whether test_fns is provided ---
  if (is.null(test_fns)) {
    # original path: theta must be a named numeric vector
    if (!is.numeric(theta_test))
      stop("prior_sampler must return a numeric vector.")
    if (is.null(names(theta_test)) || any(names(theta_test) == ""))
      stop("prior_sampler must return a named numeric vector.")
    param_names <- names(theta_test)
  }
  # if test_fns is provided, we don't validate theta structure at all —
  # it can be a list, vector, anything. test_fns defines the comparison.

  y_test <- likelihood_sampler(n_obs, theta_test)
  if (!is.numeric(y_test) && !is.matrix(y_test) && !is.list(y_test))
    stop("likelihood_sampler must return a numeric vector, matrix, or list.")

  # --- CHANGED: generator handles both cases ---
  if (is.null(test_fns)) {
    generator <- SBC::SBC_generator_function(
      .make_generator_single(prior_sampler, likelihood_sampler, n_obs,
                             prior_hyper_params)
    )
  } else {
    generator <- SBC::SBC_generator_function(
      .make_generator_single_testfns(prior_sampler, likelihood_sampler, n_obs,
                                     prior_hyper_params, test_fns)
    )
  }

  dataset <- SBC::generate_datasets(generator, n_sims)

  # --- CHANGED: backend branches on test_fns ---
  if (is.null(test_fns)) {
    # original path
    backend <- SBC::SBC_backend_function(
      func = function(generated) {
        res_raw <- posterior_sampler(
          ndraws             = n_draws,
          y                  = generated$y,
          prior_hyper_params = prior_hyper_params
        )
        thin_idx <- seq(n_burnin + 1, n_draws, by = n_thin)
        res_raw  <- res_raw[thin_idx, ]
        colnames(res_raw) <- param_names
        posterior::as_draws_matrix(res_raw)
      }
    )
  } else {
    # new path: posterior returns list of theta-lists, apply test_fns
    backend <- SBC::SBC_backend_function(
      func = function(generated) {
        draws_list <- posterior_sampler(
          ndraws             = n_draws,
          y                  = generated$y,
          prior_hyper_params = prior_hyper_params
        )
        # draws_list: list of length n_draws, each element is a theta-list
        result <- t(sapply(draws_list, function(th)
          sapply(test_fns, function(f) f(th, generated$y))))
        colnames(result) <- names(test_fns)
        posterior::as_draws_matrix(result)
      }
    )
  }

  if (parallelize == TRUE) {
    future::plan(future::multisession)
    on.exit(future::plan(future::sequential))

    # merge package-internal globals with user-supplied ones:

    internal_globals <- if (is.null(test_fns)) {
      c("posterior_sampler", "n_draws", "prior_hyper_params", "param_names", "n_burnin", "n_thin")
    } else {
      c("posterior_sampler", "n_draws", "prior_hyper_params", "test_fns", "n_burnin", "n_thin")
    }

    all_globals <- union(internal_globals, globals)

    res <- SBC::compute_SBC(dataset, backend, globals = all_globals)
  } else {
    res <- SBC::compute_SBC(dataset, backend)
  }

  structure(
    list(sbc_result = res,
         dataset    = dataset),
    class = "bayescheckr_sbc"
  )
}

# Convert a bayescheckr_ranks object to an SBC_results object for plotting
#' @export
ranks_to_sbc_results <- function(ranked) {

  # build the minimal stats dataframe the SBC plotting functions need
  n_sims   <- nrow(ranked$ranks)
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
