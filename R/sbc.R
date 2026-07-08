# Used by run_geweke (geweke_mc_draws): raw theta draws, no test functions.
.make_generator_single <- function(prior_sampler, likelihood_sampler, n_obs,
                                   prior_hyper_params) {
  function() {
    theta <- prior_sampler(prior_hyper_params)
    y     <- likelihood_sampler(n_obs, theta, prior_hyper_params)
    list(
      variables = as.list(theta),
      generated = list(y = y,
                       n = n_obs)
    )
  }
}

# Used by run_sbc: variables are test-function evaluations of the prior draw.
.make_generator_single_test_fns <- function(prior_sampler, likelihood_sampler,
                                            n_obs, prior_hyper_params, test_fns) {
  function() {
    theta <- prior_sampler(prior_hyper_params)
    y     <- likelihood_sampler(n_obs, theta, prior_hyper_params)
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
                    prior_hyper_params,
                    test_fns,
                    parallelize = FALSE,
                    globals = NULL) {

  # ---- validate test_fns ----
  if (missing(test_fns) || is.null(test_fns))
    stop("test_fns is required. Build it with make_test_functions(); ",
         "for marginal checks use identity functions, e.g. ",
         "make_test_functions(mu = function(theta, y) theta$mu).")

  # accept either a bayescheckr_tests object or a plain named list of
  # functions — in the latter case, route through make_test_functions()
  # so all validation lives in one place
  if (!inherits(test_fns, "bayescheckr_tests"))
    test_fns <- do.call(make_test_functions, as.list(test_fns))


  theta_test <- prior_sampler(prior_hyper_params)

  # ---- validate prior_sampler output: must be a named list ----
  theta_test <- prior_sampler(prior_hyper_params)
  if (!is.list(theta_test))
    stop("prior_sampler must return a named list, ",
         "e.g. list(mu = mu, sigmasq = sigma2).")
  if (is.null(names(theta_test)) || any(names(theta_test) == ""))
    stop("prior_sampler must return a *named* list: every element needs a name.")


  y_test <- likelihood_sampler(n_obs, theta_test, prior_hyper_params)
  if (!is.numeric(y_test) && !is.matrix(y_test) && !is.list(y_test))
    stop("likelihood_sampler must return a numeric vector, matrix, or list.")


  # ---- validate each test function returns a numeric scalar ----
  for (nm in names(test_fns)) {
    val <- test_fns[[nm]](theta_test, y_test)
    if (!is.numeric(val) || length(val) != 1L)
      stop("Test function '", nm, "' must return a single numeric value; ",
           "got ", class(val)[1], " of length ", length(val), ".")
  }

  # ---- generator ----
  generator <- SBC::SBC_generator_function(
    .make_generator_single_test_fns(prior_sampler, likelihood_sampler, n_obs,
                                    prior_hyper_params, test_fns)
  )

  dataset <- SBC::generate_datasets(generator, n_sims)

  # ---- backend ----
  backend <- SBC::SBC_backend_function(
    func = function(generated) {
      draws_raw <- posterior_sampler(
        n_draws            = n_draws,
        y                  = generated$y,
        prior_hyper_params = prior_hyper_params
      )

      # normalize posterior output to a list of named theta-lists:
      #  - matrix with column names -> convert row-by-row (fixed-dim models)
      #  - list of named theta-lists -> use as-is (e.g. trans-dimensional models)
      if (is.matrix(draws_raw)) {
        cn <- colnames(draws_raw)
        if (is.null(cn) || any(cn == ""))
          stop("posterior_sampler returned a matrix without full column names. ",
               "Set colnames(THETA) to the parameter names before returning.")
        draws_list <- lapply(seq_len(nrow(draws_raw)),
                             function(s) as.list(draws_raw[s, ]))
      } else if (is.list(draws_raw)) {
        draws_list <- draws_raw
      } else {
        stop("posterior_sampler must return either a matrix with named ",
             "columns or a list of named theta-lists.")
      }

      # shape-safe reduction: n_draws x n_tests, guaranteed even when
      # length(test_fns) == 1 (the old t(sapply(...)) silently produced
      # a 1 x n_draws matrix in that case)
      result <- matrix(NA_real_,
                       nrow = length(draws_list),
                       ncol = length(test_fns),
                       dimnames = list(NULL, names(test_fns)))
      for (j in seq_along(test_fns)) {
        f <- test_fns[[j]]
        result[, j] <- vapply(draws_list,
                              function(th) f(th, generated$y),
                              numeric(1))
      }
      posterior::as_draws_matrix(result)
    }
  )


  if (parallelize == TRUE) {
    future::plan(future::multisession)
    on.exit(future::plan(future::sequential))

    # merge package-internal globals with user-supplied ones:

    internal_globals <- c("posterior_sampler", "n_draws",
                          "prior_hyper_params", "test_fns")
    all_globals <- union(internal_globals, globals)

    res <- SBC::compute_SBC(dataset, backend, globals = all_globals)
  } else {
    res <- SBC::compute_SBC(dataset, backend)
  }

  structure(
    list(sbc_result = res,
         dataset    = dataset,
         test_fns = test_fns),
    class = "bayescheckr_sbc"
  )
}

#' @export
run_sbc_prebaked <- function(
    model, #options: "iidnorm", "linreg", "doucet", "primiceri"
    n_sims,
    n_obs,
    n_draws            = 1e4,
    prior_hyper_params = NULL,
    test_fns           = NULL,
    parallelize        = FALSE,
    globals            = NULL
) {
  registry <- .make_sampler_registry()

  if (!model %in% names(registry$prior))
    stop("Unknown model: '", model, "'. Available models: ",
         paste(names(registry$prior), collapse = ", "))

  prior_sampler      <- registry$prior[[model]]
  likelihood_sampler <- registry$likelihood[[model]]
  posterior_sampler  <- registry$posterior[[model]]

  if (is.null(prior_hyper_params))
    prior_hyper_params <- registry$default_hyper[[model]]

  run_sbc(
    prior_sampler      = prior_sampler,
    likelihood_sampler = likelihood_sampler,
    posterior_sampler  = posterior_sampler,
    n_sims             = n_sims,
    n_obs              = n_obs,
    n_draws            = n_draws,
    prior_hyper_params = prior_hyper_params,
    test_fns           = test_fns,
    parallelize        = parallelize,
    globals            = globals
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
