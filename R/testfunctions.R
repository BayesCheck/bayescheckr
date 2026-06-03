#' We require that the user gives test functions inside a NAMED list.
#'
#' ```r
#' test_fns <- make_test_functions(
#'   mu_sq = function(theta, y) theta["mu"]^2,
#'   product = function(theta, y) theta["mu"] * theta["sigmasq"],
#'   log_like = function(theta, y)
#'     sum(dnorm(y, theta["mu"],
#'               sqrt(theta["sigmasq"]),
#'               log = TRUE))
#' )
#' ```

#FOR BOTH GEWEKE AND SBC:

make_test_functions <- function(...) {
  fns <- list(...)

  # must be named
  if (is.null(names(fns)) || any(names(fns) == ""))
    stop("All test functions must be named, e.g. make_test_functions(mu_sq = function(theta, y) ...)")

  # each element must be a function
  if (!all(sapply(fns, is.function)))
    stop("All elements must be functions.")

  # each function must take exactly theta and y as arguments
  for (nm in names(fns)) {
    args <- names(formals(fns[[nm]]))
    if (!identical(args, c("theta", "y")))
      stop("Test function '", nm, "' must have exactly two arguments: theta and y.")
  }

  structure(fns, class = "bayescheckr_tests")
}

#' The user calls it like:
#' ```r
#' test_fns <- make_test_functions(
#'  mu_sq    = function(theta, y) theta["mu"]^2,
#'  product  = function(theta, y) theta["mu"] * theta["sigmasq"],
#'  log_like = function(theta, y) sum(dnorm(y, theta["mu"],
#'                                    sqrt(theta["sigmasq"]), log = TRUE))
#')
#' ```r

#===============================================================================
#FOR GEWEKE ONLY (might delete later):
#' @export
all_geweke_tests <- function(direct_draws, gibbs_draws, test_functions) {

  #allocate storage for results: 2 stats x n tests - same as testfunctions.R
  stats <- c("statistic", "p_value")
  tests <- c("Geweke", "KS", "Convergence")

  results_matrix <- matrix(NA, nrow = 6, ncol = length(test_functions),
                           dimnames =
                             list(
                               #rownames
                               c(paste("Geweke", stats),
                                 paste("KS", stats),
                                 paste("Convergence", stats)),
                               #colnames
                               names(test_functions)
                             ))

  #intentionally copying naming conventions from recompute_ranks()
  #to be consistent

  for (nm in names(test_functions)) {
    fn <- test_functions[[nm]]

    # apply test function to get scalar per draw
    g_mc <- apply(direct_draws$theta, 1, fn)
    g_sc <- apply(gibbs_draws$theta,  1, fn)

    #Geweke test: first 2 rows --

    #calculate: using test function written above, functions already applied
    gw <- geweke_test(g_mc, g_sc)

    #add results
    results_matrix["Geweke statistic", nm] <- gw$stat
    results_matrix["Geweke p_value",   nm] <- gw$p_value

    #KS test: next 2 rows --

    #calculate: apply test function first, then pass resulting vectors
    ks <- stats::ks.test(g_mc, g_sc)

    #add results
    results_matrix["KS statistic", nm] <- ks$statistic
    results_matrix["KS p_value",   nm] <- ks$p.value

    #MCMC convergence test: last 2 rows --

    #calculate: geweke.diag takes 1 mcmc chain object, using only Gibbs
    gibbs_chain <- coda::as.mcmc(g_sc)
    conv        <- coda::geweke.diag(gibbs_chain)

    #add results
    results_matrix["Convergence statistic", nm] <- as.numeric(conv$z)
    results_matrix["Convergence p_value",   nm] <- 2 * pnorm(-abs(conv$z))
  }

  return(data.frame(results_matrix))
}

#===============================================================================
#FOR SBC ONLY:

#' Test functions for SBC. We recompute ranks manually because the original
#'implementation in the SBC package was slow.
#'@param result is created as `result <- run_sbc(...)`
recompute_ranks <- function(result, test_fns) {

  n_sims <- length(result$sbc_result$fits)
  n_tests <- length(test_fns)

  #preallocate ranks matrix: n_sims x n_tests
  ranks <- matrix(NA,
                  nrow = n_sims,
                  ncol = n_tests,
                  dimnames = list(NULL, names(test_fns)))

  for (i in seq_len(n_sims)) {

    #extract the 3 things we need for simulation i
    theta_tilde <- as.numeric(result$dataset$variables[i, ])
    names(theta_tilde) <- colnames(result$dataset$variables)
    y <- result$dataset$generated[[i]]$y
    draws <- posterior::as_draws_matrix(result$sbc_result$fits[[i]])

    for (nm in names(test_fns)) {
      fn <- test_fns[[nm]]

      #compute f(theta_tilde, y): scalar
      tilde_val <- fn(theta_tilde, y)

      #compute f(theta_1, y) for each posterior draw: vector
      post_vals <- apply(draws, 1, function(theta) fn(theta, y))

      #rank of theta tilde among posterior draws
      ranks[i, nm] <- sum(post_vals < tilde_val)

    }
  }

  structure(
    list(ranks = ranks,
         n_draws = nrow(result$sbc_result$fits[[1]]),
         n_sims = n_sims),
    class = "bayescheckr_ranks"
  )
}

