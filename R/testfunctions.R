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
#' for both geweke and sbc
#' @export
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
#see `tabulate_geweke_tests()` in geweke.R

#===============================================================================
#FOR SBC ONLY:

#' Test functions for SBC. We recompute ranks manually because the original
#'implementation in the SBC package was slow.
#'@param result is created as `result <- run_sbc(...)`
#' @export
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

