#testfunctions

#quants <- derived_quantities(
#  product = mu * sigmasq,
#  log_like = sum(dnorm(y, mean = mu, sd = sqrt(sigmasq), log = TRUE)),
#)
#

#' We require that the user gives test functions inside a NAMED list.
#' Something like:
#' test_functions <- list(
#'  mu_sq = function(theta, y) theta["mu"]^2,
#'  product = function(theta, y) theta["mu"] * theta["sigmasq"],
#'  log_like = function(theta, y) sum(dnorm(y, theta["mu"],
#'                                sqrt(theta["sigmasq"]), log = TRUE))
#')

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
#' test_fns <- make_test_functions(
#'  mu_sq    = function(theta, y) theta["mu"]^2,
#'  product  = function(theta, y) theta["mu"] * theta["sigmasq"],
#'  log_like = function(theta, y) sum(dnorm(y, theta["mu"],
#'                                    sqrt(theta["sigmasq"]), log = TRUE))
#')

run_tests <- function(test_fns, theta_matrix, y_matrix) {

  #theta matrix: n_draws x n_params (either direct or MCMC draws)
  #y_matrix    : n_draws x n_obs

  results <- lapply(names(test_fns), function(nm) {
    fn <- test_fns[[nm]]

    #apply test function to each row (each draw):
    values <- sapply(seq_len(nrow(theta_matrix)), function(i) {
      fn(theta = theta_matrix[i, ], y = y_matrix[i, ])
    })

    data.frame(
      test_function = nm,
      values = values

    )

    })
}



