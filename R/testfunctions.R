#'make_test_functions
#'
#' creates a list of all test functions
#' @param ... a list of all used test functions.
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

#===============================================================================
#FOR GEWEKE ONLY (might delete later):
#see `tabulate_geweke_tests()` in geweke.R

#===============================================================================
#FOR SBC ONLY:

#'recompute_ranks
#'
#'For SBC only
#'
#'@export
recompute_ranks <- function(result, test_fns,
                            parallelize = FALSE,
                            n_thin = NULL) {

  n_sims <- length(result$sbc_result$fits)
  n_tests <- length(test_fns)

  if(parallelize == TRUE) {
    future::plan(future::multisession)
    on.exit(future::plan(future::sequential))
  }

  # thin if requested
  n_draws_total <- nrow(posterior::as_draws_matrix(result$sbc_result$fits[[1]]))
  thin_idx <- if (is.null(n_thin)) seq_len(n_draws_total)
              else seq(1, n_draws_total, by = n_thin)

  ranks_list <- future.apply::future_lapply(
    seq_len(n_sims),
    function(i) {
      theta_tilde <- as.numeric(result$dataset$variables[i, ])
      names(theta_tilde) <- colnames(result$dataset$variables)
      y <- result$dataset$generated[[i]]$y
      draws <- posterior::as_draws_matrix(result$sbc_result$fits[[i]])[thin_idx, ]
      row <- numeric(n_tests)
      names(row) <- names(test_fns)
      for (nm in names(test_fns)){
          fn <- test_fns[[nm]]
          tilde_val <- fn(theta_tilde, y)
          post_vals <- apply(draws, 1, function(theta) fn(theta, y))
          row[nm] <- sum(post_vals < tilde_val)
      }
      row
    },
    future.seed = TRUE
  )

  ranks <- do.call(rbind, ranks_list)
  colnames(ranks) <- names(test_fns)

  structure(
    list(ranks = ranks,
         n_draws = length(thin_idx),
         n_sims = n_sims),
    class = "bayescheckr_ranks"
  )
}

