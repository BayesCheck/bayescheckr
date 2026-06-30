#geweke.R

#inputs: the huge object outputted by simulate.R

# trying to generalize these samplers

# call the SBC pkg for all these functions and our pkg = downstream of theirs
#
# - SBC_generator_function
# - generate_datasets

#' the following two functions create the successive-conditional and
#' marginal-conditional sampling needed for comparison per John Geweke's 2004
#' JASA paper about "Getting it right."

#First, we make sure theta is constructed in a way that makes it robust to different access methods

#' @export
as_theta <- function(x) {
  if (is.list(x)) {
    return(x) #lists already support $ and [[ — nothing to do
  }
  class(x) <- c("bcheck_theta", class(x))
  x
}

#' @export
`$.bcheck_theta` <- function(x, name) {
  unclass(x)[[name]]
}


#successive conditional: sc

#helper function to determine whether to keep a draw

keep_draw <- function(draw_num, n_burn, n_thin) {
  if (draw_num < n_burn) {
    return(FALSE)
  }

  if ((draw_num - n_burn) %% n_thin == 0) {
    return(TRUE)
  }

  return(FALSE)
}

#' @export
geweke_sc_draws <- function(prior_sampler,
                            likelihood_sampler,
                            posterior_sampler,
                            n_obs,
                            n_draws = 1e4,
                            n_burn = 1000,
                            n_thin = 5,
                            prior_hyper_params) {

  #set total loop iterations using burn-in and thinning inputs
  total_draws <- n_burn + (n_draws - 1) *n_thin

  saved_draws <- n_draws #floor(n_draws / n_thin)

  #get info about parameters

  prior_draw <- prior_sampler(prior_hyper_params) #check accuracy of fct call

  n_params <- length(prior_draw)
  varnames <- names(prior_draw)
  n_row <- n_params * saved_draws #num params x num saved simulations

  #preallocate storage for all of the draws
  draws_matrix <- data.frame(matrix(0, nrow = n_row, ncol = 3,
                            dimnames = list(1:n_row, #numbered rows
                                            c("sim_id", #col names
                                              "variable",
                                              "simulated_value")
                                            )
                            )) #sim_id, variable, simulated_value

  theta_matrix <- matrix(0, nrow = saved_draws, ncol = n_params)
  y_matrix <- matrix(0, nrow = saved_draws, ncol = n_obs)

  colnames(theta_matrix) <- varnames
  colnames(y_matrix) <- paste0("y", 1:n_obs)

  #initialize
  theta <- prior_draw
  y <- runif(n_obs) #using user input

  for (m in 1:total_draws) {
    #simulate from the posterior

    theta <- posterior_sampler(n_draws = 1,
                               y = y,
                               prior_hyper_params = prior_hyper_params,
                               init = theta)
    colnames(theta) <- varnames
    theta <- theta[1, ]

    #simulate data given the latest parameter values
    y <- likelihood_sampler(n_obs, as_theta(theta))

    if (keep_draw(m, n_burn, n_thin)) { #decide whether to save the draw

      #calculate simulation ID
      n <- 1 + (m - n_burn) / n_thin

      #store draws in respective matrices
      theta_matrix[n, ] <- theta
      y_matrix[n, ] <- y

      #for long matrix: cycle through rows for one simulation
      for (var in varnames) {
        index <- n_params*(n - 1) + which(varnames == var)

        draws_matrix[index, "sim_id"] <- n #constant across variables
        draws_matrix[index, "variable"] <- var
        draws_matrix[index, "simulated_value"] <- theta[var]
      }
    }
  }

  #pivot to wide format: 1 row/simulation

  draws_matrix_wide <- draws_matrix |>
    dplyr::group_by(sim_id) |>
    tidyr::pivot_wider(names_from = variable, values_from = simulated_value)

  #return final object
  gibbs_draws <- list(draws = draws_matrix,
                      wide = draws_matrix_wide,
                      theta = theta_matrix,
                      y = y_matrix)

  return(gibbs_draws)
}

#using SBC to speed up my joint sampler
#marginal conditional: mc
#' @export
geweke_mc_draws <- function(prior_sampler,
                            likelihood_sampler,
                            n_obs,
                            n_draws = 1e4,
                            prior_hyper_params,
                            parallelize = TRUE) {

  if (parallelize) {
    future::plan(future::multisession)
    on.exit(future::plan(future::sequential))
  }

  #set up the generator function for SBC

  gen <- SBC::SBC_generator_function(
    .make_generator_single(prior_sampler, likelihood_sampler,
                           n_obs, prior_hyper_params)
  )

  gen_data <- SBC::generate_datasets(gen, n_draws) # = n_sims? check

  #store

  #theta_matrix <- as.matrix(data$variables)

  theta_matrix <- as.matrix(as.data.frame(gen_data$variables))
  colnames(theta_matrix) <- names(prior_sampler(prior_hyper_params))

  y_matrix <- do.call(rbind,
                      lapply(gen_data$generated, function(g) g$y)) #from Claude
  colnames(y_matrix) <- paste0("y", 1:n_obs)

  wide_matrix <- cbind(theta_matrix, y_matrix)

  draws_matrix <- data.frame(gen_data$variables) |>
    dplyr::mutate(sim_id = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = -sim_id,
                        names_to = "variable",
                        values_to = "simulated_value")

  direct_draws <- list(draws = draws_matrix,
                       theta = theta_matrix,
                       y = y_matrix,
                       wide = wide_matrix)

  return(direct_draws)
}

#' @export
run_geweke <- function(prior_sampler,
                       likelihood_sampler,
                       posterior_sampler,
                       n_obs,
                       n_draws = 1e4,
                       n_burn = 1000,
                       n_thin = 5,
                       prior_hyper_params,
                       test_fns = NULL,
                       parallelize = TRUE) {

  #0. Validate stuff:
  theta_test <- prior_sampler(prior_hyper_params)

  if (is.null(test_fns)) {
    if (!is.numeric(theta_test))
      stop("prior_sampler must return a numeric vector.")
    if (is.null(names(theta_test)) || any(names(theta_test) == ""))
      stop("prior_sampler must return a named numeric vector. ",
           "e.g. return(c(mu = mu, sigmasq = sigma2))")
    }

  # new check for likelihood_sampler
  y_test <- likelihood_sampler(n_obs, theta_test)
  if (!is.numeric(y_test) && !is.matrix(y_test) && !is.list(y_test))
    stop("likelihood_sampler must return a numeric vector, matrix, or list.")
  if (is.vector(y_test) && !is.list(y_test)) {
    # coerce to matrix silently so downstream code always sees consistent format
    warning("likelihood_sampler returned a plain vector - consider returning a matrix for generality.")
  }

  direct_draws <- geweke_mc_draws(prior_sampler = prior_sampler,
                                  likelihood_sampler = likelihood_sampler,
                                  n_obs = n_obs,
                                  n_draws = n_draws,
                                  prior_hyper_params = prior_hyper_params,
                                  parallelize)

  gibbs_draws <- geweke_sc_draws(prior_sampler = prior_sampler,
                                 likelihood_sampler = likelihood_sampler,
                                 posterior_sampler = posterior_sampler,
                                 n_obs = n_obs,
                                 n_draws = n_draws,
                                 n_burn = n_burn,
                                 n_thin = n_thin,
                                 prior_hyper_params = prior_hyper_params)
  if (!is.null(test_fns)) {
    return(apply_test_functions(direct_draws = direct_draws,
                         gibbs_draws = gibbs_draws,
                         test_functions = test_fns))
  }
  return(list(direct = direct_draws, gibbs = gibbs_draws))
}

#' @export
run_geweke_prebaked <- function(
    model, #string options: "iidnorm", "linreg", "doucet", "primiceri"
    n_obs,
    n_draws            = 1e4,
    n_burn             = 1000,
    n_thin             = 5,
    prior_hyper_params = NULL,
    test_fns           = NULL,
    parallelize        = TRUE
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

  run_geweke(
    prior_sampler      = prior_sampler,
    likelihood_sampler = likelihood_sampler,
    posterior_sampler  = posterior_sampler,
    n_obs              = n_obs,
    n_draws            = n_draws,
    n_burn             = n_burn,
    n_thin             = n_thin,
    prior_hyper_params = prior_hyper_params,
    test_fns           = test_fns,
    parallelize        = parallelize
  )
}

#' Now moving to test functions: difference in means testing (per Geweke 2004)
#' and KS testing, for z-score comparison, among other tests
#' @export
geweke_test <- function(g_mc,
                        g_sc) { #input: vectors, length = n_draws

  #get terms for test statistics --

  #marginal-conditional
  m1 <- length(g_mc) #length of values; m1 = m2
  mean_mc <- mean(g_mc, na.rm = TRUE) #test function means = form stat numerator
  var_mc <- var(g_mc,  na.rm = TRUE)

  #successive-conditional
  m2 <- length(g_sc)
  mean_sc <- mean(g_sc, na.rm = TRUE)
  tausq_sc <- coda::spectrum0.ar(g_sc)$spec #long-run variance

  #calculate the test stat and p-value --

  test_stat <- (mean_mc - mean_sc) / sqrt(var_mc/m1 + tausq_sc/m2)
  p_value <- 2 * pnorm(-abs(test_stat))

  return(list(stat = test_stat, p_value = p_value))
}

#' @export
apply_test_functions <- function(direct_draws, gibbs_draws, test_functions) {

  n_row <- nrow(direct_draws$theta) #same for Gibbs

  direct_test_matrix <- matrix(NA, nrow = n_row,
                               ncol = length(test_functions),
                               dimnames = list(NULL, names(test_functions)))

  gibbs_test_matrix <- matrix(NA, nrow = n_row,
                              ncol = length(test_functions),
                              dimnames = list(NULL, names(test_functions)))

  for (nm in names(test_functions)) {
    fn <- test_functions[[nm]]

    direct_test_matrix[, nm] <- as.numeric(sapply(
      seq_len(n_row), function(i) {
        fn(as_theta(direct_draws$theta[i, ]), as.numeric(direct_draws$y[i, ]))
      }
    ))

    gibbs_test_matrix[, nm] <- as.numeric(sapply(
      seq_len(n_row), function(i) {
        fn(as_theta(gibbs_draws$theta[i, ]), as.numeric(gibbs_draws$y[i, ]))
      }
    ))
  }

  return(list(direct = direct_test_matrix, gibbs = gibbs_test_matrix))
}

#' @export
tabulate_geweke_tests <- function(draws,
                                  test_functions,
                                  tests_applied = FALSE) {

  #allocate storage for results: 2 stats x n tests - same as testfunctions.R
  stats <- c("statistic", "p_value")
  tests <- c("Geweke", "KS", "Convergence")

  stat_names <- c("Geweke statistic", "Geweke p_value",
                  "KS statistic", "KS p_value",
                  "Convergence statistic", "Convergence p_value")

  if (tests_applied) {
    test_matrix <- draws
  } else {
    test_matrix <- apply_test_functions(draws$direct,
                                        draws$gibbs,
                                        test_functions)
  }

  results_matrix <- matrix(NA, nrow = length(test_functions), ncol = 6,
                           dimnames =
                             list(
                               #colnames
                               names(test_functions),
                               #rownames
                               c(paste("Geweke", stats),
                                 paste("KS", stats),
                                 paste("Convergence", stats))
                             ))

  #intentionally copying naming conventions from recompute_ranks()
  #to be consistent

  for (nm in names(test_functions)) {

    # g_mc <- as.numeric(sapply(seq_len(nrow(direct_draws$theta)), function(i) {
    #   fn(direct_draws$theta[i, ], as.numeric(direct_draws$y[i, ]))
    # }))
    #
    # g_sc <- as.numeric(sapply(seq_len(nrow(gibbs_draws$theta)), function(i) {
    #   fn(gibbs_draws$theta[i, ], as.numeric(gibbs_draws$y[i, ]))
    # }))

    g_mc <- test_matrix$direct[, nm]
    g_sc <- test_matrix$gibbs[, nm]

    #Geweke test: first 2 rows --

    #calculate: using test function written above, functions already applied
    gw <- geweke_test(g_mc, g_sc)

    #add results
    results_matrix[nm, stat_names[1]] <- gw$stat
    results_matrix[nm, stat_names[2]] <- gw$p_value

    #KS test: next 2 rows --

    #calculate: apply test function first, then pass resulting vectors
    ks <- stats::ks.test(g_mc, g_sc)

    #add results
    results_matrix[nm, stat_names[3]] <- ks$statistic
    results_matrix[nm, stat_names[4]] <- ks$p.value

    #MCMC convergence test: last 2 rows --

    #calculate: geweke.diag takes 1 mcmc chain object, using only Gibbs
    gibbs_chain <- coda::as.mcmc(g_sc)
    conv        <- coda::geweke.diag(gibbs_chain)

    #add results
    results_matrix[nm, stat_names[5]] <- as.numeric(conv$z)
    results_matrix[nm, stat_names[6]] <- 2 * pnorm(-abs(conv$z))
  }

  return(as.data.frame(results_matrix, col.names = stat_names))
}

#' @export
plot_geweke_tests <- function(draws,
                              test_functions,
                              tests_applied = FALSE,
                              probs = seq(0.05, 0.95, by = 0.05)) {

  if (tests_applied) {
    test_matrix <- draws
  } else {
    test_matrix <- apply_test_functions(draws$direct,
                                        draws$gibbs,
                                        test_functions)
  }

  #probs <- seq(0.05, 0.95, by = 0.05)
  n_fns <- length(test_functions)

  # set up plot grid
  par(mfrow = c(ceiling(n_fns / 3), 3))

  for (nm in names(test_functions)) {

    g_mc <- test_matrix$direct[, nm]
    g_sc <- test_matrix$gibbs[, nm]

    # compute quantiles and plot
    quantiles_direct <- quantile(g_mc, probs)
    quantiles_gibbs  <- quantile(g_sc, probs)

    plot(quantiles_direct, quantiles_gibbs,
         col = "red", pch = 19, cex = 0.5,
         main = nm,
         xlab = "Direct draw",
         ylab = "Gibbs draw")
    abline(a = 0, b = 1)
  }

  # reset plot layout
  par(mfrow = c(1, 1))
}

