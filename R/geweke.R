#geweke.R

#inputs: the huge object outputted by simulate.R

# trying to generalize these samplers

# call the SBC pkg for all these functions and our pkg = downstream of theirs
#
# - SBC_generator_function
# - generate_datasets
#
# what I need to bring to the table: the joint posterior sampler
# Gibbs-specific alternating thing
# (must take same/compatible inputs as other functions SBC uses)
# only "successive-conditional" sampler; SBC has direct marginal-conditional
# decouple the main work from the second step = test functions
# test stuff: table with
# - the moments/summary stats
# - p-values for KS test and OG paper test (difference in means) z-score
# - convergence diagnostics (autocorrelation? R-squared? canned R functions)
# - whatever is inputted to make that table should also be a callable object
# - test functions = up to the user
# - two function prototype:
#   - geweke_summary(direct draws, gibbs draws, test functions):
#     - returns a table (aka dataframe)
#   - geweke_plot(direct draws, gibbs draws, test functions):
#     - doesn't return anything, just gives a panel for the QQ plots
# main subtlety/thing to watch out for:
#   user-friendly way to send to computer the test functions we want it to run
# end goal: Friday meeting will have him live-testing this as "blind user"
# format: long-formatted output matrix:
# -- e.g. es_linreg_saved$stats has these columns:
#         sim_id variable simulated_value rank z_score mean

# successive conditional draws using first principles as inputs

library(tidyverse)

#' the following two functions create the successive-conditional and
#' marginal-conditional sampling needed for comparison per John Geweke's 2004
#' JASA paper about "Getting it right."

geweke_suc_cond_draws <- function(prior_sampler,
                                  likelihood_sampler,
                                  posterior_sampler,
                                  prior_hyper_params,
                                  n_draws, n_obs) {

  #get info about parameters

  prior_draw <- prior_sampler(prior_hyper_params) #check accuracy of fct call

  n_params <- length(prior_draw)
  varnames <- names(prior_draw)

  #colnames (pre-rank/z-score computation): sim_id, variable, simulated_value
  #to be mutated() later: rank, z_score, mean

  n_row <- n_params * n_draws #number of params x number of simulations

  #preallocate storage for all of the draws
  draws_matrix <- data.frame(matrix(0, nrow = n_row, ncol = 3,
                            dimnames = list(1:n_row, #numbered rows
                                            c("sim_id", #col names
                                              "variable",
                                              "simulated_value")
                                            )
                            )) #sim_id, variable, simulated_value

  theta_matrix <- matrix(0, nrow = n_draws, ncol = n_params)
  y_matrix <- matrix(0, nrow = n_draws, ncol = n_obs)

  colnames(theta_matrix) <- varnames
  colnames(y_matrix) <- paste0("y", 1:n_obs)

  #initialize
  theta <- prior_draw
  y <- runif(n_obs) #using user input

  for (m in 1:n_draws) {
    #simulate from the posterior
    theta <- posterior_sampler(ndraws = 1,
                               y = y,
                               prior_hyper_params = prior_hyper_params)
    theta <- theta[1, ]
    theta_matrix[m, ] <- theta

    #simulate data given the latest parameter values
    y <- likelihood_sampler(n_obs, theta)
    y_matrix[m, ] <- y

    #cycle through rows for one simulation
    for (var in varnames) {
      index <- n_params*(m - 1) + which(varnames == var)
      draws_matrix[index, "sim_id"] <- m #constant across variables
      draws_matrix[index, "variable"] <- var
      draws_matrix[index, "simulated_value"] <- theta[var]
    }
  }

  #pivot to wide format

  draws_matrix_wide <- draws_matrix |>
    group_by(sim_id) |>
    dplyr::pivot_wider(names_from = variable, values_from = simulated_value)

  #return final object
  gibbs_draws <- list(draws = draws_matrix,
                      wide = draws_matrix_wide,
                      theta = theta_matrix,
                      y = y_matrix)

  return(gibbs_draws)
}

geweke_marg_cond_draws <- function(prior_sampler,
                                   likelihood_sampler,
                                   prior_hyper_params,
                                   n_params, n_draws, n_obs) {

  #set up the generator function for SBC

  generator_single <- function() {
    theta <- prior_sampler(prior_hyper_params)
    y <- likelihood_sampler(n_obs, theta)
  }

  #directly sample parameters and draws
  gen <- SBC::SBC_generator_function(generator_single)
  SBC::generate_datasets(gen, n_draws) # = n_sims? check

  #store
  theta_matrix <- as.matrix(data$variables)
  y_matrix <- do.call(rbind,
                      lapply(data$generated, function(g) g$y)) #from Claude
  colnames(y_matrix) <- paste0("y", 1:n_obs)

  direct_draws <- list(theta = theta_matrix,
       y = y_matrix,
       wide = cbind(theta_matrix, y_matrix))

  return(direct_draws)
}

#' Now moving to test functions: difference in means testing (per Geweke 2004)
#' and KS testing, for z-score comparison, among other tests

geweke_test <- function(g_mc, g_sc) {

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

all_the_tests <- function(direct_draws, gibbs_draws, test_functions) {

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

  #must return a matrix so I can save the following:

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
    #!need to add stats dependency

    #calculate: apply test function first, then pass resulting vectors
    ks <- stats::ks.test(g_mc, g_sc)

    #add results
    results_matrix["KS statistic", nm] <- ks$statistic
    results_matrix["KS p_value",   nm] <- ks$p.value

    #MCMC convergence test: last 2 rows --

    #!need to add coda dependency package

    #calculate: geweke.diag takes 1 mcmc chain object, using only Gibbs
    gibbs_chain <- coda::as.mcmc(g_sc)
    conv        <- coda::geweke.diag(gibbs_chain)

    #add results
    results_matrix["Convergence statistic", nm] <- as.numeric(conv$z)
    results_matrix["Convergence p_value",   nm] <- 2 * pnorm(-abs(conv$z))
  }

  return(results_matrix)
}

# NICE TO HAVES / CUT FOR TIME
# par(mfrow = c(3, 3))
# probs <- seq(0.05, 0.95, by = 0.05)
# pvals <- numeric(n + 2)
# labels <- c(paste0("beta", 1:p), "sigma", paste0("y", 1:n))
# for (i in 1:(n + p + 1)) {
#   quantiles_direct <- quantile(direct_draws[, i], probs)
#   quantiles_gibbs <- quantile(gibbs_draws[, i], probs)
#   plot(quantiles_direct, quantiles_gibbs, col = "red", pch = 19, cex = 0.5,
#        main = labels[i], xlab = "Direct draw", ylab = "Gibbs draw")
#   abline(a = 0, b = 1)
#   pvals[i] <- ks.test(direct_draws[, i], gibbs_draws[, i])$p.value
# }
#
# #Rafael's test functions
#
# par(mfrow = c(1, 1))
# probs <- seq(0.05, 0.95, by = 0.05)
#
# g_direct <- direct_draws[, 1] * direct_draws[, 2] * direct_draws[, 3]
# g_gibbs <- gibbs_draws[, 1] * gibbs_draws[, 2] * gibbs_draws[, 3]
#
# quantiles_direct <- quantile(g_direct, probs)
# quantiles_gibbs <- quantile(g_gibbs, probs)
# plot(quantiles_direct, quantiles_gibbs, col = "red", pch = 19, cex = 0.5,
#      main = "beta_1 * beta_2 * sigma^2")
# abline(a = 0, b = 1)
#
# #Histogram of the product:
#
# hist(g_direct[abs(g_direct) < 50], breaks = "Scott")
# hist(g_gibbs[abs(g_gibbs) < 50], breaks = "Scott")
#
# #Let me also use different test functions. First, let's look at *second moments*
#
# par(mfrow = c(1, 3))
# probs <- seq(0.05, 0.95, by = 0.05)
# pvals_second <- numeric(p + 1)
# labels_second <- c(paste("beta", 1:p, sep = "_"), "sigma2")
#
# for(i in 1:(p + 1)){
#   # compute second moments by squaring
#   second_direct <- direct_draws[, i]^2
#   second_gibbs <- gibbs_draws[, i]^2
#
#   quantiles_direct <- quantile(second_direct, probs)
#   quantiles_gibbs <- quantile(second_gibbs, probs)
#
#   plot(quantiles_direct, quantiles_gibbs,
#        col = "blue", pch = 19, cex = 0.5,
#        main = paste(labels_second[i], "squared"),
#        xlab = "direct", ylab = "gibbs")
#   abline(a = 0, b = 1)
# }
#
# #Finally let's look at the test function $h(\theta, y) = |\beta_1 - \beta_2|$:
#
# par(mfrow = c(1, 1))
# probs <- seq(0.05, 0.95, by = 0.05)
# h_direct <- abs(direct_draws[, 1] - direct_draws[, 2])
# h_gibbs <- abs(gibbs_draws[, 1] - gibbs_draws[, 2])
# quantiles_direct <- quantile(h_direct, probs)
# quantiles_gibbs <- quantile(h_gibbs, probs)
# plot(quantiles_direct, quantiles_gibbs, col = "red", pch = 19, cex = 0.5,
#      main = "|beta_1 - beta_2|",
#      xlab = "direct", ylab = "gibbs")
# abline(a = 0, b = 1)
