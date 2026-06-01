#geweke.R

#inputs: the huge object outputted by simulate.R

# trying to generalize these samplers

# call the SBC pkg for all these functions and our pkg = downstream of theirs
#
# - SBC_generator_function
# - generate_datasets
#
# what I need to bring to the table: the Gibbs-specific alternating thing
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


direct_joint_sampler <- function(prior_sampler,
                                 likelihood_sampler,
                                 prior_hyperparams,
                                 ndraws, n_obs) {

  draws_ncol <- 1 #some number dependent on length of n_obs, length of
    # preallocate storage for all of the draws

  #initialize to 0
  direct_draws <- matrix(0, nrow = ndraws, ncol = draws_ncol)

  colnames(direct_draws) <- c(paste0("beta", 1:p), "sigma", paste0("y", 1:n))

  for (m in 1:ndraws) {
    # simulate exactly from the prior
    theta <- prior_sampler(bbar0, v0, a0, b0)
    # simulate exactly from the likelihood given the draw from the prior
    5
    y <- likelihood_sampler(x, n, theta)
    direct_draws[m, ] <- c(theta, y)
  }
  # return an iid sample from the joint distribution
  return(direct_draws)
}

# trying to generalize
gibbs_joint_sampler <- function(prior_sampler,
                                likelihood_sampler,
                                prior_hyperparams,
                                ndraws, n_obs) {
  # preallocate storage for all of the draws
  p <- ncol(x)
  ## column number: n for y, p for beta, + 1 for sigma
  gibbs_draws <- matrix(0, nrow = ndraws, ncol = n + p + 1)
  colnames(gibbs_draws) <- c(paste0("beta", 1:p), "sigma", paste0("y", 1:n))
  # initialize
  theta <- c(rep(0, p), 1) #default beta (p-length vector) and sigma
  y <- runif(n)
  for(m in 1:ndraws){
    # simulate from the posterior using my sampler from above
    # the initialization is very important here!
    theta <- posterior_sampler(1, x, y, bbar0, v0, a0, b0, init = theta)
    theta <- as.vector(theta)
    # simulate data given the latest parameter values
    y <- likelihood_sampler(x, n, theta)
    6
    gibbs_draws[m, ] <- c(theta, y)
  }
  return(gibbs_draws)
}

get_geweke_draws <- function(prior_sampler,
                             likelihood_sampler,
                             posterior_sampler,
                             n_sims        = 100,
                             n_obs         = 50,
                             n_draws       = 1000,
                             prior_hyper_params) {

  direct_draws <- direct_joint_sampler(prior_sampler = prior_sampler,
                                       likelihood_sampler = likelihood_sampler,
                                       n_obs = n_obs)

  gibbs_draws <- gibbs_joint_sampler(posterior_sampler = posterior_sampler,
                                     n_obs = n_obs)

}

#check the marginals

par(mfrow = c(3, 3))
probs <- seq(0.05, 0.95, by = 0.05)
pvals <- numeric(n + 2)
labels <- c(paste0("beta", 1:p), "sigma", paste0("y", 1:n))
for(i in 1:(n + p + 1)){
  quantiles_direct <- quantile(direct_draws[, i], probs)
  quantiles_gibbs <- quantile(gibbs_draws[, i], probs)
  plot(quantiles_direct, quantiles_gibbs, col = "red", pch = 19, cex = 0.5,
       main = labels[i], xlab = "Direct draw", ylab = "Gibbs draw")
  abline(a = 0, b = 1)
  pvals[i] <- ks.test(direct_draws[, i], gibbs_draws[, i])$p.value
}

#Rafael's test functions

par(mfrow = c(1, 1))
probs <- seq(0.05, 0.95, by = 0.05)

g_direct <- direct_draws[, 1] * direct_draws[, 2] * direct_draws[, 3]
g_gibbs <- gibbs_draws[, 1] * gibbs_draws[, 2] * gibbs_draws[, 3]

quantiles_direct <- quantile(g_direct, probs)
quantiles_gibbs <- quantile(g_gibbs, probs)
plot(quantiles_direct, quantiles_gibbs, col = "red", pch = 19, cex = 0.5,
     main = "beta_1 * beta_2 * sigma^2")
abline(a = 0, b = 1)

#Histogram of the product:

hist(g_direct[abs(g_direct) < 50], breaks = "Scott")
hist(g_gibbs[abs(g_gibbs) < 50], breaks = "Scott")


#Let me also use different test functions. First, let's look at *second moments*:


par(mfrow = c(1, 3))
probs <- seq(0.05, 0.95, by = 0.05)
pvals_second <- numeric(p + 1)
labels_second <- c(paste("beta", 1:p, sep = "_"), "sigma2")

for(i in 1:(p + 1)){
  # compute second moments by squaring
  second_direct <- direct_draws[, i]^2
  second_gibbs <- gibbs_draws[, i]^2

  quantiles_direct <- quantile(second_direct, probs)
  quantiles_gibbs <- quantile(second_gibbs, probs)

  plot(quantiles_direct, quantiles_gibbs,
       col = "blue", pch = 19, cex = 0.5,
       main = paste(labels_second[i], "squared"),
       xlab = "direct", ylab = "gibbs")
  abline(a = 0, b = 1)
}

#Finally let's look at the test function $h(\theta, y) = |\beta_1 - \beta_2|$:


par(mfrow = c(1, 1))
probs <- seq(0.05, 0.95, by = 0.05)
h_direct <- abs(direct_draws[, 1] - direct_draws[, 2])
h_gibbs <- abs(gibbs_draws[, 1] - gibbs_draws[, 2])
quantiles_direct <- quantile(h_direct, probs)
quantiles_gibbs <- quantile(h_gibbs, probs)
plot(quantiles_direct, quantiles_gibbs, col = "red", pch = 19, cex = 0.5,
     main = "|beta_1 - beta_2|",
     xlab = "direct", ylab = "gibbs")
abline(a = 0, b = 1)
