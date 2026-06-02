library(bayescheckr)   # devtools::load_all(".") if working from source
library(ggplot2)

set.seed(42)

# -----------------------------------------------------------------------------
# 0. Hyperparameters
# -----------------------------------------------------------------------------
prior_hyper_params <- list(
  mu_0    = 0,     # prior mean
  sigma_0 = 1,     # prior SD
  sigma   = 2      # known likelihood SD (fixed)
)

n_obs   <- 30     # observations per simulation
n_sims  <- 300    # number of SBC simulations
n_draws <- 500    # posterior draws per simulation


# =============================================================================
# 1. Define the prior, likelihood, and posterior samplers
# =============================================================================

# --- Prior ---
# Returns a *named* vector so bayescheckr can label parameters.
my_prior_sampler <- function(hyper) {
  mu <- rnorm(1, mean = hyper$mu_0, sd = hyper$sigma_0)
  c(mu = mu)
}

# --- Likelihood ---
# Given theta (a named vector) and n, return a length-n vector of observations.
my_likelihood_sampler <- function(n, theta) {
  rnorm(n, mean = theta["mu"], sd = prior_hyper_params$sigma)
}

# --- Correct posterior sampler ---
# Returns an n_draws x n_params matrix; each row is one posterior draw.
correct_posterior_sampler <- function(ndraws, y, prior_hyper_params) {
  n       <- length(y)
  sigma   <- prior_hyper_params$sigma
  mu_0    <- prior_hyper_params$mu_0
  sigma_0 <- prior_hyper_params$sigma_0

  sigma_n_sq <- 1 / (n / sigma^2 + 1 / sigma_0^2)
  mu_n       <- sigma_n_sq * (sum(y) / sigma^2 + mu_0 / sigma_0^2)

  mu_draws <- rnorm(ndraws, mean = mu_n, sd = sqrt(sigma_n_sq))
  matrix(mu_draws, ncol = 1, dimnames = list(NULL, "mu"))
}

# --- Broken posterior sampler (ignores data — draws from prior only) ---
# A sampler that forgets to condition on y; SBC should catch this.
broken_posterior_sampler <- function(ndraws, y, prior_hyper_params) {
  mu_draws <- rnorm(ndraws,
                    mean = prior_hyper_params$mu_0,
                    sd   = prior_hyper_params$sigma_0)
  matrix(mu_draws, ncol = 1, dimnames = list(NULL, "mu"))
}


# =============================================================================
# 2. SBC with the CORRECT posterior sampler
# =============================================================================
cat("Running SBC with correct posterior sampler...\n")

sbc_correct <- run_sbc(
  prior_sampler      = my_prior_sampler,
  likelihood_sampler = my_likelihood_sampler,
  posterior_sampler  = correct_posterior_sampler,
  n_sims             = n_sims,
  n_obs              = n_obs,
  n_draws            = n_draws,
  prior_hyper_params = prior_hyper_params,
  parallelize        = FALSE
)

# --- Inspect with native SBC plots ---
p_correct <- plot_rank_hist(sbc_correct$sbc_result)
print(p_correct + ggtitle("Correct sampler — rank histogram (should be flat)"))

p_ecdf_correct <- plot_ecdf_diff(sbc_correct$sbc_result)
print(p_ecdf_correct + ggtitle("Correct sampler — ECDF difference (should hug 0)"))

# --- Rank-based diagnostic statistics ---
cat("\n--- Correct sampler diagnostics ---\n")

mean_test <- rank_mean_test(sbc_correct$sbc_result, L = n_draws)
cat(sprintf("Rank mean:  observed = %.2f,  expected = %.2f,  diff = %.2f\n",
            mean_test$observed_mean, mean_test$expected_mean, mean_test$difference))

var_test <- rank_variance_test(sbc_correct$sbc_result, L = n_draws)
cat(sprintf("Rank variance ratio (obs/exp): %.3f  (1.0 = perfect)\n",
            var_test$ratio))

ks_test <- rank_ks_test(sbc_correct$sbc_result, L = n_draws)
cat(sprintf("KS test p-value: %.4f  (large p = no evidence of non-uniformity)\n",
            ks_test$p.value))

kl_div <- rank_kl_divergence(sbc_correct$sbc_result, L = n_draws)
cat(sprintf("KL divergence: %.4f  (near 0 = empirical ≈ uniform)\n",
            kl_div$KL_divergence))


# =============================================================================
# 3. SBC with the BROKEN posterior sampler
# =============================================================================
cat("\nRunning SBC with broken posterior sampler...\n")

sbc_broken <- run_sbc(
  prior_sampler      = my_prior_sampler,
  likelihood_sampler = my_likelihood_sampler,
  posterior_sampler  = broken_posterior_sampler,
  n_sims             = n_sims,
  n_obs              = n_obs,
  n_draws            = n_draws,
  prior_hyper_params = prior_hyper_params,
  parallelize        = FALSE
)

p_broken <- plot_rank_hist(sbc_broken$sbc_result)
print(p_broken + ggtitle("Broken sampler — rank histogram (should show a pattern)"))

cat("\n--- Broken sampler diagnostics ---\n")

mean_test_b <- rank_mean_test(sbc_broken$sbc_result, L = n_draws)
cat(sprintf("Rank mean:  observed = %.2f,  expected = %.2f,  diff = %.2f\n",
            mean_test_b$observed_mean, mean_test_b$expected_mean,
            mean_test_b$difference))

ks_test_b <- rank_ks_test(sbc_broken$sbc_result, L = n_draws)
cat(sprintf("KS test p-value: %.6f  (small p = sampler is miscalibrated)\n",
            ks_test_b$p.value))

kl_div_b <- rank_kl_divergence(sbc_broken$sbc_result, L = n_draws)
cat(sprintf("KL divergence: %.4f\n", kl_div_b$KL_divergence))


# =============================================================================
# 4. Geweke test for joint sampler validity
# =============================================================================
# The Geweke test compares the marginal distribution of test functions
# under (a) direct draws (theta ~ prior, y ~ likelihood(theta)) versus
# (b) the stationary Gibbs / MCMC chain.
#
# For the Normal-Normal model we have a closed-form posterior, so we
# simulate the "Gibbs" chain by forward-sampling as well — the two
# distributions should be identical and the QQ plots should fall on y = x.
# -----------------------------------------------------------------------------

cat("\nRunning Geweke QQ-plot check...\n")

n_geweke <- 2000   # number of draws for each distribution

# Direct draws: sample (theta, y) jointly from the prior-predictive
direct_theta <- matrix(NA, nrow = n_geweke, ncol = 1,
                       dimnames = list(NULL, "mu"))
direct_y     <- matrix(NA, nrow = n_geweke, ncol = n_obs)

for (i in seq_len(n_geweke)) {
  theta_i           <- my_prior_sampler(prior_hyper_params)
  direct_theta[i, ] <- theta_i
  direct_y[i, ]     <- my_likelihood_sampler(n_obs, theta_i)
}

# "Gibbs" draws: same procedure here (both should match marginals)
# In a real Gibbs scenario you would run your MCMC chain here.
gibbs_theta <- matrix(NA, nrow = n_geweke, ncol = 1,
                      dimnames = list(NULL, "mu"))
gibbs_y     <- matrix(NA, nrow = n_geweke, ncol = n_obs)

for (i in seq_len(n_geweke)) {
  theta_i          <- my_prior_sampler(prior_hyper_params)
  gibbs_theta[i, ] <- theta_i
  gibbs_y[i, ]     <- my_likelihood_sampler(n_obs, theta_i)
}

direct_draws <- list(theta = direct_theta, y = direct_y)
gibbs_draws  <- list(theta = gibbs_theta,  y = gibbs_y)

# Test functions: scalar summaries of (theta, y) to compare
test_fns <- make_test_functions(
  mu_identity = function(theta, y) theta["mu"],
  mu_sq       = function(theta, y) theta["mu"]^2,
  y_mean      = function(theta, y) mean(y),
  log_like    = function(theta, y)
    sum(dnorm(y, mean = theta["mu"],
              sd   = prior_hyper_params$sigma, log = TRUE))
)

# QQ plots: points should lie on the y = x diagonal
geweke_plot(
  direct_draws   = direct_draws,
  gibbs_draws    = gibbs_draws,
  test_functions = test_fns
)
title(main = "Geweke QQ plots — correct sampler (points should lie on y = x)",
      outer = TRUE, line = -1)


# =============================================================================
# 5. Accessing raw results for custom analyses
# =============================================================================
# The object returned by run_sbc() exposes the full SBC machinery.

# Posterior draws for simulation 1  (n_draws x n_params matrix)
post_draws_sim1 <- sbc_correct$sbc_result$fits[[1]]
cat("\nFirst 5 posterior draws for mu (simulation 1):\n")
print(head(post_draws_sim1, 5))

# True parameter value used in simulation 1
theta_tilde_sim1 <- sbc_correct$dataset$variables[1, ]
cat(sprintf("\nTrue mu in simulation 1: %.4f\n", theta_tilde_sim1$mu))

# Observed data for simulation 1
y_sim1 <- sbc_correct$dataset$generated[[1]]$y
cat(sprintf("Data summary — mean: %.4f, sd: %.4f\n",
            mean(y_sim1), sd(y_sim1)))
