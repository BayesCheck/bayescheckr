# =============================================================================
# bayescheckr example: Normal-Normal conjugate model
# =============================================================================
#
# Model:
#   Prior:      mu ~ Normal(mu_0, sigma_0^2),   sigma^2 fixed and known
#   Likelihood: y_i | mu ~ Normal(mu, sigma^2),  i = 1, ..., n
#
# Closed-form posterior:
#   mu | y ~ Normal(mu_n, sigma_n^2)
#   where:
#     sigma_n^2 = 1 / (n / sigma^2 + 1 / sigma_0^2)
#     mu_n      = sigma_n^2 * (sum(y) / sigma^2 + mu_0 / sigma_0^2)
#
# Because the posterior is available in closed form, this is an ideal
# validation target: any miscalibration the package detects is unambiguously
# a sampler bug, not a model approximation.
#
# Sections:
#   1. Define samplers
#   2. SBC — correct sampler  (ranks should be uniform)
#   3. SBC — broken sampler   (ranks should be non-uniform)
#   4. Geweke — correct Gibbs (QQ plots should hug y = x; all_the_tests passes)
#   5. Geweke — broken Gibbs  (QQ plots off diagonal; all_the_tests fails)
# =============================================================================

set.seed(42)

# -----------------------------------------------------------------------------
# 0. Shared configuration
# -----------------------------------------------------------------------------
prior_hyper_params <- list(
  mu_0    = 0,   # prior mean
  sigma_0 = 1,   # prior SD
  sigma   = 2    # known likelihood SD (fixed)
)

n_obs    <- 30    # observations per simulation
n_sims   <- 300   # SBC simulations
n_draws  <- 500   # posterior draws per SBC simulation
n_geweke <- 2000  # draws for each Geweke arm


# =============================================================================
# 1. Samplers
# =============================================================================

# --- Prior ---
# Must return a *named* vector; names become the parameter labels everywhere.
my_prior_sampler <- function(hyper) {
  c(mu = rnorm(1, mean = hyper$mu_0, sd = hyper$sigma_0))
}

# --- Likelihood ---
# Signature: (n_obs, theta) -> length-n numeric vector.
my_likelihood_sampler <- function(n, theta) {
  rnorm(n, mean = theta["mu"], sd = prior_hyper_params$sigma)
}

# --- Correct posterior sampler (conjugate update) ---
# Signature required by bayescheckr: (ndraws, y, prior_hyper_params)
#   -> ndraws x n_params matrix, column named "mu".
correct_posterior_sampler <- function(ndraws, y, prior_hyper_params) {
  n       <- length(y)
  sigma   <- prior_hyper_params$sigma
  mu_0    <- prior_hyper_params$mu_0
  sigma_0 <- prior_hyper_params$sigma_0

  sigma_n_sq <- 1 / (n / sigma^2 + 1 / sigma_0^2)
  mu_n       <- sigma_n_sq * (sum(y) / sigma^2 + mu_0 / sigma_0^2)

  matrix(rnorm(ndraws, mean = mu_n, sd = sqrt(sigma_n_sq)),
         ncol = 1, dimnames = list(NULL, "mu"))
}

# --- Broken posterior sampler ---
# Ignores the data and draws from the prior — SBC and Geweke should both
# flag this clearly.
broken_posterior_sampler <- function(ndraws, y, prior_hyper_params) {
  matrix(rnorm(ndraws,
               mean = prior_hyper_params$mu_0,
               sd   = prior_hyper_params$sigma_0),
         ncol = 1, dimnames = list(NULL, "mu"))
}

# --- Test functions ---
# Passed to both run_geweke_tests / all_the_tests (Geweke) and
# recompute_ranks (SBC). Each must have exactly (theta, y) as arguments.
test_fns <- make_test_functions(
  mu          = function(theta, y) theta["mu"],
  mu_sq       = function(theta, y) theta["mu"]^2,
  y_mean      = function(theta, y) mean(y),
  log_like    = function(theta, y)
    sum(dnorm(y, mean = theta["mu"],
              sd = prior_hyper_params$sigma, log = TRUE))
)


# =============================================================================
# 2. SBC — correct posterior sampler
# =============================================================================
cat("=== SBC: correct sampler ===\n")

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

# 2a. Native SBC plots on raw parameter ranks --------------------------------
p1 <- SBC::plot_rank_hist(sbc_correct$sbc_result)
print(p1 + ggtitle("Correct sampler — rank histogram (should be flat)"))

p2 <- SBC::plot_ecdf_diff(sbc_correct$sbc_result)
print(p2 + ggtitle("Correct sampler — ECDF difference (should hug 0)"))

# 2b. Recompute ranks on *test functions*, then plot -------------------------
# recompute_ranks() lets you check SBC calibration for derived quantities,
# not just raw parameters.
ranked_correct <- recompute_ranks(sbc_correct, test_fns)
sbc_shell_correct <- .ranks_to_sbc_results(ranked_correct)

p3 <- SBC::plot_rank_hist(sbc_shell_correct)
print(p3 + ggtitle("Correct sampler — test-function rank histograms"))

p4 <- SBC::plot_ecdf_diff(sbc_shell_correct)
print(p4 + ggtitle("Correct sampler — test-function ECDF diff"))

# 2c. Scalar rank diagnostics ------------------------------------------------
cat("\n--- Correct sampler: rank diagnostics (raw mu) ---\n")

m <- rank_mean_test(sbc_correct$sbc_result, L = n_draws)
cat(sprintf("  Mean:     observed = %.2f  |  expected = %.2f  |  diff = %.2f\n",
            m$observed_mean, m$expected_mean, m$difference))

v <- rank_variance_test(sbc_correct$sbc_result, L = n_draws)
cat(sprintf("  Variance: obs/exp ratio = %.3f  (1.0 = perfect)\n", v$ratio))

k <- rank_ks_test(sbc_correct$sbc_result, L = n_draws)
cat(sprintf("  KS p-value: %.4f  (large = no evidence of non-uniformity)\n",
            k$p.value))

kl <- rank_kl_divergence(sbc_correct$sbc_result, L = n_draws)
cat(sprintf("  KL divergence: %.4f  (near 0 = empirical ≈ uniform)\n",
            kl$KL_divergence))


# =============================================================================
# 3. SBC — broken posterior sampler
# =============================================================================
cat("\n=== SBC: broken sampler ===\n")

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

p5 <- SBC::plot_rank_hist(sbc_broken$sbc_result)
print(p5 + ggtitle("Broken sampler — rank histogram (should show non-uniformity)"))

ranked_broken <- recompute_ranks(sbc_broken, test_fns)
sbc_shell_broken <- .ranks_to_sbc_results(ranked_broken)

p6 <- SBC::plot_rank_hist(sbc_shell_broken)
print(p6 + ggtitle("Broken sampler — test-function rank histograms"))

cat("\n--- Broken sampler: rank diagnostics (raw mu) ---\n")

m_b <- rank_mean_test(sbc_broken$sbc_result, L = n_draws)
cat(sprintf("  Mean:     observed = %.2f  |  expected = %.2f  |  diff = %.2f\n",
            m_b$observed_mean, m_b$expected_mean, m_b$difference))

k_b <- rank_ks_test(sbc_broken$sbc_result, L = n_draws)
cat(sprintf("  KS p-value: %.6f  (small = miscalibrated)\n", k_b$p.value))

kl_b <- rank_kl_divergence(sbc_broken$sbc_result, L = n_draws)
cat(sprintf("  KL divergence: %.4f\n", kl_b$KL_divergence))


# =============================================================================
# 4. Geweke — correct Gibbs chain
# =============================================================================
# geweke_marg_cond_draws() samples (theta, y) jointly from the joint
# prior-predictive — these are the "direct" (marginal-conditional) draws.
#
# geweke_suc_cond_draws() runs the Gibbs/MCMC chain: alternates between
# drawing theta | y from the posterior and y | theta from the likelihood.
# Under a correct sampler the two arms share the same stationary distribution,
# so all test-function QQ plots should fall on y = x and all_the_tests()
# should return large p-values.
# =============================================================================
cat("\n=== Geweke: correct Gibbs ===\n")

direct_correct <- geweke_marg_cond_draws(
  prior_sampler      = my_prior_sampler,
  likelihood_sampler = my_likelihood_sampler,
  prior_hyper_params = prior_hyper_params,
  n_params           = 1,
  n_draws            = n_geweke,
  n_obs              = n_obs
)

gibbs_correct <- geweke_suc_cond_draws(
  prior_sampler      = my_prior_sampler,
  likelihood_sampler = my_likelihood_sampler,
  posterior_sampler  = correct_posterior_sampler,
  prior_hyper_params = prior_hyper_params,
  n_draws            = n_geweke,
  n_obs              = n_obs
)

# 4a. QQ plots — points should lie on the y = x diagonal --------------------
geweke_plot(
  direct_draws   = direct_correct,
  gibbs_draws    = gibbs_correct,
  test_functions = test_fns
)
title("Geweke QQ — correct sampler (should hug y = x)", outer = TRUE, line = -1)

# 4b. run_geweke_tests(): scalar values per draw for each test function ------
# Returns a long data.frame: columns test_function and values.
# Useful if you want to build your own summaries or plots downstream.
g_direct_vals <- run_geweke_tests(test_fns,
                                  theta_matrix = direct_correct$theta,
                                  y_matrix     = direct_correct$y)

g_gibbs_vals  <- run_geweke_tests(test_fns,
                                  theta_matrix = gibbs_correct$theta,
                                  y_matrix     = gibbs_correct$y)

cat("run_geweke_tests() output (first 6 rows, direct arm):\n")
print(head(g_direct_vals))

# 4c. all_the_tests(): Geweke z-test + KS test + MCMC convergence -----------
# Returns a data.frame with rows:
#   Geweke statistic / p_value
#   KS statistic     / p_value
#   Convergence statistic / p_value   (coda::geweke.diag on the Gibbs chain)
# Columns = one per test function.
cat("\n--- all_the_tests() output (correct sampler) ---\n")
results_correct <- all_the_tests(
  direct_draws   = direct_correct,
  gibbs_draws    = gibbs_correct,
  test_functions = test_fns
)
print(results_correct)
# Expect: Geweke and KS p-values all large (> 0.05).


# =============================================================================
# 5. Geweke — broken Gibbs chain
# =============================================================================
# Using broken_posterior_sampler: the Gibbs chain never updates mu from the
# data, so it stays near the prior while the direct arm reflects the
# prior-predictive. The test-function distributions will differ noticeably.
# =============================================================================
cat("\n=== Geweke: broken Gibbs ===\n")

# Direct arm is the same prior-predictive regardless of the sampler.
# We can reuse direct_correct.

gibbs_broken <- geweke_suc_cond_draws(
  prior_sampler      = my_prior_sampler,
  likelihood_sampler = my_likelihood_sampler,
  posterior_sampler  = broken_posterior_sampler,
  prior_hyper_params = prior_hyper_params,
  n_draws            = n_geweke,
  n_obs              = n_obs
)

# 5a. QQ plots — should deviate from y = x for data-dependent test functions
geweke_plot(
  direct_draws   = direct_correct,   # reuse direct arm
  gibbs_draws    = gibbs_broken,
  test_functions = test_fns
)
title("Geweke QQ — broken sampler (deviations expected for y_mean, log_like)",
      outer = TRUE, line = -1)

# 5b. all_the_tests() — expect small p-values for data-dependent functions
cat("\n--- all_the_tests() output (broken sampler) ---\n")
results_broken <- all_the_tests(
  direct_draws   = direct_correct,
  gibbs_draws    = gibbs_broken,
  test_functions = test_fns
)
print(results_broken)
# Expect: Geweke and KS p-values near 0 for y_mean and log_like,
#         since those functions depend on y which the broken sampler ignores.


# =============================================================================
# 6. Accessing raw internals for custom analyses
# =============================================================================

# --- SBC internals ---
# Posterior draw matrix for simulation 1 (n_draws x n_params)
post_sim1 <- sbc_correct$sbc_result$fits[[1]]
cat("\nFirst 5 posterior draws for mu, simulation 1:\n")
print(head(post_sim1, 5))

# True parameter value and data for simulation 1
theta_tilde_1 <- sbc_correct$dataset$variables[1, ]
y_1           <- sbc_correct$dataset$generated[[1]]$y
cat(sprintf("\nTrue mu (sim 1): %.4f\n", theta_tilde_1$mu))
cat(sprintf("Data: mean = %.4f, sd = %.4f\n", mean(y_1), sd(y_1)))

# --- Geweke internals ---
# theta_matrix (n_draws x n_params) and y_matrix (n_draws x n_obs) are
# stored on both the direct and Gibbs draw objects.
cat(sprintf("\nGibbs theta_matrix: %d draws x %d params\n",
            nrow(gibbs_correct$theta), ncol(gibbs_correct$theta)))
cat(sprintf("Gibbs y_matrix:     %d draws x %d obs\n",
            nrow(gibbs_correct$y), ncol(gibbs_correct$y)))
