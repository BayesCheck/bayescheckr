#library(bayeschekr)
# Simulation-Based Calibration example using a Normal-Normal conjugate model.
#
# Model
#   Prior:      mu ~ Normal(mu0, tau^2)          [tau^2 known]
#   Likelihood: y_i ~ Normal(mu, sigma^2)        [sigma^2 known]
#   Posterior:  mu | y ~ Normal(mu_n, tau_n^2)   [analytical]
#
#     where:
#       tau_n^2 = 1 / (1/tau^2 + n/sigma^2)
#       mu_n    = tau_n^2 * (mu0/tau^2 + n*ybar/sigma^2)
#
# The posterior sampler below draws directly from this closed-form posterior,
# so the SBC ranks should be uniform if the implementation is correct.
# ─────────────────────────────────────────────────────────────────────────────

source("simulate.R")   # loads run_simulations(), rank_* tests, etc.

# ── 0. Hyper-parameters ───────────────────────────────────────────────────────

prior_hyper_params <- list(
  mu0    = 0,    # prior mean
  tau2   = 1,    # prior variance
  sigma2 = 4     # likelihood variance (treated as known)
)

# ── 1. Prior sampler ──────────────────────────────────────────────────────────
# Must return a *named* numeric vector.

my_prior_sampler <- function(hyper) {
  mu <- rnorm(1, mean = hyper$mu0, sd = sqrt(hyper$tau2))
  c(mu = mu)
}

# ── 2. Likelihood sampler ─────────────────────────────────────────────────────
# Receives n_obs (integer) and the named draw from the prior.

my_likelihood_sampler <- function(n_obs, theta_tilde) {
  rnorm(n_obs, mean = theta_tilde["mu"], sd = sqrt(prior_hyper_params$sigma2))
}

# ── 3. Posterior sampler ──────────────────────────────────────────────────────
# Must return a matrix with n_draws rows and one named column per parameter.

my_posterior_sampler <- function(ndraws, y, prior_hyper_params) {
  n      <- length(y)
  ybar   <- mean(y)
  mu0    <- prior_hyper_params$mu0
  tau2   <- prior_hyper_params$tau2
  sigma2 <- prior_hyper_params$sigma2

  # Conjugate update
  tau_n2 <- 1 / (1/tau2 + n/sigma2)
  mu_n   <- tau_n2 * (mu0/tau2 + n*ybar/sigma2)

  # Draw from exact posterior
  draws <- rnorm(ndraws, mean = mu_n, sd = sqrt(tau_n2))

  # Return as named matrix (one column per parameter)
  matrix(draws, ncol = 1, dimnames = list(NULL, "mu"))
}

# ── 4. Run simulations ────────────────────────────────────────────────────────

set.seed(42)

sims <- run_simulations(
  prior_sampler      = my_prior_sampler,
  likelihood_sampler = my_likelihood_sampler,
  posterior_sampler  = my_posterior_sampler,
  n_sims             = 200,
  n_obs              = 50,
  n_draws            = 999,    # L; rank ∈ {0, 1, ..., L}
  prior_hyper_params = prior_hyper_params
)

cat("Simulations complete.\n")
cat("  n_sims  :", sims$n_sims,  "\n")
cat("  n_obs   :", sims$n_obs,   "\n")
cat("  n_draws :", sims$n_draws, "\n\n")

# ── 5. Compute ranks manually from the sims object ───────────────────────────
# rank_i = #{posterior draws < theta_tilde}  for each simulation

param_names <- names(sims$sims[[1]]$theta_tilde)   # "mu"

rank_matrix <- matrix(
  NA_integer_,
  nrow     = sims$n_sims,
  ncol     = length(param_names),
  dimnames = list(NULL, param_names)
)

for (i in seq_len(sims$n_sims)) {
  theta_tilde <- sims$sims[[i]]$theta_tilde
  draws       <- sims$sims[[i]]$draws          # [n_draws x n_params] matrix
  for (p in param_names) {
    rank_matrix[i, p] <- sum(draws[, p] < theta_tilde[p])
  }
}

cat("Rank summary for 'mu':\n")
print(summary(rank_matrix[, "mu"]))
cat("\n")

# ── 6. Wrap ranks in an SBC_results-like object expected by the test fns ──────
# The test functions check  inherits(results, "SBC_results")  and read
# results$stats$rank, so we construct a minimal compatible object.

make_sbc_results <- function(ranks_vec) {
  structure(
    list(stats = data.frame(rank = ranks_vec)),
    class = "SBC_results"
  )
}

L <- sims$n_draws   # number of posterior draws = upper bound of rank

results_mu <- make_sbc_results(rank_matrix[, "mu"])

# ── 7. Run the diagnostic tests ───────────────────────────────────────────────

tests <- list(
  mean     = rank_mean_test,
  variance = rank_variance_test,
  ks       = rank_ks_test,
  kl       = rank_kl_divergence
)

test_output <- run_tests(results_mu, L = L, tests = tests)

# ── 8. Report results ─────────────────────────────────────────────────────────

cat("═══════════════════════════════════════\n")
cat("  SBC Diagnostic Results  —  mu\n")
cat("═══════════════════════════════════════\n\n")

cat("── Mean test ──\n")
cat(sprintf("  Observed mean : %.3f\n", test_output$mean$observed_mean))
cat(sprintf("  Expected mean : %.3f\n", test_output$mean$expected_mean))
cat(sprintf("  Difference    : %+.3f\n\n", test_output$mean$difference))

cat("── Variance test ──\n")
cat(sprintf("  Observed var  : %.3f\n", test_output$variance$observed_variance))
cat(sprintf("  Expected var  : %.3f\n", test_output$variance$expected_variance))
cat(sprintf("  Ratio         : %.3f  (1.0 = perfect)\n\n",
            test_output$variance$ratio))

cat("── KS test ──\n")
cat(sprintf("  D statistic   : %.4f\n", test_output$ks$statistic))
cat(sprintf("  p-value       : %.4f  (%s)\n\n",
            test_output$ks$p.value,
            ifelse(test_output$ks$p.value > 0.05, "PASS", "FAIL")))

cat("── KL divergence ──\n")
cat(sprintf("  KL(obs||unif) : %.4f  (0 = perfect)\n\n",
            test_output$kl$KL_divergence))

# ── 9. Quick rank histogram ───────────────────────────────────────────────────

cat("── Rank histogram (mu) ──\n")
hist(
  rank_matrix[, "mu"],
  breaks = 20,
  main   = "SBC Rank Histogram — mu",
  xlab   = paste0("Rank  (0–", L, ")"),
  col    = "steelblue",
  border = "white"
)
abline(h = sims$n_sims / 20, col = "red", lty = 2, lwd = 2)   # expected flat line
legend("topright", legend = "Expected (uniform)", col = "red", lty = 2, lwd = 2)


#NOTE: CHANGE ME WHEN THE PACKAGE IS COMPLETE SO THAT THE EXAMPLE IS ACTUALY INDICATIVE OF WHAT THE CODE DOES.
