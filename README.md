------------------------------------------------------------------------

# bayescheckr

**bayescheckr** is an R package for validating Bayesian posterior samplers. It several approaches so you can be confident your sampler is drawing from the right posterior before you use it for inference.

------------------------------------------------------------------------

## Installation

``` r
# Install from GitHub (once public)
devtools::install_github("BayesCheck/bayescheckr")

# Or load from source during development
devtools::load_all(".")
```

**Dependencies:** `SBC`, `posterior`, `coda`, `future`, `tidyverse`

------------------------------------------------------------------------

## Quick start

The three things you always need to provide are a **prior sampler**, a **likelihood sampler**, and a **posterior sampler**. 
You might also want to set some prior hyperparameters. `bayescheckr` does the rest.

``` r
library(bayescheckr)

# Hyperparameters passed to all three samplers
prior_hyper_params <- list(mu_0 = 0, sigma_0 = 1, sigma = 2)

# 1. Prior: returns a *named* vector — names become parameter labels
my_prior <- function(hyper) {
  c(mu = rnorm(1, mean = hyper$mu_0, sd = hyper$sigma_0))
}

# 2. Likelihood: (n_obs, theta) -> returns a matrix or list (best to make output general).
my_likelihood <- function(n, theta) {
  rnorm(n, mean = theta["mu"], sd = prior_hyper_params$sigma)
}

# 3. Posterior: (ndraws, y, prior_hyper_params) -> ndraws x n_params matrix
my_posterior <- function(ndraws, y, prior_hyper_params) {
  n         <- length(y)
  sigma_n_sq <- 1 / (n / prior_hyper_params$sigma^2 + 1 / prior_hyper_params$sigma_0^2)
  mu_n       <- sigma_n_sq * (sum(y) / prior_hyper_params$sigma^2 +
                              prior_hyper_params$mu_0 / prior_hyper_params$sigma_0^2)
  matrix(rnorm(ndraws, mu_n, sqrt(sigma_n_sq)), ncol = 1, dimnames = list(NULL, "mu"))
}
```

------------------------------------------------------------------------

## SBC

SBC runs your sampler across many simulated datasets and checks that the true parameter ranks uniformly among the posterior draws. If your sampler is correct, the rank histogram should be flat.

``` r
result <- run_sbc(
  prior_sampler      = my_prior,
  likelihood_sampler = my_likelihood,
  posterior_sampler  = my_posterior,
  n_sims             = 300,
  n_obs              = 30,
  n_draws            = 500,
  prior_hyper_params = prior_hyper_params
)

# Rank histogram — should be flat
SBC::plot_rank_hist(result$sbc_result)

# ECDF difference — should hug zero
SBC::plot_ecdf_diff(result$sbc_result)
```

### Test functions

By default SBC checks raw parameter ranks. You can also check ranks of *derived quantities* using test functions — useful for catching bugs that only show up in nonlinear summaries of the posterior.

``` r
test_fns <- make_test_functions(
  mu_sq    = function(theta, y) theta["mu"]^2,
  log_like = function(theta, y)
    sum(dnorm(y, theta["mu"], prior_hyper_params$sigma, log = TRUE))
)

# Recompute ranks on test functions rather than raw parameters
ranked  <- recompute_ranks(result, test_fns)
shell   <- ranks_to_sbc_results(ranked)

SBC::plot_rank_hist(shell)
SBC::plot_ecdf_diff(shell)
```

------------------------------------------------------------------------

## Geweke test

The Geweke test compares two samplers that should share the same stationary distribution:

- **Direct (marginal-conditional)** draws: sample θ from the prior, then y from the likelihood. These are i.i.d. draws from the joint distribution p(θ, y).
- **Successive-conditional (Gibbs)** draws: starting from any initial value, alternate between drawing θ \| y from the posterior and y \| θ from the likelihood. At stationarity these are also draws from p(θ, y).

If your posterior sampler is correct, both arms have the same marginal distribution for every test function h(θ, y). The QQ plots should fall on the y = x diagonal and the formal tests should return large p-values.

``` r
direct_draws <- geweke_mc_draws(
  prior_sampler      = my_prior,
  likelihood_sampler = my_likelihood,
  prior_hyper_params = prior_hyper_params,
  n_draws            = 2000,
  n_obs              = 30
)

gibbs_draws <- geweke_sc_draws(
  prior_sampler      = my_prior,
  likelihood_sampler = my_likelihood,
  posterior_sampler  = my_posterior,
  prior_hyper_params = prior_hyper_params,
  n_draws            = 2000,
  n_obs              = 30
)

# QQ plots — points should lie on y = x
plot_geweke_tests(direct_draws, gibbs_draws, test_fns)
```

### Formal tests

`tabulate_geweke_tests()` returns a data frame with three tests per test function:

``` r
results <- tabulate_geweke_tests(direct_draws, gibbs_draws, test_fns)
print(results)
```

| Row | What it tests |
|------------------------------------|------------------------------------|
| `Geweke statistic / p_value` | Difference-in-means z-test (Geweke 2004) using long-run variance for the Gibbs arm |
| `KS statistic / p_value` | Kolmogorov–Smirnov two-sample test |
| `Convergence statistic / p_value` | `coda::geweke.diag` on the Gibbs chain — checks internal chain convergence |

Large p-values (\> 0.05) across all three indicate no evidence of miscalibration.

------------------------------------------------------------------------

## Function reference

### SBC

| Function | Description |
|------------------------------------|------------------------------------|
| `run_sbc()` | Main entry point. Runs SBC and returns a `bayescheckr_sbc` object |
| `recompute_ranks()` | Recomputes ranks using user-supplied test functions |
| `rank_mean_test()` | Tests whether the rank mean equals the expected value `L/2` |
| `rank_variance_test()` | Tests whether rank variance matches the uniform expectation |
| `rank_ks_test()` | KS test of rank uniformity |
| `rank_kl_divergence()` | KL divergence between empirical and uniform rank distribution |

### Geweke

| Function | Description |
|------------------------------------|------------------------------------|
| `geweke_mc_draws()` | Generates direct (marginal-conditional) draws |
| `geweke_sc_draws()` | Runs the successive-conditional (Gibbs) chain |
| `tabulate_geweke_tests()` | Returns dataframe with Geweke z-test, KS test, and convergence diagnostics in one call |
| `plot_geweke_tests()` | QQ plots for each test function |

### Shared

| Function | Description |
|------------------------------------|------------------------------------|
| `make_test_functions()` | Validates and packages user-supplied test functions |

------------------------------------------------------------------------

## Sampler interface

All three samplers must follow these exact signatures:

``` r
# Prior: no arguments other than hyperparameters; returns a named vector
prior_sampler <- function(prior_hyper_params) { ... }

# Likelihood: (n_obs, theta) -> numeric vector of length n_obs,
#             or matrix of dimensions n_obs x p for multivariate observations
likelihood_sampler <- function(n_obs, theta) { ... }

# Posterior: (ndraws, y, prior_hyper_params) -> ndraws x n_params matrix
#            column names must match the names returned by prior_sampler
posterior_sampler <- function(ndraws, y, prior_hyper_params) { ... }
```

Test functions must take exactly `(theta, y)` and return a scalar:

``` r
# theta: named numeric vector (one draw)
# y:     numeric vector of observations
my_test <- function(theta, y) { ... }  # must return a single number
```
