# bayescheckr

**bayescheckr** is an R package for validating Bayesian posterior samplers. It provides multiple complementary validation methods—including **Simulation-Based Calibration (SBC)**, **Geweke's joint-distribution test**, and **distributional-shift diagnostics**—to help ensure that posterior samplers produce draws from the correct target distribution before they are used for inference.

Rather than relying on a single diagnostic, **bayescheckr** offers several independent validation approaches. Different methods are sensitive to different implementation errors, making it easier to identify subtle bugs that may otherwise go unnoticed.

------------------------------------------------------------------------

## Installation

Install the package directly from GitHub:

``` r
devtools::install_github("BayesCheck/bayescheckr")
```

When developing the package locally:

``` r
devtools::load_all(".")
```

# Quick Start

Every validation method in **bayescheckr** relies on three user-supplied functions that define your Bayesian model:

1.  **Prior sampler** – draws parameters from the prior distribution.
2.  **Likelihood sampler** – simulates observations given parameters.
3.  **Posterior sampler** – produces posterior draws given observed data.

Once these three functions are supplied, the same model specification can be used throughout the package for SBC, Geweke diagnostics, and Classification testing.

------------------------------------------------------------------------

## Example model

``` r
library(bayescheckr)

# Hyperparameters passed to all three samplers
prior_hyper_params <- list(
  mu_0    = 0,
  sigma_0 = 1,
  sigma   = 2
)

#----------------------------------------------------------
# Prior sampler
#----------------------------------------------------------

my_prior <- function(hyper) {

  list(
    mu = rnorm(
      1,
      mean = hyper$mu_0,
      sd   = hyper$sigma_0
    )
  )

}

#----------------------------------------------------------
# Likelihood sampler
#----------------------------------------------------------

my_likelihood <- function(n_obs, theta, hyper) {

  rnorm(
    n_obs,
    mean = theta["mu"],
    sd   = prior_hyper_params$sigma
  )

}

#----------------------------------------------------------
# Posterior sampler
#----------------------------------------------------------

my_posterior <- function(
  ndraws,
  y,
  prior_hyper_params,
  init = NULL
) {

  n <- length(y)

  sigma_n_sq <-
    1 /
    (
      n / prior_hyper_params$sigma^2 +
      1 / prior_hyper_params$sigma_0^2
    )

  mu_n <-
    sigma_n_sq *
    (
      sum(y) / prior_hyper_params$sigma^2 +
      prior_hyper_params$mu_0 /
      prior_hyper_params$sigma_0^2
    )

  matrix(
    rnorm(
      ndraws,
      mu_n,
      sqrt(sigma_n_sq)
    ),
    ncol = 1,
    dimnames = list(NULL, "mu")
  )

}
```
## Test functions

Both SBC and Geweke's validation framework only validate on the scalar level. This means validating each model parameter on its own is no issue, but to check the interaction between parameters, we must set up some test functions that transform many parameters into a scalar, for example, multiplying them together.
Our package contains a `make_test_functions` function that makes this process easy and checks if the test functions are valid:

```{r}
test_fns <- make_test_functions(

  mu_marginal      = function(theta, y) theta$mu,
  sigmasq_marginal = function(theta, y) theta$sigmasq,
  product          = function(theta, y) theta$mu * theta$sigma,
  log_like         = function(theta, y) sum(
                                          dnorm(
                                          y, 
                                          theta$mu, 
                                          sqrt(theta$sigmasq), 
                                          log = TRUE))
)

```

------------------------------------------------------------------------

# Simulation-Based Calibration (SBC)

Simulation-Based Calibration (SBC) validates Bayesian posterior samplers by repeatedly simulating parameters from the prior, generating synthetic datasets, and checking whether the true parameter values are uniformly ranked among posterior draws.

If the posterior sampler is correct, every parameter should have a **uniform rank histogram**, and the empirical cumulative distribution function (ECDF) differences should fluctuate around zero.

Running SBC is straightforward:

``` r
result <- run_sbc(

  prior_sampler      = my_prior,
  likelihood_sampler = my_likelihood,
  posterior_sampler  = my_posterior,

  n_sims  = 1000,
  n_obs   = 100,
  n_draws = 10000,

  prior_hyper_params = prior_hyper_params,
  test_fns = test_fns,

  parallelize = TRUE,
  globals = NULL

)
```

`run_sbc()` returns a `bayescheckr_sbc` object containing

- SBC output from the `SBC` package
- simulation metadata
- posterior draws
- true parameter values
- diagnostic summaries from the established test functions
(NOTE: they must include "marginal" test functions if you wish to verify parameters on their own.)

------------------------------------------------------------------------

## Visual diagnostics

Rank histograms should appear approximately uniform.

``` r
SBC::plot_rank_hist(result$sbc_result)
```

The ECDF difference plot should remain close to zero.

``` r
SBC::plot_ecdf_diff(result$sbc_result)
```

Together, these provide a quick visual assessment of posterior calibration.

------------------------------------------------------------------------

## Parallel execution

Large SBC experiments can be computationally intensive.

`run_sbc()` supports parallel execution using the **future** framework. Setting

``` r
parallelize = TRUE
```

automatically distributes independent simulations across available workers.

------------------------------------------------------------------------

# Geweke Joint-Distribution Test

The Geweke test compares two different sampling procedures that should produce draws from the same joint distribution (p(\theta, y)).

If the posterior sampler is implemented correctly, both procedures generate samples from the same stationary distribution, even though they obtain those samples in different ways.

The two procedures are:

- **Marginal–conditional (direct) sampler** – independently draw parameters from the prior, then simulate data from the likelihood.
- **Successive–conditional (Gibbs) sampler** – repeatedly alternate between sampling from the posterior (p(\theta \mid y)) and the likelihood (p(y \mid \theta)).

Agreement between these two sampling procedures provides evidence that the posterior sampler is correct.

------------------------------------------------------------------------

## Generating draws

The marginal–conditional sampler generates independent draws from the joint distribution.

``` r
direct_draws <- geweke_mc_draws(

  prior_sampler      = my_prior,
  likelihood_sampler = my_likelihood,

  prior_hyper_params = prior_hyper_params,

  n_draws = 2000,
  n_obs   = 30

)
```

The successive–conditional sampler repeatedly alternates between posterior sampling and likelihood simulation.

``` r
gibbs_draws <- geweke_sc_draws(

  prior_sampler      = my_prior,
  likelihood_sampler = my_likelihood,
  posterior_sampler  = my_posterior,

  prior_hyper_params = prior_hyper_params,

  n_draws = 2000,
  n_obs   = 30

)
```

------------------------------------------------------------------------

## QQ plots

For every supplied test function, **bayescheckr** compares the two sets of draws using QQ plots.

``` r
plot_geweke_tests(
  direct_draws,
  gibbs_draws,
  test_fns
)
```

If both sampling procedures target the same joint distribution, the plotted points should closely follow the 45° line.

Systematic departures from the diagonal indicate discrepancies between the two distributions and suggest an error in the posterior implementation.

------------------------------------------------------------------------

# Formal Geweke Tests

In addition to visual diagnostics, **bayescheckr** performs several statistical tests.

``` r
results <- tabulate_geweke_tests(
  direct_draws,
  gibbs_draws,
  test_fns
)

print(results)
```

The returned data frame contains one row for each test function and one column for each diagnostic.

| Statistic | Purpose |
|----------------------|--------------------------------------------------|
| Geweke statistic | Difference-in-means z-test comparing the two sampling procedures |
| Geweke p-value | Evidence against equality of expectations |
| KS statistic | Two-sample Kolmogorov–Smirnov test |
| KS p-value | Evidence against equality of distributions |
| Convergence statistic | Geweke convergence diagnostic computed on the Gibbs chain |
| Convergence p-value | Evidence against within-chain convergence |

Ideally,

- QQ plots follow the diagonal,
- Geweke p-values are large,
- KS p-values are large,
- convergence diagnostics show no evidence of nonstationarity.

No single diagnostic should be interpreted in isolation. The combination of graphical and formal tests provides substantially stronger evidence of correctness than any individual statistic.

------------------------------------------------------------------------

#  Validation by classification

Words

------------------------------------------------------------------------

## Running a Classification test

``` r
shift_result <- run_distributional_shift(
  direct_draws = thetaA,
  gibbs_draws = thetaB,
  train_frac = 0.8,
  nrounds = 100
)
```

words

``` r
c2st <- shift_result$influential$c2st
cat(sprintf(
  "C2ST accuracy = %.3f  (n_test = %d)\n  z = %.2f, p = %.3g\n",
  c2st$statistic, c2st$n_test, c2st$z, c2st$p_value
))
```

more words

``` r
print(tabulate_shift_divergences(shift_result))
```

------------------------------------------------------------------------

## Supported classifiers: Not True but might be a good idea to make it do this

Several classification algorithms are available.

| Method | Typical use |
|--------------------|----------------------------------------------------|
| `"logistic"` | Fast baseline diagnostic |
| `"random_forest"` | Captures nonlinear interactions |
| `"xgboost"` | High predictive accuracy for subtle distributional differences |
| `"svm"` | Effective for complex decision boundaries |

Different classifiers may detect different types of posterior discrepancies. When computationally feasible, comparing several classifiers is recommended.

------------------------------------------------------------------------

## Interpreting the results

Suppose a classifier achieves approximately **50% accuracy** on held-out data and has relatively equal feature importance for all features.

This indicates that the classifier cannot reliably distinguish between the two posterior sample sets, providing evidence that the samplers generate similar distributions.

On the other hand, classification accuracy substantially above chance suggests a detectable distributional shift.

While this does not identify the source of the discrepancy, it provides strong evidence that the candidate sampler differs from the reference distribution.

------------------------------------------------------------------------

# Function Reference

## Simulation-Based Calibration

| Function | Description |
|---------------------|---------------------------------------------------|
| `run_sbc()` | Runs Simulation-Based Calibration and returns a `bayescheckr_sbc` object. |
| `recompute_ranks()` | Recomputes SBC ranks using user-defined test functions. |
| `rank_mean_test()` | Tests whether the empirical rank mean equals its theoretical value. |
| `rank_variance_test()` | Tests whether the empirical rank variance matches the uniform expectation. |
| `rank_ks_test()` | Performs a Kolmogorov–Smirnov test for rank uniformity. |
| `rank_kl_divergence()` | Computes the KL divergence between empirical and uniform rank distributions. |

------------------------------------------------------------------------

## Geweke Diagnostics

| Function | Description |
|----------------------|-------------------------------------------------|
| `geweke_mc_draws()` | Generates marginal–conditional (direct) draws. |
| `geweke_sc_draws()` | Runs the successive–conditional Gibbs sampler. |
| `plot_geweke_tests()` | Produces QQ plots for each supplied test function. |
| `tabulate_geweke_tests()` | Computes Geweke, KS, and convergence diagnostics in a single table. |

------------------------------------------------------------------------

## Distributional Shift

| Function | Description |
|------------------------|------------------------------------------------|
| `run_distributional_shift()` | Compares posterior samples using supervised classification. |
| `summary()` | Summarizes classification performance and diagnostic statistics. |
| `plot()` | Visualizes classifier performance and distributional separation. |

------------------------------------------------------------------------

## Shared Utilities

| Function | Description |
|----------------------|--------------------------------------------------|
| `make_test_functions()` | Validates and packages user-defined scalar test functions. |
| `ranks_to_sbc_results()` | Converts recomputed ranks into an SBC-compatible object for plotting. |

# Variable-Dimension Samplers (RJMCMC)

Some Bayesian samplers produce parameter vectors whose dimension changes during sampling.

Since many validation procedures assume a fixed-dimensional parameter space, these samplers require an additional representation before validation.

## Fixed-dimension embedding

Each posterior draw should be embedded into a fixed-dimensional representation corresponding to the largest model considered.

For example,

``` text
Maximum dimension = 5

Model k = 2

(theta1, theta2)

↓

(theta1, theta2, NA, NA, NA)
```

Inactive parameters should be padded using a consistent placeholder (typically `NA`).

Users may also include additional variables such as

- the current model dimension,
- indicator variables describing which parameters are active,
- any other quantities needed for downstream diagnostics.

Once embedded into a fixed-dimensional representation, the resulting samples can be analyzed using the standard **bayescheckr** workflow.

------------------------------------------------------------------------

# Recommended Validation Workflow

Although each diagnostic can be used independently, we recommend combining multiple validation methods.

A typical workflow is

``` text
Prior sampler
      │
      ▼
Likelihood sampler
      │
      ▼
Posterior sampler
      │
      ▼
 Simulation-Based Calibration
      │
      ▼
 Geweke Joint Test
      │
      ▼
 Distributional Shift
```

add statement here
