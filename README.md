# bayescheckr

**bayescheckr** is an R package for validating Bayesian posterior samplers. It provides multiple complementary validation methods—including [**Simulation-Based Calibration (SBC)**](https://github.com/hyunjimoon/SBC), [**Geweke's joint-distribution test**](https://www.tandfonline.com/doi/abs/10.1198/016214504000001132), and [**distributional-shift diagnostics**](https://arxiv.org/abs/1610.06545)—to help ensure that posterior samplers produce draws from the correct target distribution before they are used for inference.

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
  mu_0 = 0,
  sigma2_0 = 1,
  sigma2 = 2
)

#----------------------------------------------------------
# Prior sampler
#----------------------------------------------------------

my_prior <- function(prior_hyper_params) {

  return(list(
    mu = rnorm(
      1,
      mean = hyper$mu_0,
      sd = hyper$sigma2_0
    )
  ))

}

#----------------------------------------------------------
# Likelihood sampler
#----------------------------------------------------------

my_likelihood <- function(n_obs, theta, prior_hyper_params) {

  return(rnorm(
    n_obs,
    mean = theta["mu"],
    sd = prior_hyper_params$sigma2
  ))

}

#----------------------------------------------------------
# Posterior sampler
#----------------------------------------------------------

my_posterior <- function(
  n_draws,
  y,
  prior_hyper_params,
  init = NULL
) {

  n <- length(y)

  sigma2_n <-
    1 /
    (
      n / prior_hyper_params$sigma^2 +
      1 / prior_hyper_params$sigma2_0^2
    )

  mu_n <-
    sigma2_n *
    (
      sum(y) / prior_hyper_params$sigma^2 +
      prior_hyper_params$mu_0 /
      prior_hyper_params$sigma_0^2
    )

  return(matrix(
    rnorm(
      ndraws,
      mu_n,
      sqrt(sigma2_n)
    ),
    ncol = 1,
    dimnames = list(NULL, "mu")
  ))

}
```
# How to format the samplers

To work well with our package, format the inputs and outputs of the samplers in a certain way. 

The prior: 
- **input**: `prior_hyper_params`
- **return**: a named object (list or named vector) where `names()` gives parameter names and `length()` gives parameter count.

Other notes: The example above wraps scalar parameters in a simple list. To include vector parameters, the recommended pattern is as follows: `c(setNames(vector_params, names), list(scalar_param = value))`.

```{r}
#for parameters beta (a vector) and sigmasq (a scalar)

beta_list <- setNames(as.list(beta), paste0("beta", seq_along(beta)))
return(c(beta_list, list(sigmasq = sig2)))
```

The likelihood: 
- **input**: `(n_obs, theta, prior_hyper_params)`
- **return**: a plain numeric vector of length n_obs.

Other notes: Access theta defensively. Use theta[["name"]] or theta[[index]] with a double bracket not single bracket, and `unlist(theta[1:p])` to extract sub-vectors.

The posterior:
- **input**: `(n_draws, y, prior_hyper_params, init = NULL)`. The `init` argument allows the posterior to initialize to better and better update points as the MCMC chain converges. First line of the function body should **always be** `theta <- if (is.null(init)) <your_default> else unlist(init)`. 
- **return**: a matrix with `n_draws` rows and named columns matching the names from `prior_sampler`.

## Test functions

Both SBC and Geweke's validation framework only validate on the scalar level. This means validating each model parameter on its own is no issue, but to check the interaction between parameters, we must set up some test functions that transform many parameters into a scalar, for example, multiplying them together.

Test functions must always access parameters within the inputted `theta` using `[["parameter"]]` or `$parameter`.

Our package contains a `make_test_functions` function that makes this process easy and checks if the test functions are valid:

```{r}
test_fns <- make_test_functions(

  mu_marginal = function(theta, y) theta$mu,
  sigmasq_marginal = function(theta, y) theta$sigmasq,
  product = function(theta, y) theta$mu * theta$sigma,
  log_like = function(theta, y) sum(dnorm(y, 
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

NOTE TO SELF TRISHA: make all Geweke functions intake draws in the same format (either one object or 2, the gibbs/direct draws separately)

The Geweke test compares two different sampling procedures that should produce draws from the same joint distribution ($p(\theta, y)$).

If the posterior sampler is implemented correctly, both procedures generate samples from the same stationary distribution, even though they obtain those samples in different ways.

The two procedures are:

- **Marginal–conditional (direct) sampler** – independently draw parameters from the prior, then simulate data from the likelihood.
- **Successive–conditional (Gibbs) sampler** – repeatedly alternate between sampling from the posterior ($p(\theta \mid y)$) and the likelihood ($p(y \mid \theta)$).

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

The classifier 2-sample test (C2ST), another validation method, operates on the same principle as the Geweke comparison method: if the posterior sampler is correct, then draws from it should come from the same distribution as joint direct draws (using marginal-conditional samplers). A classifier trained on the two sets of draws should then ideally be unable to tell them apart and thus have around 50% classification accuracy. Further, no parameter should significantly contribute to the classifier’s performance if all are sampled and updated correctly. 

Our package offers a single method, `run_distributional_shift()`, that generates the draws to compare, runs a classifier of your choice, and outputs the classifier accuracy and the relative contributions of each parameter (feature importance) to that accuracy. 

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

`run_distributional_shift()` trains a classifier to distinguish direct draws from Gibbs draws, using `theta` and `y`jointly. It splits each set into a training and testing group according to `train_frac`, fits the specified `classifier` (logistic regression set as default; xgboost is available for high-dimensional data), and evaluates classification accuracy on the held-out testing set. It returns a list with three elements: 

- `classifier_accuracy`: test-set accuracy; 
- `c2st`: a list containing the test statistic, sample size, standard error under the null, z-score, p-value, and Reject/Fail-to-Reject decision at the given `alpha`; 
- `feature_importance`: a data frame of per-feature permutation importance, giving each feature's observed accuracy drop and permutation-based p-value. 

Under a correctly implemented sampler, `classifier_accuracy` should sit close to 0.5 and the C2ST p-value should exceed `alpha`; a low p-value indicates the classifier can distinguish the direct draws from the Gibbs reliably well, which is potential evidence of an erroneous sampler. 

``` r
c2st <- shift_result$influential$c2st
cat(sprintf(
  "C2ST accuracy = %.3f  (n_test = %d)\n  z = %.2f, p = %.3g\n",
  c2st$statistic, c2st$n_test, c2st$z, c2st$p_value
))
```

You can use the $ operator to access the c2st dataframe outputted by `run_distributional_shift()`. It gives results for the test *overall*; if you are interested in determining the relative importance of any one parameter in a model, access the feature importance data frame instead.

``` r
c2st <- shift_result$feature_importance
```
In this data frame, each row represents one parameter (e.g., “beta1” or “mu”) and lists specific statistics for each: `observed_drop` and a corresponding `p_value` ($H_0$: the parameter does not contribute much to classifier accuracy). `observed_drop` is computed as the drop in classifier accuracy from leaving each parameter out iteratively; this determines each parameter’s relative influence. A high `observed_drop` and low `p_value` might indicate that a parameter is sampled incorrectly in the posterior, nudging a researcher to start their debugging there.   

------------------------------------------------------------------------

## Supported classifiers

Two classification algorithms are currently available, with additional methods planned.

| Method            | Typical use                                                     | Status             |
|-------------------|------------------------------------------------------------------|--------------------|
| `"logistic"`      | Fast baseline diagnostic                                        | Available          |
| `"xgboost"`       | High predictive accuracy for subtle distributional differences  | Available           |
| `"random_forest"` | Captures nonlinear interactions                                 | Not yet available  |
| `"svm"`           | Effective for complex decision boundaries                       | Not yet available  |

Different classifiers may detect different types of posterior discrepancies. When computationally feasible, comparing several classifiers is recommended. Any classifier can be supplied directly via the `classifier` argument by implementing the `train`/`predict` interface described below, independent of this table.

The classifier you define must return a list of two functions, `train` and `predict`. The `train` function must intake parameters `(X, y)` and output a fitted model object. The `predict` function must intake parameters `(fit, newX)` with `fit` being `train`'s output and `newX` being the testing data.

```{r}
my_classifier <- function() {
  list(
    train = function(X, y) {
      #run and output the trained model
    },
    predict = function(fit, newX) {
      #run and output predictions
    }
  )
}
```

------------------------------------------------------------------------

## Interpreting the results

Suppose a classifier achieves approximately **50% accuracy** on held-out data and has relatively equal feature importance for all features.

This indicates that the classifier cannot reliably distinguish between the two posterior sample sets, providing evidence that the samplers generate similar distributions.

On the other hand, classification accuracy substantially above chance suggests a detectable distributional shift.

While this does not identify the source of the discrepancy, it provides strong evidence that the candidate sampler differs from the reference distribution.

------------------------------------------------------------------------
# Visualization

Besides the SBC plots and Geweke's QQ plots, all diagnostics so far are numerical, including all of those pertaining to distributional shift.

Numerical summaries can be difficult to interpret directly, particularly for models with many parameters or a trans-dimensional structure (e.g. reversible-jump samplers), where the raw parameter space cannot be inspected visually. In these settings, projecting the classifier's input space into two dimensions can make patterns of separation between direct and Gibbs draws visible that would otherwise only be apparent from a table of numbers.

## `fit_c2st_for_viz()`

`fit_c2st_for_viz()` refits a classifier on a train/test split of `direct_draws$theta` and `gibbs_draws$theta`, mirroring the internal logic of `run_distributional_shift()`, but returns the fitted classifier, the held-out feature matrix, the true labels, and the predicted probabilities rather than only the summary statistics.

This function is provided for computational convenience — it avoids re-implementing the same train/test split and classifier call already used internally by `run_distributional_shift()`. It does not introduce a dependency on C2ST or on any classifier for the purpose of visualization itself: the resulting embeddings below are generated directly from the held-out `theta` values and the true draw labels (`ytest`), both of which are available directly from Geweke's output, independent of whether a classifier is ever fit. 
`fit_c2st_for_viz()` is used only as a convenient source of an already-split test set.

```r
fit_c2st_for_viz(theta_A, theta_B, classifier = NULL, train_frac = 0.8, seed = NULL)
```

**Returns** a list with:

- `fit` — the fitted classifier object.
- `Xtest` — the held-out feature data.frame (theta values from the test split).
- `ytest` — true labels for the test set (`0` = direct, `1` = gibbs).
- `p_test` — the classifier's predicted probability of class 1 (gibbs) for each test row.

`Xtest` and `ytest` are the two elements needed to generate a visualization; `fit` and `p_test` are made available for optional diagnostics such as coloring an embedding by predicted probability instead of true source, or plotting a decision boundary.

## Example: t-SNE embedding

```r
library(Rtsne)
library(ggplot2)

viz <- fit_c2st_for_viz(direct_draws$theta, gibbs_draws$theta,
                        classifier = default_xgb_classifier(), seed = 1)

X_mat <- as.matrix(viz$Xtest)

tsne_out <- Rtsne(X_mat, dims = 2, perplexity = 30, check_duplicates = FALSE)

tsne_df <- data.frame(
  tsne1 = tsne_out$Y[, 1],
  tsne2 = tsne_out$Y[, 2],
  source = factor(viz$ytest, levels = c(0, 1), labels = c("direct", "gibbs"))
)

ggplot(tsne_df, aes(x = tsne1, y = tsne2, color = source)) +
  geom_point(alpha = 0.4, size = 0.8) +
  theme_minimal()
```

Under a correctly implemented sampler, direct and Gibbs draws should be visually indistinguishable in the embedding. Visible clustering by `source` indicates a distributional shift between the two draw sets, consistent with a C2ST rejection.

t-SNE is one option among several. `Xtest` can be passed to any dimensionality-reduction or visualization method that accepts a numeric matrix — UMAP and PaCMAP are two alternatives worth considering, particularly for models with many parameters, since results can vary by technique and comparing more than one is recommended when computationally feasible.

------------------------------------------------------------------------

# Function Reference

## Simulation-Based Calibration

| Function | Description |
|---------------------|---------------------------------------------------|
| `run_sbc()` | Runs Simulation-Based Calibration and returns a `bayescheckr_sbc` object. |
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
| `default_glm_classifier()` | Logistic regression classifier - set as default in `run_distributional_shift()`. |
| `default_xgb_classifier()` | XG Boost classifier - input as `classifier` argument in `run_distributional_shift()`. |
| `summary()` | Summarizes classification performance and diagnostic statistics. |
| `plot()` | Visualizes classifier performance and distributional separation. |
| `fit_c2st_for_viz()` | Creates list object for visualization. | 

------------------------------------------------------------------------

## Shared Utilities

| Function | Description |
|----------------------|--------------------------------------------------|
| `make_test_functions()` | Validates and packages user-defined scalar test functions. |
| `ranks_to_sbc_results()` | Converts recomputed ranks into an SBC-compatible object for plotting. |

# Variable-Dimension Samplers (RJMCMC)

Some Bayesian samplers produce parameter vectors whose dimension changes during sampling.

Since many validation procedures assume a fixed-dimensional parameter space, these samplers require an additional representation before validation. Each posterior draw should be embedded into a fixed-dimensional representation corresponding to the largest model considered.

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
Define Prior sampler
      │
      ▼
Define Likelihood sampler
      │
      ▼
Define Posterior sampler
      │
      ▼
Simulation-Based Calibration
      │
      ▼
Geweke Joint Test
      │
      ▼
Distributional Shift
      │
      ▼
Visualization
 
```

Thank you for engaging with our package! We hope it helps you in your endeavors.
If you have questions, feedback, or encounter any issues, we would love to hear and help you out! Please contact:

```text
rafael.mautner@duke.edu
trisha.iyer@duke.edu
lance.browne@duke.edu
```
