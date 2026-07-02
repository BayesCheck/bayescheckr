# DistributionShift.R
#
# Implements "What is Different Between These Datasets?" (Babbar, Guo & Rudin,
# JMLR 26(180), 2025) as a bayescheckr module.
#
# Sits in R/ alongside sbc.R, geweke.R, and testfunctions.R.
# Takes the same draw objects that geweke_mc_draws() / geweke_sc_draws() return.
#
# Public functions
# ----------------
#   run_distributional_shift()     main entry point; returns a named list
#   tabulate_shift_divergences()   Method 2 — numeric divergence table
#   plot_shift_densities()         Method 2 — marginal density overlays
#
# Dependencies: e1071, ggplot2, tidyr, patchwork


# ------------------------------------------------------------------------------
# Unexported helpers
# ------------------------------------------------------------------------------


.kl_kde <- function(x, y, n_grid = 512) {
  r  <- range(c(x, y))
  lo <- r[1] - diff(r) * 0.1
  hi <- r[2] + diff(r) * 0.1
  pA <- density(x, from = lo, to = hi, n = n_grid)$y + 1e-10
  pB <- density(y, from = lo, to = hi, n = n_grid)$y + 1e-10
  pA <- pA / sum(pA)
  pB <- pB / sum(pB)
  sum(pA * log(pA / pB))
}


.js_kde <- function(x, y, n_grid = 512) {
  r  <- range(c(x, y))
  lo <- r[1] - diff(r) * 0.1
  hi <- r[2] + diff(r) * 0.1
  pA <- density(x, from = lo, to = hi, n = n_grid)$y + 1e-10
  pB <- density(y, from = lo, to = hi, n = n_grid)$y + 1e-10
  pA <- pA / sum(pA)
  pB <- pB / sum(pB)
  M  <- (pA + pB) / 2
  0.5 * sum(pA * log(pA / M)) + 0.5 * sum(pB * log(pB / M))
}


.wass_1d <- function(x, y) {
  n  <- max(length(x), length(y))
  xs <- quantile(x, probs = seq(0, 1, length.out = n))
  ys <- quantile(y, probs = seq(0, 1, length.out = n))
  mean(abs(xs - ys))
}


# ------------------------------------------------------------------------------
# Method 1 — influential-example explanations (classifier two-sample test)
# ------------------------------------------------------------------------------
#
# The A-vs-B classifier accuracy doubles as a Classifier Two-Sample Test
# (C2ST; Lopez-Paz & Oquab, ICLR 2017). Under H0: A and B are drawn from the
# same distribution, the held-out accuracy t_hat is approximately
# Normal(1/2, 1/(4*n_te)). We report the accuracy as the C2ST statistic, a
# z-score against that null, and a one-sided p-value for "accuracy > chance".


.run_influential <- function(theta_A,
                             theta_B,
                             varnames,
                             train_frac = 0.8,
                             nrounds = 100){

  ##-----------------------
  ## Split A
  ##-----------------------

  idx_A <- sample(seq_len(nrow(theta_A)),
                  floor(train_frac*nrow(theta_A)))

  A_train <- theta_A[idx_A,,drop=FALSE]
  A_test  <- theta_A[-idx_A,,drop=FALSE]

  ##-----------------------
  ## Split B
  ##-----------------------

  idx_B <- sample(seq_len(nrow(theta_B)),
                  floor(train_frac*nrow(theta_B)))

  B_train <- theta_B[idx_B,,drop=FALSE]
  B_test  <- theta_B[-idx_B,,drop=FALSE]



  ############################################################
  ## Train A vs B
  ############################################################

  Xtrain <- rbind(A_train,B_train)
  ytrain <- c(rep(0,nrow(A_train)),
              rep(1,nrow(B_train)))

  dtrain <- xgboost::xgb.DMatrix(
    as.matrix(Xtrain),
    label=ytrain
  )

  fit_AB <- xgboost::xgb.train(
    params=list(
      objective="binary:logistic",
      eval_metric="logloss",
      max_depth=4,
      eta=.1,
      subsample=.8,
      colsample_bytree=.8
    ),
    data=dtrain,
    nrounds=nrounds,
    verbose=0
  )



  ############################################################
  ## Test
  ############################################################

  Xtest <- rbind(A_test,B_test)
  ytest <- c(rep(0,nrow(A_test)),
             rep(1,nrow(B_test)))

  n_te <- length(ytest)

  p <- predict(
    fit_AB,
    xgboost::xgb.DMatrix(as.matrix(Xtest))
  )

  pred <- ifelse(p>.5,1,0)

  accuracy <- mean(pred==ytest)



  ############################################################
  ## Classifier two-sample test (C2ST) statistics
  ##
  ## t_hat = accuracy is the C2ST statistic. Under H0 ("A and B come from
  ## the same distribution"), n_te * t_hat ~ Binomial(n_te, 1/2), which for
  ## large n_te is approximated by Normal(1/2, 1/(4*n_te)) (Lopez-Paz &
  ## Oquab, 2017, Sec. 3.1). We use that null to build a z-score and a
  ## one-sided p-value against "accuracy > chance".
  ############################################################

  c2st_se      <- sqrt(1 / (4 * n_te))
  c2st_z       <- (accuracy - 0.5) / c2st_se
  c2st_pvalue  <- stats::pnorm(c2st_z, lower.tail = FALSE)

  c2st <- list(
    statistic = accuracy,
    n_test    = n_te,
    se_null   = c2st_se,
    z         = c2st_z,
    p_value   = c2st_pvalue
  )

  list(

    classifier_accuracy = accuracy,

    c2st = c2st

  )
}

##Trisha's code

#label: load-packages #!pending: need to add these dependencies

library(tidymodels)
library(tidyverse)
library(bayescheckr)

#label: make-dataset

#make the dataset - each should have a label with distribution 0 or 1
# draws from the same distribution so x or z or whatever it's called = draw from dist.
# then column label = 0 or 1 depending on which set of draws it belongs to
# shuffle the order of the rows randomly to get the final full dataset to split

alpha <- 0.05
n <- 5e3

#sampler draws - will eventually add in more complex posterior samplers
#like sc versus mc geweke draws
set0 <- rnorm(n/2)
set1 <- rnorm(n/2)

set0 <- ungroup(geweke_mc_draws(prior_sampler = linreg_prior_sampler,
                                likelihood_sampler = linreg_likelihood_sampler,
                                #posterior_sampler = linreg_posterior_sampler,
                                n_draws = n/2,
                                n_obs = n_obs,
                                prior_hyper_params = prior_hyper_params)$wide)

set1 <- ungroup(geweke_sc_draws(prior_sampler = linreg_prior_sampler,
                                likelihood_sampler = linreg_likelihood_sampler,
                                posterior_sampler = linreg_posterior_sampler,
                                n_draws = n/2,
                                n_obs = n_obs,
                                prior_hyper_params = prior_hyper_params)$wide)

#wide version of Geweke draw outputs

#combine the draws into dataframe

# log_data <- data.frame(z = c(set0, set1),
#                        label = as.factor(c(rep(0, n/2), rep(1, n/2)))
#                        )

log_data <- bind_rows(set0, set1, .id = "label") |>
  mutate(label = as.factor(as.numeric(label) - 1)) |>
  select(-sim_id)

log_data <- log_data[sample(nrow(log_data)), ] #shuffle rows

#label: split-dataset

set.seed(701) #make this a parameter
log_split <- initial_split(log_data, prop = 0.8, #make this a parameter
                           strata = label)

data_tr <- training(log_split)
data_te <- testing(log_split)
n_te <- nrow(data_te)

#label: train-logreg

#fit the model
log_fit <- logistic_reg() |>
  #fit(label ~ z, data = data_tr)
  fit(label ~ ., data = data_tr)

#display tibble
tidy(log_fit)

#test-logreg-and-test-statistic

data_aug <- data_te |>
  mutate(prob_label1 = predict(log_fit,
                               new_data = data.frame(data_te),
                               type = "prob")$.pred_1
  ) |>
  mutate(indicator_inner = as.numeric(prob_label1 > 0.5)) |>
  mutate(indicator_outer = as.numeric(indicator_inner == label))
#make this robust to any kind of classifier output
#all you need are 2 columns - prob(label = 1) and then actual label
#need to add in model accuracy and feature importance analysis - could this reveal flawed dependencies in s & theta?

#tibble(data_aug) #display

data_aug |>
  summarize(t_stat = mean(indicator_outer, na.rm = TRUE), #get t-statistic
            p_value = 2*(1 - pnorm(abs(t_stat - 0.5) + 0.5, #get p-value
                                   mean = 0.5,
                                   sd = sqrt(1/(4*n_te))
            )),
            alpha = alpha, #display significance level #make alpha a parameter
            h0 = ifelse(p_value <= alpha, "Reject", "Fail to Reject")
  )


#input: output of classifier run on testing data - need real labels and the prob(label1)
#param: data_input = the data object from the classifier
#param: correct_label = name of the column that's the label (e.g., here, label)
#param: label_prob = name of column with probability of label 1 (here, prob_label1)
#param: alpha = significance level desired
#param: model_fit = log_fit (if a GLM) that we use to get feature importance calcs
#param: is_glm = used to determine type of feature importance (glm yes/no)

tabulate_classifier_tests <- function(data_input, #a dataframe or list
                                      correct_label, #a string
                                      label_prob, #a string
                                      alpha = 0.05, #a number: sig. level
                                      model_fit = NULL, #log_fit or some other
                                      is_glm #a logical: TRUE/FALSE
                                      ) {

  n_te <- nrow(data_input)

  #get the C2ST test statistic
  data_tstat <- data_input |>
    mutate(correct_label = .data[[correct_label]],
           label_prob = .data[[label_prob]],
           indicator_inner = label_prob > 0.5,
           indicator_outer = as.numeric(indicator_inner == correct_label)) |>
    summarize(t_stat = mean(indicator_outer, na.rm = TRUE), #t-statistic
              p_value = 2*(1 - pnorm(abs(t_stat - 0.5) + 0.5, #2-sided p-value
                                     mean = 0.5,
                                     sd = sqrt(1/(4*n_te))
              )),
              alpha = alpha, #significance level
              h0 = ifelse(p_value <= alpha, "Reject", "Fail to Reject")
    )

  #add the feature importance
  data_features <- data.frame(feature = tidy(model_fit)$term,
                              z_score = tidy(model_fit)$statistic,
                              p_value = tidy(model_fit)$p.value)

  #return t-stat and feature importance dataframes together
  return(list(t_statistic = data_tstat,
       feature_importance = data_features))
}

# ------------------------------------------------------------------------------
# Method 2 — divergence metrics
# ------------------------------------------------------------------------------


.run_divergences <- function(theta_A, theta_B, varnames) {

  rows <- lapply(varnames, function(var) {
    x      <- theta_A[, var]
    y      <- theta_B[, var]
    ks_res <- ks.test(x, y)

    data.frame(
      variable      = var,
      KL_A_to_B     = .kl_kde(x, y),
      JS_divergence = .js_kde(x, y),
      Wasserstein1  = .wass_1d(x, y),
      KS_statistic  = unname(ks_res$statistic),
      KS_pvalue     = ks_res$p.value,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}


# ------------------------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------------------------

#' @export
run_distributional_shift <- function(
    direct_draws,
    gibbs_draws,
    train_frac = 0.8,
    nrounds = 100
) {

  theta_A  <- direct_draws$theta
  theta_B  <- gibbs_draws$theta
  varnames <- colnames(theta_A)

  if (!identical(varnames, colnames(theta_B))) {
    stop("direct_draws$theta and gibbs_draws$theta must have the same column names.")
  }

  message("Method 1: training A-vs-B classifier ...")

  m1 <- .run_influential(
    theta_A,
    theta_B,
    varnames,
    train_frac = train_frac,
    nrounds = nrounds
  )

  message("Method 2: computing divergence metrics ...")

  m2 <- .run_divergences(
    theta_A,
    theta_B,
    varnames
  )

  list(
    influential = m1,
    divergences = m2,
    varnames    = varnames,
    n_A         = nrow(theta_A),
    n_B         = nrow(theta_B)
  )
}


# ------------------------------------------------------------------------------
# Tabulate divergences — Method 2 numeric output
# ------------------------------------------------------------------------------

#' @export
tabulate_shift_divergences <- function(shift_result) {
  shift_result$divergences
}
# ------------------------------------------------------------------------------
# Plot functions — each returns a ggplot object
# ------------------------------------------------------------------------------

# NOTE:
# plot_shift_importance() has been removed because feature importance
# is no longer computed. Method 1 now reports only the classifier
# accuracy and C2ST statistics.

#' @export
plot_shift_densities <- function(direct_draws,
                                 gibbs_draws,
                                 shift_result,
                                 label_A = "Direct (MC)",
                                 label_B = "Gibbs (SC)") {

  theta_A  <- direct_draws$theta
  theta_B  <- gibbs_draws$theta
  varnames <- shift_result$varnames

  plots <- lapply(varnames, function(var) {

    df <- data.frame(
      value      = c(theta_A[, var], theta_B[, var]),
      population = rep(
        c(label_A, label_B),
        c(nrow(theta_A), nrow(theta_B))
      )
    )

    ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = value,
        fill = population,
        colour = population
      )
    ) +
      ggplot2::geom_density(
        alpha = 0.3,
        linewidth = 0.8
      ) +
      ggplot2::scale_fill_manual(
        values = setNames(
          c("#378ADD", "#D85A30"),
          c(label_A, label_B)
        )
      ) +
      ggplot2::scale_colour_manual(
        values = setNames(
          c("#185FA5", "#993C1D"),
          c(label_A, label_B)
        )
      ) +
      ggplot2::labs(
        title = var,
        x = NULL,
        y = "density"
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(
        legend.position = "none",
        plot.title = ggplot2::element_text(size = 12)
      )

  })

  patchwork::wrap_plots(
    plots,
    ncol = 2
  ) +
    patchwork::plot_annotation(
      title = "Method 2: Marginal density overlays",
      subtitle = sprintf(
        "Blue = %s  |  Red = %s",
        label_A,
        label_B
      ),
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(size = 14)
      )
    )
}
