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
