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
#   tabulate_shift_divergences()   Method 3 — numeric divergence table
#   plot_shift_importance()        Method 1 — permutation importance bar chart
#   plot_shift_prototypes()        Method 2 — neighbourhood density bar chart
#   plot_shift_densities()         Method 3 — marginal density overlays
#
# Dependencies: e1071, cluster, ggplot2, tidyr, patchwork


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


.neighbourhood_density <- function(prototype, samples, radius) {
  dists <- sqrt(rowSums(sweep(as.matrix(samples), 2,
                              as.numeric(prototype), "-")^2))
  mean(dists < radius)
}


.auto_radius <- function(theta_matrix, max_rows = 500) {
  idx <- sample(nrow(theta_matrix), min(max_rows, nrow(theta_matrix)))
  2* median(dist(theta_matrix[idx, , drop = FALSE]))
}


# ------------------------------------------------------------------------------
# Method 1 — influential-example explanations
# ------------------------------------------------------------------------------


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

  p <- predict(
    fit_AB,
    xgboost::xgb.DMatrix(as.matrix(Xtest))
  )

  pred <- ifelse(p>.5,1,0)

  accuracy <- mean(pred==ytest)



  ############################################################
  ## Importance
  ############################################################

  importance_AB <-
    xgboost::xgb.importance(
      feature_names=varnames,
      model=fit_AB
    )
  importance_AB <- importance_AB[, c("Feature", "Gain")]

  colnames(importance_AB) <- c(
    "variable",
    "importance"
  )


  ############################################################
  ## Reverse labels
  ############################################################

  ytrain_rev <- 1-ytrain

  fit_BA <- xgboost::xgb.train(
    params=list(
      objective="binary:logistic",
      eval_metric="logloss",
      max_depth=4,
      eta=.1,
      subsample=.8,
      colsample_bytree=.8
    ),
    data=xgboost::xgb.DMatrix(
      as.matrix(Xtrain),
      label=ytrain_rev
    ),
    nrounds=nrounds,
    verbose=0
  )

  importance_BA <-
    xgboost::xgb.importance(
      feature_names=varnames,
      model=fit_BA
    )
  importance_BA <- importance_BA[, c("Feature", "Gain")]

  colnames(importance_BA) <- c(
    "variable",
    "importance"
  )


  list(

    classifier_accuracy = accuracy,

    importance = importance_AB,

    importance_A_to_B = importance_AB,

    importance_B_to_A = importance_BA

  )
}


# ------------------------------------------------------------------------------
# Method 2 — prototype-neighbourhood explanations
# ------------------------------------------------------------------------------


.run_prototypes <- function(theta_A, theta_B, varnames, k, radius) {

  pam_fit <- cluster::pam(theta_A, k = k, metric = "euclidean")
  protos  <- as.data.frame(pam_fit$medoids)

  rows <- lapply(seq_len(k), function(i) {
    proto  <- protos[i, ]
    dens_A <- .neighbourhood_density(proto, theta_A, radius)
    dens_B <- .neighbourhood_density(proto, theta_B, radius)
    data.frame(
      prototype      = paste0("P", i),
      density_A      = dens_A,
      density_B      = dens_B,
      ratio_B_over_A = dens_B / max(dens_A, 1e-9),
      verdict        = ifelse(
        dens_B < dens_A * 0.7, "Under-represented in B",
        ifelse(dens_B > dens_A * 1.3, "Over-represented in B", "Similar")
      ),
      stringsAsFactors = FALSE
    )
  })

  list(
    prototypes  = protos,
    comparison  = do.call(rbind, rows),
    radius_used = radius,
    k_used      = k
  )
}


# ------------------------------------------------------------------------------
# Method 3 — divergence metrics
# ------------------------------------------------------------------------------


.run_divergences <- function(theta_A, theta_B, varnames) {

  rows <- lapply(varnames, function(var) {
    x      <- theta_A[, var]
    y      <- theta_B[, var]
    ks_res <- ks.test(x, y)   # call once; reuse both statistic and p.value
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
    nrounds = 100,
    k_prototypes = 4L,
    radius = NULL
) {

  theta_A  <- direct_draws$theta
  theta_B  <- gibbs_draws$theta
  varnames <- colnames(theta_A)

  if (!identical(varnames, colnames(theta_B))) {
    stop("direct_draws$theta and gibbs_draws$theta must have the same column names.")
  }

  if (is.null(radius)) {
    radius <- .auto_radius(theta_A)
  }

  k <- min(k_prototypes, nrow(theta_A) - 1L)

  message("Method 1: training A-vs-B classifier ...")
  m1 <- .run_influential(
    theta_A,
    theta_B,
    varnames,
    train_frac = train_frac,
    nrounds = nrounds
  )

  message("Method 2: computing prototype neighbourhoods ...")
  m2 <- .run_prototypes(theta_A, theta_B, varnames, k = k, radius = radius)

  message("Method 3: computing divergence metrics ...")
  m3 <- .run_divergences(theta_A, theta_B, varnames)

  synthesis <- merge(
    m1$importance_A_to_B,
    m3,
    by = "variable",
    all.x = TRUE
  )
  synthesis <- synthesis[order(-synthesis$importance), ]
  rownames(synthesis) <- NULL

  list(
    influential = m1,
    prototypes  = m2,
    divergences = m3,
    synthesis   = synthesis,
    varnames    = varnames,
    n_A         = nrow(theta_A),
    n_B         = nrow(theta_B)
  )
}


# ------------------------------------------------------------------------------
# Tabulate divergences — Method 3 numeric output
# ------------------------------------------------------------------------------

#' @export
tabulate_shift_divergences <- function(shift_result) {
  shift_result$divergences
}


# ------------------------------------------------------------------------------
# Plot functions — each returns a ggplot object
# ------------------------------------------------------------------------------

#' @export
plot_shift_importance <- function(shift_result,
                                  label_A = "Direct (MC)",
                                  label_B = "Gibbs (SC)") {

  imp <- shift_result$influential$importance_A_to_B
  acc <- shift_result$influential$classifier_accuracy

  # coerce importance > 0 to character so scale_fill_manual keys are unambiguous
  imp$positive <- ifelse(imp$importance > 0, "positive", "negative")

  ggplot2::ggplot(imp,
                  ggplot2::aes(x = reorder(variable, importance),
                               y = importance,
                               fill = positive)) +
    ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
    ggplot2::scale_fill_manual(
      values = c(positive = "#D85A30", negative = "#AAAAAA")
    ) +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title    = "Method 1: Which parameters drive the distributional difference?",
      subtitle = sprintf(
        "Permutation importance from %s-vs-%s classifier  |  accuracy = %.3f",
        label_A, label_B, acc),
      x = NULL,
      y = "Importance (drop in classifier accuracy)"
    ) +
    ggplot2::theme_minimal(base_size = 13)
}


#' @export
plot_shift_prototypes <- function(shift_result,
                                  label_A = "Direct (MC)",
                                  label_B = "Gibbs (SC)") {

  comp <- shift_result$prototypes$comparison[, c("prototype",
                                                 "density_A",
                                                 "density_B")]

  long <- tidyr::pivot_longer(comp,
                              cols      = c("density_A", "density_B"),
                              names_to  = "population",
                              values_to = "density"
  )
  long$population <- ifelse(long$population == "density_A", label_A, label_B)

  ggplot2::ggplot(long,
                  ggplot2::aes(x = prototype, y = density, fill = population)) +
    ggplot2::geom_col(position = "dodge", width = 0.6) +
    ggplot2::scale_fill_manual(
      values = setNames(c("#378ADD", "#D85A30"), c(label_A, label_B))
    ) +
    ggplot2::labs(
      title    = "Method 2: Prototype-neighbourhood density comparison",
      subtitle = sprintf(
        "Fraction of draws within radius %.3g of each A-prototype",
        shift_result$prototypes$radius_used),
      x    = "Prototype",
      y    = "Local density",
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 13)
}


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
      population = rep(c(label_A, label_B),
                       c(nrow(theta_A), nrow(theta_B)))
    )
    ggplot2::ggplot(df,
                    ggplot2::aes(x = value, fill = population, colour = population)) +
      ggplot2::geom_density(alpha = 0.3, linewidth = 0.8) +
      ggplot2::scale_fill_manual(
        values = setNames(c("#378ADD", "#D85A30"), c(label_A, label_B))
      ) +
      ggplot2::scale_colour_manual(
        values = setNames(c("#185FA5", "#993C1D"), c(label_A, label_B))
      ) +
      ggplot2::labs(title = var, x = NULL, y = "density") +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(legend.position = "none",
                     plot.title      = ggplot2::element_text(size = 12))
  })

  patchwork::wrap_plots(plots, ncol = 2) +
    patchwork::plot_annotation(
      title    = "Method 3: Marginal density overlays",
      subtitle = sprintf("Blue = %s  |  Red = %s", label_A, label_B),
      theme    = ggplot2::theme(plot.title = ggplot2::element_text(size = 14))
    )
}
