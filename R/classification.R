# ------------------------------------------------------------------------------
# PATCH: pluggable classifier for Method 1 (replaces .run_influential() in
# DistributionShift.R). Drop these functions in; remove the old xgboost-only
# .run_influential(). run_distributional_shift() gains a `classifier` arg.
#
# Classifier interface -- any classifier is just a 2-element list:
#   train(X, y)        -> a fitted model object, where y is 0/1 numeric
#   predict(fit, newX)  -> numeric vector of P(y = 1) for each row of newX
#
# X / newX are passed as plain matrices/data.frames of theta values (same
# columns as theta_A / theta_B). Nothing else about the classifier's internals
# matters to the rest of the pipeline.
# ------------------------------------------------------------------------------


# ---- default classifier: plain logistic regression (base R, no extra deps) --

.default_glm_classifier <- function() {
  list(
    train = function(X, y) {
      df <- data.frame(X, .y = y)
      stats::glm(.y ~ ., data = df, family = stats::binomial())
    },
    predict = function(fit, newX) {
      stats::predict(fit, newdata = as.data.frame(newX), type = "response")
    }
  )
}


# ---- xgboost classifier: opt-in, matches the old default behavior ----------
#' default_xgb_classifier
#'
#' Runs XG boost classifier instead of the default logistic regression.
#'
#' @param nrounds How many "cuts" XGBoost.
#' @return returns a list of two functions: train and test.
#' @export
default_xgb_classifier <- function(nrounds = 100) {
  list(
    train = function(X, y) {
      dtrain <- xgboost::xgb.DMatrix(as.matrix(X), label = y)
      xgboost::xgb.train(
        params = list(
          objective        = "binary:logistic",
          eval_metric      = "logloss",
          max_depth        = 4,
          eta              = .1,
          subsample        = .8,
          colsample_bytree = .8
        ),
        data    = dtrain,
        nrounds = nrounds,
        verbose = 0
      )
    },
    predict = function(fit, newX) {
      predict(fit, xgboost::xgb.DMatrix(as.matrix(newX)))
    }
  )
}


# ---- generic permutation feature importance, works for any classifier ------

.permutation_importance <- function(fit, predict_fn, Xtest, ytest,
                                    baseline_acc, n_perm = 99) {

  feature_cols <- colnames(Xtest)

  purrr::map_dfr(feature_cols, function(f) {
    perm_acc <- replicate(n_perm, {
      Xperm      <- Xtest
      Xperm[, f] <- sample(Xperm[, f])
      p_perm     <- predict_fn(fit, Xperm)
      mean(ifelse(p_perm > 0.5, 1, 0) == ytest, na.rm = TRUE)
    })

    data.frame(
      feature       = f,
      observed_drop = baseline_acc - mean(perm_acc),
      p_value       = (sum(perm_acc >= baseline_acc) + 1) / (n_perm + 1)
    )
  })
}


# ---- Method 1, now classifier-agnostic --------------------------------------

.run_influential <- function(theta_A,
                             theta_B,
                             varnames,
                             classifier  = NULL,   # <- pluggable
                             train_frac  = 0.8,
                             nrounds     = 100,     # only used by the default xgb classifier
                             n_perm      = 20,
                             alpha       = 0.05) {  # significance level for h0 decision

  if (is.null(classifier)) {
    classifier <- .default_glm_classifier()
  }
  stopifnot(is.function(classifier$train), is.function(classifier$predict))

  ##-----------------------
  ## Split A / Split B (unchanged)
  ##-----------------------

  idx_A <- sample(seq_len(nrow(theta_A)), floor(train_frac * nrow(theta_A)))
  A_train <- theta_A[idx_A, , drop = FALSE]
  A_test  <- theta_A[-idx_A, , drop = FALSE]

  idx_B <- sample(seq_len(nrow(theta_B)), floor(train_frac * nrow(theta_B)))
  B_train <- theta_B[idx_B, , drop = FALSE]
  B_test  <- theta_B[-idx_B, , drop = FALSE]

  Xtrain <- rbind(A_train, B_train)
  ytrain <- c(rep(0, nrow(A_train)), rep(1, nrow(B_train)))

  Xtest <- rbind(A_test, B_test)
  ytest <- c(rep(0, nrow(A_test)), rep(1, nrow(B_test)))
  n_te  <- length(ytest)

  ############################################################
  ## Train / predict via the pluggable interface
  ############################################################

  fit_AB <- classifier$train(Xtrain, ytrain)
  p      <- classifier$predict(fit_AB, Xtest)
  pred   <- ifelse(p > 0.5, 1, 0)

  accuracy <- mean(pred == ytest)

  ############################################################
  ## C2ST statistics
  ##
  ## Two-sided, to match tabulate_classifier_tests(): under H0 ("A and B
  ## come from the same distribution"), accuracy ~ Normal(0.5, 1/(4*n_te)).
  ## We test |accuracy - 0.5| against that null rather than only the
  ## one-sided "accuracy > chance" direction, since an accuracy well below
  ## chance is also evidence against H0 (e.g. a mis-specified or inverted
  ## classifier) and Trisha's version flags that case too.
  ############################################################

  c2st_se     <- sqrt(1 / (4 * n_te))
  c2st_z      <- (accuracy - 0.5) / c2st_se
  c2st_pvalue <- 2 * stats::pnorm(abs(c2st_z), lower.tail = FALSE)
  c2st_h0     <- ifelse(c2st_pvalue <= alpha, "Reject", "Fail to Reject")

  c2st <- data.frame(
    statistic = accuracy,
    n_test    = n_te,
    se_null   = c2st_se,
    z         = c2st_z,
    p_value   = c2st_pvalue,
    alpha     = alpha,
    h0        = c2st_h0
  )

  ############################################################
  ## Feature importance -- same classifier, any type
  ############################################################

  feature_importance <- .permutation_importance(
    fit          = fit_AB,
    predict_fn   = classifier$predict,
    Xtest        = Xtest,
    ytest        = ytest,
    baseline_acc = accuracy,
    n_perm       = n_perm
  )

  list(
    classifier_accuracy = accuracy,
    c2st                = c2st,
    feature_importance  = feature_importance
  )
}


# ---- run_distributional_shift(): pass `classifier` through -----------------
#
# Only the signature and the .run_influential() call change; everything else
# (Method 2, message()s, output list) stays the same.

#' run_distributional_shift
#'
#' Runs an input classifier to determine whether the two distributions can be determined as different by the machine or not.
#' Ideally you want 50% accuracy with a uniform feature importance.
#' In the case of RJMCMC or other sampling methods that output samples in varying dimensions you may have trouble with some classification tools. The solution usualy is to put all samples in a max dimensional matrix and just leave the empty dimentions blank for low dimention samples.
#' In the case of high dimensional samples (I.E. 15+ dimensions) sometimes your classifier will be unable to work properly, such as SVM. In this case you will need a different classification method. We recommend xgboost.
#'
#' @param direct_draws Direct draws from the joint distribution p(y, theta)
#' @param gibbs_draws Gibbs sampler draws from a simulated joint made up of successive conditional draws
#' @param classifier The algorithm used for classifying which proportion of the data is from which distribution. Should output a matrix of predictions.
#' @param train_frac The fraction of draws from your sampler that you want your classifier to train on. The remainder will be tested on for accuracy.
#' @param n_rounds Used specifically for certain types of classifiers that are iterative, like random forests or xgboost. Its the number of times it itterates.
#' @param n_perm Used to determine feature importance, does not effect the classifier.
#' @param alpha Used to bound the P value of the C2ST test on our classifier.
#' @return returns a list of the form: list(classifier_accuracy = accuracy,c2st = c2st, feature_importance = feature_importance).
#' @export
run_distributional_shift <- function(
    direct_draws,
    gibbs_draws,
    classifier  = NULL,
    train_frac  = 0.8,
    nrounds     = 100,
    n_perm      = 20,
    alpha       = 0.05
) {

  theta_A  <- cbind(direct_draws$theta, direct_draws$y)
  theta_B  <- cbind(gibbs_draws$theta, gibbs_draws$y)
  varnames <- colnames(theta_A)

  if (!identical(varnames, colnames(theta_B))) {
    stop("direct_draws$theta and gibbs_draws$theta must have the same column names.")
  }

  message("Training A-vs-B classifier ...")

  m1 <- .run_influential(
    theta_A,
    theta_B,
    varnames,
    classifier = classifier,
    train_frac = train_frac,
    nrounds    = nrounds,
    n_perm     = n_perm,
    alpha      = alpha
  )
}


# ==============================================================================
# USAGE EXAMPLES
# ==============================================================================

# 1) Default behavior -- now plain logistic regression (base R, no extra
#    dependency required):
#
#   result <- run_distributional_shift(direct_draws, gibbs_draws)


# 2) Opt back into xgboost (the old default):
#
#   result <- run_distributional_shift(
#     direct_draws, gibbs_draws,
#     classifier = .default_xgb_classifier(nrounds = 100)
#   )


# 3) Plug in a tidymodels/parsnip spec, e.g. a random forest via ranger
#    (matches the tidymodels style already used in Trisha's
#    tabulate_classifier_tests()):
#
#   rf_classifier <- list(
#     train = function(X, y) {
#       spec <- parsnip::rand_forest(trees = 500) |>
#         parsnip::set_engine("ranger") |>
#         parsnip::set_mode("classification")
#       df <- data.frame(X, .y = factor(y, levels = c(0, 1), labels = c("A", "B")))
#       parsnip::fit(spec, .y ~ ., data = df)
#     },
#     predict = function(fit, newX) {
#       predict(fit, new_data = as.data.frame(newX), type = "prob")$.pred_B
#     }
#   )
#
#   result <- run_distributional_shift(
#     direct_draws, gibbs_draws,
#     classifier = rf_classifier
#   )


# 4) Plug in glmnet (regularized logistic regression):
#
#   glmnet_classifier <- list(
#     train = function(X, y) {
#       glmnet::cv.glmnet(as.matrix(X), y, family = "binomial")
#     },
#     predict = function(fit, newX) {
#       as.numeric(predict(fit, as.matrix(newX), s = "lambda.min", type = "response"))
#     }
#   )
