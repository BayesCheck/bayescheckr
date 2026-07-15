#' Fit a classifier for C2ST visualization
#'
#' Mirrors the train/test split and classifier-fitting steps used internally
#' by \code{run_distributional_shift()}, but returns the fitted classifier
#' and held-out predictions so they can be used for visualization (e.g.
#' predicted-probability histograms, decision boundary plots, or as input
#' to a dimensionality-reduction technique like t-SNE/UMAP/PaCMAP).
#'
#' @param theta_A Matrix/data.frame of theta values from the direct draws.
#' @param theta_B Matrix/data.frame of theta values from the Gibbs draws.
#'   Must have the same column names as \code{theta_A}.
#' @param classifier A pluggable classifier list with \code{train}/\code{predict}
#'   elements (see \code{run_distributional_shift()}). Defaults to
#'   \code{default_glm_classifier()}.
#' @param train_frac Fraction of each of theta_A/theta_B used for training.
#'   The remainder is held out as the test set. Default 0.8.
#' @param seed Optional integer seed for reproducibility. If supplied, the
#'   function saves and restores the global RNG state on exit, so calling
#'   this with a seed does not affect subsequent random draws in the session.
#'
#' @return A list with:
#' \describe{
#'   \item{fit}{The fitted classifier object.}
#'   \item{Xtest}{The held-out feature data.frame (theta values).}
#'   \item{ytest}{True labels for the test set (0 = A/direct, 1 = B/gibbs).}
#'   \item{p_test}{Classifier's predicted probability of class 1 (gibbs) per test row.}
#' }
#'
#' @export
fit_c2st_for_viz <- function(theta_A, theta_B, classifier = NULL,
                             train_frac = 0.8, seed = NULL) {

  varnames <- colnames(theta_A)
  if (!identical(varnames, colnames(theta_B))) {
    stop("theta_A and theta_B must have the same column names.")
  }
  if (train_frac <= 0 || train_frac >= 1) {
    stop("train_frac must be strictly between 0 and 1.")
  }
  if (nrow(theta_A) < 2 || nrow(theta_B) < 2) {
    stop("theta_A and theta_B must each have at least 2 rows.")
  }

  if (is.null(classifier)) classifier <- default_glm_classifier()
  stopifnot(is.function(classifier$train), is.function(classifier$predict))

  # seed handling: reproducible if requested, but doesn't leak into the
  # caller's subsequent random draws
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    } else {
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

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

  fit    <- classifier$train(Xtrain, ytrain)
  p_test <- classifier$predict(fit, Xtest)

  list(
    fit    = fit,
    Xtest  = as.data.frame(Xtest),
    ytest  = ytest,
    p_test = p_test
  )
}
