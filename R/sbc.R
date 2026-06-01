library(future)
library(future.apply)

# =============================================================================
# CONFIG
# =============================================================================

N_SIM   <- 200000
L_DRAWS <- 1000

PRIOR <- list(
  name        = "Normal(5, 16)",
  hyperparams = c(mean = 5, sd = 16),
  draw        = function(hp) rnorm(1, hp["mean"], hp["sd"])
)

LIKELIHOOD <- list(
  name     = "Normal obs, n=30",
  simulate = function(n_obs, theta) rnorm(n_obs, mean = theta, sd = 1)
)

POSTERIOR <- list(
  name = "Normal approximation",
  draw = function(y, n_obs, hp, L) {
    post_mean <- mean(y) - 0.1
    post_sd   <- 1 / sqrt(n_obs)
    rnorm(L, mean = post_mean, sd = post_sd)
  }
)

# =============================================================================
# RANKS
# =============================================================================

compute_ranks <- function(sims) {
  param_names <- names(sims$sims[[1]]$theta_tilde)
  n_params    <- length(param_names)

  ranks <- matrix(NA,
                  nrow     = sims$n_sims,
                  ncol     = n_params,
                  dimnames = list(NULL, param_names))

  for (i in seq_len(sims$n_sims)) {
    theta_tilde <- sims$sims[[i]]$theta_tilde
    draws       <- sims$sims[[i]]$draws
    for (p in param_names) {
      ranks[i, p] <- sum(draws[, p] < theta_tilde[p])
    }
  }

  structure(
    list(ranks = ranks, n_draws = sims$n_draws, n_sims = sims$n_sims),
    class = "bayescheckr_ranks"
  )
}

# =============================================================================
# DIAGNOSTICS
# =============================================================================

plot_sbc_diagnostics <- function(ranks, L, n_sim,
                                 prior_name, likelihood_name, posterior_name) {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mfrow = c(1, 3), mar = c(4.5, 4.5, 3.5, 1.5),
      oma = c(0, 0, 4, 0), bg = "#F7F7F2")

  n_bins   <- min(30, L + 1)
  breaks   <- seq(0, L, length.out = n_bins + 1)
  expected <- n_sim / n_bins
  ci_lo    <- qbinom(0.005, n_sim, 1 / n_bins)
  ci_hi    <- qbinom(0.995, n_sim, 1 / n_bins)

  h        <- hist(ranks, breaks = breaks, plot = FALSE)
  bar_cols <- ifelse(h$counts < ci_lo | h$counts > ci_hi, "#E05C5C", "#4A90D9")

  hist(ranks, breaks = breaks, col = bar_cols, border = "white",
       main = "Rank histogram", xlab = paste0("Rank  (0-", L, ")"),
       ylab = "Count", xlim = c(0, L), las = 1)
  abline(h = expected, lty = 2, lwd = 1.8, col = "grey30")
  rect(0, ci_lo, L, ci_hi, col = adjustcolor("grey60", 0.15), border = NA)

  u    <- (ranks + runif(length(ranks))) / (L + 1)
  ks_p <- ks.test(u, "punif")$p.value
  xs   <- seq(0, 1, length.out = 500)
  eps  <- 1.36 / sqrt(n_sim)

  plot(xs, ecdf(u)(xs), type = "l", lwd = 2, col = "#4A90D9",
       main = "ECDF vs Uniform", xlab = "(Rank + U) / (L + 1)",
       ylab = "Cumulative probability", las = 1)
  abline(0, 1, lty = 2, lwd = 1.5, col = "grey40")
  lines(xs, pmin(xs + eps, 1), lty = 3, col = "grey60")
  lines(xs, pmax(xs - eps, 0), lty = 3, col = "grey60")
  mtext(sprintf("KS p-value = %.3f", ks_p), side = 3, cex = 0.75, col = "grey30")

  z_scores <- (h$counts - expected) / sqrt(expected * (1 - 1 / n_bins))
  bar_col2 <- ifelse(abs(z_scores) > 2.58, "#E05C5C", "#7BC8A4")

  barplot(z_scores, col = bar_col2, border = "white", main = "Bin z-scores",
          xlab = "Bin", ylab = "z-score",
          ylim = c(min(-3.5, min(z_scores) - 0.5), max(3.5, max(z_scores) + 0.5)),
          las = 1)
  abline(h = c(-2.58, 0, 2.58), lty = c(2, 1, 2), lwd = c(1.5, 1, 1.5),
         col = c("grey40", "grey30", "grey40"))

  mtext(paste0("SBC  |  Prior: ", prior_name, "   Likelihood: ", likelihood_name,
               "\nPosterior: ", posterior_name,
               "   (", n_sim, " sims, L = ", L, ")"),
        outer = TRUE, cex = 0.85, font = 2, col = "grey20")
}

summarise_sbc <- function(ranks, L, n_sim) {
  u  <- (ranks + runif(length(ranks))) / (L + 1)
  ks <- ks.test(u, "punif")

  cat(sprintf("\nSimulations: %d  |  L: %d  |  Rank range: %d-%d\n",
              n_sim, L, min(ranks), max(ranks)))
  cat(sprintf("Mean rank: %.2f (ideal %.2f)  |  SD: %.2f (ideal %.2f)\n",
              mean(ranks), L / 2, sd(ranks), sqrt(L * (L + 2) / 12)))
  cat(sprintf("KS p-value: %.4f  %s\n", ks$p.value,
              ifelse(ks$p.value < 0.05, "<- possible miscalibration",
                     "consistent with uniform")))
}

# =============================================================================
# RUN
# =============================================================================

plan(multisession)

sims   <- run_simulations(prior_sampler      = PRIOR$draw,
                          likelihood_sampler  = LIKELIHOOD$simulate,
                          posterior_sampler   = POSTERIOR$draw,
                          n_sims             = N_SIM,
                          n_obs              = 30,
                          n_draws            = L_DRAWS,
                          prior_hyper_params = PRIOR$hyperparams)

ranked <- compute_ranks(sims)

ranks_vec <- ranked$ranks[, 1]

summarise_sbc(ranks_vec, L_DRAWS, N_SIM)
plot_sbc_diagnostics(ranks_vec, L_DRAWS, N_SIM,
                     PRIOR$name, LIKELIHOOD$name, POSTERIOR$name)

plan(sequential)
