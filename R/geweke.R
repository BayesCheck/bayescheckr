#geweke.R

#inputs: the huge object outputted by simulate.R

compare_geweke <- function(...) {

  direct_draws <- sims[[i]][theta_tilde:y]

  gibbs_draws <- sims[[i]][draws]
}

#check the marginals

par(mfrow = c(3, 3))
probs <- seq(0.05, 0.95, by = 0.05)
pvals <- numeric(n + 2)
labels <- c(paste0("beta", 1:p), "sigma", paste0("y", 1:n))
for(i in 1:(n + p + 1)){
  quantiles_direct <- quantile(direct_draws[, i], probs)
  quantiles_gibbs <- quantile(gibbs_draws[, i], probs)
  plot(quantiles_direct, quantiles_gibbs, col = "red", pch = 19, cex = 0.5,
       main = labels[i], xlab = "Direct draw", ylab = "Gibbs draw")
  abline(a = 0, b = 1)
  pvals[i] <- ks.test(direct_draws[, i], gibbs_draws[, i])$p.value
}

#Rafael's test functions

par(mfrow = c(1, 1))
probs <- seq(0.05, 0.95, by = 0.05)

g_direct <- direct_draws[, 1] * direct_draws[, 2] * direct_draws[, 3]
g_gibbs <- gibbs_draws[, 1] * gibbs_draws[, 2] * gibbs_draws[, 3]

quantiles_direct <- quantile(g_direct, probs)
quantiles_gibbs <- quantile(g_gibbs, probs)
plot(quantiles_direct, quantiles_gibbs, col = "red", pch = 19, cex = 0.5,
     main = "beta_1 * beta_2 * sigma^2")
abline(a = 0, b = 1)

#Histogram of the product:

hist(g_direct[abs(g_direct) < 50], breaks = "Scott")
hist(g_gibbs[abs(g_gibbs) < 50], breaks = "Scott")


#Let me also use different test functions. First, let's look at *second moments*:


par(mfrow = c(1, 3))
probs <- seq(0.05, 0.95, by = 0.05)
pvals_second <- numeric(p + 1)
labels_second <- c(paste("beta", 1:p, sep = "_"), "sigma2")

for(i in 1:(p + 1)){
  # compute second moments by squaring
  second_direct <- direct_draws[, i]^2
  second_gibbs <- gibbs_draws[, i]^2

  quantiles_direct <- quantile(second_direct, probs)
  quantiles_gibbs <- quantile(second_gibbs, probs)

  plot(quantiles_direct, quantiles_gibbs,
       col = "blue", pch = 19, cex = 0.5,
       main = paste(labels_second[i], "squared"),
       xlab = "direct", ylab = "gibbs")
  abline(a = 0, b = 1)
}

#Finally let's look at the test function $h(\theta, y) = |\beta_1 - \beta_2|$:


par(mfrow = c(1, 1))
probs <- seq(0.05, 0.95, by = 0.05)
h_direct <- abs(direct_draws[, 1] - direct_draws[, 2])
h_gibbs <- abs(gibbs_draws[, 1] - gibbs_draws[, 2])
quantiles_direct <- quantile(h_direct, probs)
quantiles_gibbs <- quantile(h_gibbs, probs)
plot(quantiles_direct, quantiles_gibbs, col = "red", pch = 19, cex = 0.5,
     main = "|beta_1 - beta_2|",
     xlab = "direct", ylab = "gibbs")
abline(a = 0, b = 1)
