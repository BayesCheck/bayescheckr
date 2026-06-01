#geweke.R

#inputs: the huge object outputted by simulate.R

compare_geweke <- function(...) {

  direct_draws <- sims[[i]][theta_tilde:y]

  gibbs_draws <- sims[[i]][draws]
}
