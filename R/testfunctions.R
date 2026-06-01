#testfunctions

quants <- derived_quantities(
  product = mu * sigmasq,
  log_like = sum(dnorm(y, mean = mu, sd = sqrt(sigmasq), log = TRUE)),
)


#-----Geweke Graphs-----


#------SBC Graphs-------

