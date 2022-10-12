

##-------- parameters --------
# eps_property, vector, abundance ??
# sigma_property, scalar, ??
# sigma_short, scalar, ??
# z_shortR, matrix, B-spline spatial coefficient

# TODO explicit indexing of all vectors and matrices - put params in for loop
# TODO check that all log transformed variables are correct
# TODO add uninformed priors for capture/effect coefficents
# TODO finish list at top explaining the function and dimenstion of each param


modelName <- "saturating"

modelCode <- nimbleCode({
  # priors
  beta ~ dnorm(0, tau = 1)
  beta_p ~ dnorm(0, tau = 1)
  p_unique ~ dbeta(2, 2)
  sigma_property ~ dnorm(0, tau = 1)
  eps_propertyR ~ dnorm(0, tau = 1)
  sigma_property_p ~ dnorm(0, tau = 1)
  eps_property_pR ~ dnorm(0, tau = 1)
  sigma_st ~ dnorm(0, tau = 1)
  sigma_st0 ~ dnorm(0, tau = 1)
  eta ~ dbeta(9, 1)

  for (i in 1:n_timestep) {
    eps_stR[i, 1:n_county] ~ dnorm(0, tau = 1)
  }

  sigma_short ~ dnorm(0, tau = 1)

  for(i in 1:m_short){
    z_shortR[i, 1:n_county] ~ dcar_normal(
      adj = adj[1:n_edges],
      weights = weights[1:n_edges],
      num = num[1:n_county],
      tau = tau_dcar,
      zero_mean = 1
    )
  }


  gamma[1] ~ dgamma(7.704547, 4.41925)       # snare
  gamma[2] ~ dgamma(3.613148, 3.507449)      # trap
  log_rho[1] ~ dnorm(-1.654013, 0.8092821)   # firearms
  log_rho[2] ~ dnorm(2.74274, 0.5976205)     # fixed wing
  log_rho[3] ~ dnorm(2.74274, 0.5976205)     # helicopter
  log_rho[4] ~ dnorm(0.2777449, 0.2255966)   # snare
  log_rho[5] ~ dnorm(0.6438721, 0.2225168)   # trap

  # likelihood
  for(i in 1:n_survey){
    lambda[i] <- exp(log_lambda[survey_idx[i]] + log_p)
    y ~ dpois(lambda[i])
  }

  eps_property <- eps_propertyR * sigma_property

  for(i in 1:n_survey){

    log_potential_area <- log(1 + p_unique[method[i]] * n_trap_m1[i]) +
      log_effort_per[i] +
      log_rho[method[i]] * trap_snare_indicator[i] + # shooting and aerial only
      (log_pi + 2 * log_rho[method[i]] - log(exp(log_gamma[method[i] - 3]) + exp(log_effort_per[i]))) * trap_snare_indicator[i] # traps and snares only

    log_pr_area_sampled[i] <- min(0, log_potential_area - log_survey_area_km2[i])

    log_p[i] <- log_theta[i] +
      sum(log(1 - exp(log_theta[start[i]:end[i]])))*not_first_survey[i] # if > 1st survey
  }

  log_theta <- log(ilogit(X_p[1:n_survey, 1:m_p] * beta_p +
                            eps_property_pR[p_property_idx[1:n_survey]] * sigma_property_p)) +
    log_pr_area_sampled

  z_short[1, 1:n_county] <- z_shortR[1, 1:n_county]
  for(i in 2:m_short){
    z_short[i, 1:n_county] <- z_shortR[i - 1, 1:n_county] + z_shortR[i, 1:n_county] * sigma_short
  }

  eps_st[1] <- eps_stR[1] * sigma_st0
  for(t in 2:n_timestep){
    eps_st[t] <- eta * eps_st[t - 1] + eps_stR[t] * sigma_st
  }

  for (i in 1:n_st) {
    log_lambda[i] <- X[i, 1:m_n] * beta[1:m_n] +
      X_short[i, 1:m_short] * z_short[1:m_short, county_idx[i]] +
      eps_property[property[i]] +
      eps_st[timestep[i]][county_idx[i]] +
      log_area_km2[i]
  }

})
