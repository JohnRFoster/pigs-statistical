

##-------- parameters --------
# beta, vector, non time varying coefficients - process model
# beta_p, vector, non time varying coefficients - observation model
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
  ### priors

  # non time varying coefficients - process model
  for(i in 1:m_n){
    beta[i] ~ dnorm(0, tau = 1)
  }

  # non time varying coefficients - observation model
  for(i in 1:m_p){
    beta_p[i] ~ dnorm(0, tau = 1)
  }

  for(i in 1:n_method){
    p_unique[i] ~ dbeta(2, 2)
  }

  for(i in 1:n_property){
    eps_propertyR[i] ~ dnorm(0, tau = 1)
    eps_property_pR[i] ~ dnorm(0, tau = 1)
  }

  tau_dcar ~ dgamma(0.001, 0.001)
  sigma_property ~ dnorm(0, tau = 1)
  sigma_property_p ~ dnorm(0, tau = 1)
  sigma_st0 ~ dnorm(0, tau = 1) # scale parameter at time 1 for spatially uncorrelated temporal autocorrelation
  sigma_st ~ dnorm(0, tau = 1) # scale parameter at time > 1 for spatially uncorrelated temporal autocorrelation
  eta ~ dbeta(9, 1) # temporal autocorrelation parameter (AR1)

  for (i in 1:n_timestep) {
    for(j in 1:n_county){
      eps_stR[i, j] ~ dnorm(0, tau = 1)
    }
  }

  sigma_short ~ dnorm(0, tau = 1)


  for(i in 1:3){
    gamma[i] ~ dgamma(1, 1)
  }
  gamma[4] ~ dgamma(7.704547, 4.41925)       # snare
  gamma[5] ~ dgamma(3.613148, 3.507449)      # trap

  log_rho[1] ~ dnorm(-1.654013, 0.8092821)   # firearms
  log_rho[2] ~ dnorm(2.74274, 0.5976205)     # fixed wing
  log_rho[3] ~ dnorm(2.74274, 0.5976205)     # helicopter
  log_rho[4] ~ dnorm(0.2777449, 0.2255966)   # snare
  log_rho[5] ~ dnorm(0.6438721, 0.2225168)   # trap

  eps_property[1:n_property] <- eps_propertyR[1:n_property] * sigma_property

  for(i in 1:n_survey){

    log_potential_area[i] <- log(1 + p_unique[method[i]] * n_trap_m1[i]) +
      log_effort_per[i] +
      log_rho[method[i]] * shooting_ind[i] + # shooting and aerial only
      (log_pi + 2 * log_rho[method[i]] - log(gamma[method[i]] + exp(log_effort_per[i]))) * trap_snare_ind[i] # traps and snares only

    area_diff[i] <- log_potential_area[i] - log_survey_area_km2[i]
    log_pr_area_sampled[i] <- min(0, area_diff[i])

    logit(theta[i]) <- X_p[i, 1:m_p] %*% beta_p[1:m_p] +
      eps_property_pR[p_property_idx[i]] * sigma_property_p

    log_theta[i] <- log(theta[i])

    log_p[i] <- log_theta[i] + log_pr_area_sampled[i] +
      sum(log(1 - exp(log_theta[start[i]:end[i]])))*not_first_survey[i] # if > 1st survey

    # likelihood
    lambda[i] <- exp(log_lambda[survey_idx[i]] + log_p[i])
    y[i] ~ dpois(lambda[i])
  }


  # auto regressive prior on the first basis vector
  for(i in 1:m_short){
    z_shortR[i, 1:n_county] ~ dcar_normal(
      adj = adj[1:n_edges],
      weights = weights[1:n_edges],
      num = num[1:n_county],
      tau = tau_dcar,
      zero_mean = 1
    )
  }

  # loop through autocorrelation
  z_short[1, 1:n_county] <- z_shortR[1, 1:n_county]
  for(i in 2:m_short){
    z_short[i, 1:n_county] <- z_shortR[i - 1, 1:n_county] +
      z_shortR[i, 1:n_county] * sigma_short
  }

  # create the spatio-temporal adjustment
  eps_st[1, 1:n_county] <- eps_stR[1, 1:n_county] * sigma_st0
  for(t in 2:n_timestep){
    # at subsequent time steps, there is autocorrelation (eta) between time steps
    eps_st[t, 1:n_county] <- eta * eps_st[t - 1, 1:n_county] +
      eps_stR[t, 1:n_county] * sigma_st
  }

  # expected pig abundance
  for (i in 1:n_st) {
    log_lambda[i] <- X[i, 1:m_n] %*% beta[1:m_n] + # non time varying covariates
      X_short[i, 1:m_short] %*% z_short[1:m_short, county_idx[i]] + # B-spline vectors * basis vector coefficients
      eps_property[property[i]] + # property adjustment
      eps_st[timestep[i], county_idx[i]] + # spatio-temporal adjustment
      log_area_km2[i] # area offset
  }

})
