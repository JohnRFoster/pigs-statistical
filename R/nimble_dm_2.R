modelCode <- nimbleCode({

  # priors
  for(i in 1:n_method){
    log_rho[i] ~ dnorm(0, tau = 0.1)
  }

  for(i in 1:3){
    p_mu[i] ~ dnorm(0, tau = 1)
    logit(p_unique[i]) <- p_mu[i]
  }

  for(i in 1:2){
    log_gamma[i] ~ dnorm(0, tau = 0.1)
  }

  # non time varying coefficients - observation model
  for(m in 1:n_method){
    beta1[m] ~ dnorm(0, tau = 1)
    for(i in 1:m_p){
      beta_p[m, i] ~ dnorm(0, tau = 1)
    }
  }

  # estimate apparent survival
  phi_mu ~ dbeta(phi_mu_a, phi_mu_b)
  psi_phi ~ dgamma(1, 0.001)
  a_phi <- phi_mu * psi_phi
  b_phi <- (1 - phi_mu) * psi_phi

  log_mean_ls ~ dnorm(2, tau = 1)  # mean litter size
  log(mean_ls) <- log_mean_ls

  ## convert to expected number of pigs per primary period
  log_zeta_mu <- log(pp_len) + log_mean_ls - log(365)
  log(zeta) <- log_zeta_mu
  for(i in 1:n_ls){
    J[i] ~ dpois(mean_ls)
  }

  for(i in 1:n_survey){

    log_potential_area[i] <- calc_log_potential_area(
      log_rho = log_rho[1:n_method],
      log_gamma = log_gamma[1:2],
      p_unique = p_unique[1:3],
      log_effort_per = log_effort_per[i],
      effort_per = effort_per[i],
      n_trap_m1 = n_trap_m1[i],
      log_pi = log_pi,
      method = method[i]
    )

    # probability of capture, given that an individual is in the surveyed area
    log_theta[i] <- log(
      ilogit(beta1[method[i]] + inprod(X_p[i, 1:m_p], beta_p[method[i], 1:m_p]))
    ) +
      min(0, log_potential_area[i] - log_survey_area_km2[i])

    # likelihood
    y[i] ~ dpois(p[i] * (N[p_property_idx[i], p_pp_idx[i]] - y_sum[i]))

  }

  # the probability an individual is captured on the first survey
  for(i in 1:n_first_survey){
    log(p[first_survey[i]]) <- log_theta[first_survey[i]]
  }

  # the probability an individual is captured after the first survey
  for(i in 1:n_not_first_survey){
    log(p[not_first_survey[i]]) <- log_theta[start[not_first_survey[i]]] +
      sum(log(1 - exp(log_theta[start[not_first_survey[i]]:end[not_first_survey[i]]])))
  }

  if(spatial){
    sigma_car ~ dunif(0, 100)   # prior for variance components based on Gelman (2006)
    tau_car <- 1 / sigma_car^2
    s_car[1:n_non_islands] ~ dcar_normal(
      adj = adj[1:n_edges],
      weights = weights[1:n_edges],
      num = num[1:n_non_islands],
      tau = tau_car,
      # tau = 1,
      zero_mean = 1
    )

    for(i in 1:n_non_islands){
      s[non_islands[i]] <- s_car[i]
    }
    for(i in 1:n_islands){
      s_island[i] ~ dnorm(0, tau = 1)
      s[islands[i]] <- s_island[i]
    }
  } else {
    s[1:n_county] <- 0
  }

  for(i in 1:n_property){

    # eps_property_pR[i] ~ dnorm(0, tau = 1) # property effect in observation model

    log_lambda_1[i] ~ dunif(0, 10)
    log(N[i, PPNum[i, 1]]) <- log_lambda_1[i]

    # population growth across time steps
    for(j in 2:n_pp_prop[i]){ # loop through every PP, including missing ones

      lambda[i, all_pp[i, j]] <- N[i, all_pp[i, j-1]] * zeta / 2 +
        N[i, all_pp[i, j-1]] * phi[pH[i, j-1]]
      N[i, all_pp[i, j]] ~ dpois(lambda[i, all_pp[i, j]])
      phi[pH[i, j-1]] ~ dbeta(a_phi, b_phi)

    }

  }


  # for easier monitoring of abundance - long format
  for(i in 1:n_units){
    xn[i] <- N[property_x[i], pp_x[i]]

    # omega[i] ~ dunif(0, 1)
    # z[i] ~ dbern(omega[i])
  }

  # for(i in 1:n_county_units){ # county unit is a county x PP combination
  #
  #   ### county dispersion ###
  #   # county level abundance
  #   # a priori we know it cannot me less than total abundance across properties
  #   if(likelihood_nb){
  #     M[i] ~ T(dnegbin(pM[i], size[county[i]]), N_mu[i], 0)
  #
  #     # re-parametrization with mean abundance
  #     pM[i] <- size[county[i]] / (M_mu[i] + size[county[i]])
  #   } else {
  #     M[i] ~ dpois(M_mu[i])
  #   }
  #
  #   # convert abundance to density (across properties), scale to county abundance
  #   M_mu[i] <- N_mu[i] / sum_prop_area[i] * county_area[i]
  #
  #   # all pigs in county_timestep i
  #   N_mu[i] <- sum_properties(property = M_lookup[i, 1:n_prop[i]],
  #                             z = xn[1:n_units])
  #
  # }

})
