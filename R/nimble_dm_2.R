modelCode <- nimbleCode({

  # priors
  for(i in 1:n_method){
    log_rho[i] ~ dnorm(0, tau = 0.1)
    p_mu[i] ~ dnorm(0, tau = 1)
    logit(p_unique[i]) <- p_mu[i]
  }

  for(i in 1:2){
    log_gamma[i] ~ dnorm(0, tau = 0.1)
  }

  # non time varying coefficients - observation model
  for(i in 1:m_p){
    beta_p[i] ~ dnorm(0, tau = 1)
  }

  # estimate apparent survival
  logit_mean_phi ~ dnorm(phi_prior_mean, tau = phi_prior_tau)
  sigma_phi ~ dgamma(0.01, 0.01)
  # log(sigma_phi) <- log_sigma_phi
  tau_phi <- 1/sigma_phi^2
  for(i in 1:n_property){
    for(t in 1:(n_pp_prop[i]-1)){
      # alpha_phi[pH[i, t]] ~ dnorm(0, tau = tau_phi)
      logit_phi[pH[i, t]] ~ dnorm(logit_mean_phi, tau = tau_phi)
      logit(phi[pH[i, t]]) <- logit_phi[pH[i, t]]
    }
  }

  # estimate per capita recruitment
  # mean_lpy <- exp(log_mean_lpy)  # mean litters per year
  # log_mean_lpy ~ dnorm(0, tau = 1)

  # mean_ls <- exp(log_mean_ls)  # mean litter size
  mean_ls ~ dgamma(1, 0.1)

  # for(i in 1:n_lpy){
  #   K[i] ~ dexp(mean_lpy)
  # }

  for(i in 1:n_ls){
    J[i] ~ dpois(mean_ls)
  }

  # convert to expected number of pigs per primary period
  zeta_mu <- 1 / 365 * pp_len * mean_ls
  for(i in 1:n_pp){
    zeta[i] <- zeta_mu
  }

  # property dispersion
  if(likelihood_nb){
    for(i in 1:n_county){
      size[i] ~ dunif(0, 20)
    }
  }

  # for(i in 1:n_trap_snare){
  #   log_potential_area[ts_idx[i]] <- log_pi +
  #     (2 * (log_rho[method[ts_idx[i]]] + log_effort_per[ts_idx[i]] -
  #             log(exp(log_gamma[method[ts_idx[i]]-3]) + effort_per[ts_idx[i]]))) +
  #     log(1 + (p_unique[method[ts_idx[i]]] * n_trap_m1[ts_idx[i]]))
  # }
  #
  # for(i in 1:n_shooting){
  #   log_potential_area[shooting_idx[i]] <- log_rho[method[shooting_idx[i]]] +
  #     log_effort_per[shooting_idx[i]] -
  #     # log(exp(log_gamma[method[shooting_idx[i]]]) + effort_per[shooting_idx[i]]) +
  #     log(1 + (p_unique[method[shooting_idx[i]]] * n_trap_m1[shooting_idx[i]]))
  # }

  for(i in 1:n_survey){

    log_potential_area[i] <- calc_log_potential_area(
      log_rho = log_rho[1:n_method],
      log_gamma = log_gamma[1:2],
      p_unique = p_unique[1:n_method],
      log_effort_per = log_effort_per[i],
      effort_per = effort_per[i],
      n_trap_m1 = n_trap_m1[i],
      log_pi = log_pi,
      method = method[i]
    )

    # log_pr_area_sampled[i] <- min(0, log_potential_area[i] - log_survey_area_km2[i])

    # probability of capture, given that an individual is in the surveyed area
    log_theta[i] <- log(ilogit(inprod(X_p[i, 1:m_p], beta_p[1:m_p]))) +
      min(0, log_potential_area[i] - log_survey_area_km2[i])

    # data model
    if(likelihood_binom){
      y[i] ~ dbinom(p[i], z[i])
    }

    if(likelihood_nb){
      y[i] ~ dnegbin(py[i], size[p_county_idx[i]])
      py[i] <- size[p_county_idx[i]] / (y_mu[i] + size[p_county_idx[i]])
      y_mu[i] <- p[i] * z[i]
    }

    if(likelihood_poisson){
      y[i] ~ dpois(p[i] * N[p_property_idx[i], p_pp_idx[i]])
    }

    # z[i] <- N[p_property_idx[i], p_pp_idx[i]] - y_sum[i]

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

    N[i, PPNum[i, 1]] ~ dpois(lambda_1[i])
    # log(lambda_1[i]) <- log_lambda_1[i]
    lambda_1[i] ~ dgamma(1, 0.001)

    for(t in 2:n_timesteps[i]){ # loop through sampled PP only
      N[i, PPNum[i, t]] ~ dpois(dm[i, PPNum[i, t]])
    }

    # population growth across time steps
    dm[i, all_pp[i, 1]] <- N[i, PPNum[i, 1]]
    for(j in 2:n_pp_prop[i]){ # loop through every PP, including missing ones

      Z[i, j-1] <- dm[i, all_pp[i, j-1]] - rem[i, j-1]

      if(demographic_stochasticity){
        S[i, j-1] ~ dbinom(phi[pH[i, j-1]], Z[i, j-1])
        R[i, j-1] ~ dpois(zeta[all_pp[i, j-1]] * Z[i, j-1] / 2)
      } else {
        S[i, j-1] <- phi[pH[i, j-1]] * Z[i, j-1]
        R[i, j-1] <- zeta[all_pp[i, j-1]] * Z[i, j-1] / 2
      }
      dm[i, all_pp[i, j]] <- S[i, j-1] + R[i, j-1] + s[n_county_idx[i]]
      lambda[pH[i, j-1]] <- dm[i, all_pp[i, j]] / dm[i, all_pp[i, j-1]]
    }

  }

  # for easier monitoring of abundance - long format
  for(i in 1:n_units){
    xn[i] <- N[property_x[i], pp_x[i]]
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
