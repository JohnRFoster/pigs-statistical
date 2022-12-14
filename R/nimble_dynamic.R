modelCode <- nimbleCode({

  # priors
  for(i in 1:3){
    gamma[i] ~ dgamma(1, 1)
  }
  gamma[4] ~ dgamma(0.1, 0.1)         # snare
  gamma[5] ~ dgamma(0.1, 0.1)         # trap

  for(i in 1:n_method){
    rho[i] ~ dgamma(0.1, 0.1)
    p_unique[i] ~ dbeta(2, 2)

    # non time varying coefficients - observation model
    for(j in 1:n_beta_p){
      beta_p[i, j] ~ dnorm(0, tau = 1)
    }
  }

  # sigma_property_p ~ dexp(0.01)
  sigma_phi ~ dunif(0, 100)
  tau_phi <- 1 / sigma_phi^2

  for(i in 1:n_survey){

    log(potential_area[i]) <- log(1 + p_unique[method[i]] * n_trap_m1[i]) + # all methods
      ((log_pi + 2 * (log(rho[method[i]]) + log_effort_per[i] - log(gamma[method[i]] + effort_per[i]))) * trap_snare_ind[i]) + # traps and snares only
      ((log(rho[method[i]]) + log_effort_per[i]) * shooting_ind[i]) # shooting and aerial only

    # area_constraint[i] ~ dconstraint(potential_area[i] <= survey_area_km2[i])
    pr_area_sampled[i] <- min(survey_area_km2[i], potential_area[i])

    # probability of capture, given that an individual is in the surveyed area
    # theta_star[i] ~ dunif(0, 1)
    logit(theta_star[i]) <- beta_p[method[i], 1] + inprod(X_p[i, 1:m_p], beta_p[method[i], 2:n_beta_p])
                               # eps_property_pR[p_property_idx[i]] * sigma_property_p

    theta[i] <- (pr_area_sampled[i] / survey_area_km2[i]) * theta_star[i]

    # method[i] ~ dcat(method_prob[1:n_method])

    # the probability an individual is captured
    log(p[i]) <- log(theta[i]) +
      sum(log(1 - theta[start[i]:end[i]])) * not_first_survey[i] # if > 1st survey

  }

  # sigma_car ~ dunif(0, 100)   # prior for variance components based on Gelman (2006)
  # tau_car <- 1 / sigma_car^2
  # s[1:n_county] ~ dcar_normal(
  #   adj = adj[1:n_edges],
  #   weights = weights[1:n_edges],
  #   num = num[1:n_county],
  #   tau = tau_car,
  #   # tau = 1,
  #   zero_mean = 1
  # )

  # estimate population growth (constant)
  lambda_mu ~ dexp(1)
  for(i in 1:(n_pp)){
    lambda_t[i] <- lambda_mu
  }

  for(i in 1:n_property){

    # eps_property_pR[i] ~ dnorm(0, tau = 1) # property effect in observation model

    y[i, 1, 1] ~ dbinom(p[H[i, 1, 1]], N[i, 1])
    N[i, 1] ~ dpois(lambda_N[i, 1])
    log(lambda_N[i, 1]) <- phi[i, 1]
    phi[i, 1] ~ dnorm(0, tau = 1)

    for(t in 2:n_timesteps[i]){

      y[i, t, 1] ~ dbinom(p[H[i, t, 1]], N[i, t])
      z[i, t, 1] <- N[i, t] - y[i, t, 1]

      for(j in 2:n_reps[i, t]){

        y[i, t, j] ~ dbinom(p[H[i, t, j]], z[i, t, j-1])
        z[i, t, j] <- N[i, t] - sum(y[i, t, 1:j])

      }

      # sampling unit is the primary period (one month) in a property
      N[i, t] ~ dpois(lambda_N[i, t])
      log(lambda_N[i, t]) <- phi[i, t]
      phi[i, t] ~ dnorm(phi_ex[i, t], tau = tau_phi)
      phi_ex[i, t] <- phi[i, t-1] * lambda[i, t] #+ s[county_idx[i]]

      # population growth across time steps
      log(lambda[i, t]) <- sum(log(lambda_t[(timestep[i, t-1]+1):timestep[i, t]]))

    }
  }




})
