modelCode <- nimbleCode({

  # priors
  for(i in 1:n_method){
    log_rho[i] ~ dnorm(0, tau = 0.1)
    log_gamma[i] ~ dnorm(0, tau = 0.1)
    p_mu[i] ~ dnorm(0, tau = 1)
    logit(p_unique[i]) <- p_mu[i]

    # non time varying coefficients - observation model
    for(j in 1:n_beta_p){
      beta_p[i, j] ~ dnorm(0, tau = 1)
    }
  }

  tau_proc ~ dgamma(0.01, 0.01)


  for(i in 1:n_trap_snare){
    log_potential_area[ts_idx[i]] <- log_pi +
      (2 * (log_rho[method[ts_idx[i]]] + log_effort_per[ts_idx[i]] - log(exp(log_gamma[method[ts_idx[i]]]) + effort_per[ts_idx[i]]))) +
      log(1 + (p_unique[method[ts_idx[i]]] * n_trap_m1[ts_idx[i]]))
  }

  for(i in 1:n_shooting){
    log_potential_area[shooting_idx[i]] <- log_rho[method[shooting_idx[i]]] +
      log_effort_per[shooting_idx[i]] -
      log(exp(log_gamma[method[shooting_idx[i]]]) + effort_per[shooting_idx[i]]) +
      log(1 + (p_unique[method[shooting_idx[i]]] * n_trap_m1[shooting_idx[i]]))
  }

  for(i in 1:n_survey){

    # area_constraint[i] ~ dconstraint(potential_area[i] <= survey_area_km2[i])
    log_pr_area_sampled[i] <- min(log_survey_area_km2[i], log_potential_area[i])

    # probability of capture, given that an individual is in the surveyed area
    logit(theta_star[i]) <- beta_p[method[i], 1] + inprod(X_p[i, 1:m_p], beta_p[method[i], 2:n_beta_p])
                               # eps_property_pR[p_property_idx[i]] * sigma_property_p

    log_theta[i] <- log_pr_area_sampled[i] - log_survey_area_km2[i] + log(theta_star[i])

  }

  # the probability an individual is captured on the first survey
  for(i in 1:n_first_survey){
    log(p[first_survey[i]]) <- log_theta[first_survey[i]]
  }

  # the probability an individual is captured after the first survey
  for(i in 1:n_not_first_survey){
    log(p[not_first_survey[i]]) <- log_theta[not_first_survey[i]] +
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
      s[non_islands[i]] <- log(s_car[i])
    }
    for(i in 1:n_islands){
      s_island[i] ~ dnorm(0, tau = 1)
      s[islands[i]] <- s_island[i]
    }
  } else {
    s[1:n_county] <- 1
  }


  # estimate population growth (constant)
  log_r_mu ~ dnorm(0, tau = 1)
  for(i in 1:(n_pp)){
    log_r_t[i] <- log_r_mu
  }

  for(i in 1:n_property){

    # eps_property_pR[i] ~ dnorm(0, tau = 1) # property effect in observation model

    y[i, 1, 1] ~ dbinom(p[H[i, 1, 1]], N[i, 1])
    z[i, 1, 1] <- N[i, 1] - y[i, 1, 1]
    N[i, 1] ~ dpois(lambda_N[i, 1])
    log(lambda_N[i, 1]) <- x[i, 1]
    x[i, 1] ~ dnorm(0, tau = 1)

    # y[i, 1, 1:n_reps[i, 1]] ~ dNmixture_v(lambda_N, p[H[i, 1, 1:n_reps[i, 1]]], Nmin = -1, Nmax = -1, len = length(1:n_reps[i, 1]))

    for(j in 2:n_reps[i, 1]){

      y[i, 1, j] ~ dbinom(p[H[i, 1, j]], z[i, 1, j-1])
      z[i, 1, j] <- N[i, 1] - sum(y[i, 1, 1:j])

    }


    for(t in 2:n_timesteps[i]){

      y[i, t, 1] ~ dbinom(p[H[i, t, 1]], N[i, t])
      z[i, t, 1] <- N[i, t] - y[i, t, 1]

      for(j in 2:n_reps[i, t]){

        y[i, t, j] ~ dbinom(p[H[i, t, j]], z[i, t, j-1])
        z[i, t, j] <- N[i, t] - sum(y[i, t, 1:j])

      }

      # sampling unit is the primary period (one month) in a property
      N[i, t] ~ dpois(lambda_N[i, t])
      log(lambda_N[i, t]) <- x[i, t]
      x[i, t] ~ dnorm(ex[i, t], tau = tau_proc)

      if(exponential){
        ex[i, t] <- x[i, t-1] +
          log_r[i, t] +
          s[county_idx[i]]
      }

      if(ricker){
        ex[i, t] <- x[i, t-1] +
          exp(log_r[i, t] + log(1 - exp(x[i, t-1] - log_k[i]))) +
          s[county_idx[i]]
      }

      if(gompertz){
        ex[i, t] <- x[i, t-1] +
          exp(log_r[i, t] +
                log(1 - exp(log(log(exp(x[i, t-1]) + 1)) - log(log(exp(log_k[i]) + 1))))) +
          s[county_idx[i]]
      }

      # population growth across time steps
      log_r[i, t] <- sum(log_r_t[(timestep[i, t-1]+1):timestep[i, t]])

    }
  }

  for(i in 1:n_units){
    xn[i] <- x[property_x[i], timestep_x[i]]
  }

})
