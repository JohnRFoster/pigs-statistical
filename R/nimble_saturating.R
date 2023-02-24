

##-------- parameters --------
# beta, vector, non time varying coefficients - process model
# beta_p, vector, non time varying coefficients - observation model
# eps_property, vector, abundance ??
# sigma_property, scalar, ??
# sigma_short, scalar, ??
# z_shortR, matrix, B-spline spatial coefficient

# TODO finish list at top explaining the function and dimenstion of each param


modelName <- "saturating"

modelCode <- nimbleCode({
  ### priors

  # non time varying coefficients - process model
  for(i in 1:m_n){
    beta[i] ~ dnorm(0, sd = 1)
  }

  # non time varying coefficients - observation model
  for(i in 1:m_p){
    beta_p[i] ~ dnorm(0, sd = 1)
  }

  for(i in 1:n_method){
    p_unique[i] ~ dbeta(2, 2)
    rho[i] ~ dgamma(0.1, 0.1)
    # rho[i] ~ dlnorm(log_rho_prior[i, 1], log_rho_prior[i, 2])
  }

  for(i in 1:n_property){
    eps_property_pR[i] ~ dnorm(0, sd = 1)
    eps_propertyR[i] ~ dnorm(0, sd = 1)
  }

  eps_property[1:n_property] <- eps_propertyR[1:n_property] * sigma_property

  # sigma_car ~ dunif(0, 100)   # prior for variance components based on Gelman (2006)
  # tau_car <- 1 / sigma_car^2
  tau_car ~ dexp(0.01)
  # sigma_short ~ dexp(0.01)
  sigma_property ~ dexp(0.01)
  sigma_property_p ~ dexp(0.01)
  # sigma_st0 ~ dexp(1) # scale parameter at time 1 for spatially uncorrelated temporal autocorrelation
  # sigma_st ~ dexp(1) # scale parameter at time > 1 for spatially uncorrelated temporal autocorrelation
  eta ~ dbeta(9, 1) # temporal autocorrelation parameter (AR1)



  for(i in 1:3){
    gamma[i] ~ dgamma(1, 1)
  }
  gamma[4] ~ dgamma(gamma_prior[1, 1], gamma_prior[1, 2])         # snare
  gamma[5] ~ dgamma(gamma_prior[2, 1], gamma_prior[2, 2])         # trap

  for(i in 1:n_survey){

    log(potential_area[i]) <- log(1 + p_unique[method[i]] * n_trap_m1[i]) + # all methods
      ((log_pi + 2 * (log(rho[method[i]]) + log_effort_per[i] - log(gamma[method[i]] + effort_per[i]))) * trap_snare_ind[i]) + # traps and snares only
      ((log(rho[method[i]]) + log_effort_per[i]) * shooting_ind[i]) # shooting and aerial only

    # area_constraint[i] ~ dconstraint(potential_area[i] <= survey_area_km2[i])
    pr_area_sampled[i] <- min(survey_area_km2[i], potential_area[i])

    # probability of capture, given that an individual is in the surveyed area
    logit(theta_star[i]) <- (X_p[i, 1:m_p] %*% beta_p[1:m_p] + # beta_p[1:5] are random intercepts by method
      eps_property_pR[p_property_idx[i]] * sigma_property_p)[1,1]
      # eps_property_pR[p_property_idx[i]]

    theta[i] <- (pr_area_sampled[i] / survey_area_km2[i]) * theta_star[i]

    # the probability an individual is captured
    log(p[i]) <- log(theta[i]) +
      sum(log(1 - theta[start[i]:end[i]])) * not_first_survey[i] # if > 1st survey

    # likelihood
    y[i] ~ dpois(lambda[survey_idx[i]] * p[i])
  }


  # auto regressive prior on the first basis vector
  for(i in 1:m_short){
    z_shortR[i, 1:n_county] ~ dcar_normal(
      adj = adj[1:n_edges],
      weights = weights[1:n_edges],
      num = num[1:n_county],
      tau = tau_car,
      # tau = 1,
      zero_mean = 1
    )
  }

  # loop through autocorrelation
  z_short[1, 1:n_county] <- z_shortR[1, 1:n_county]
  for(i in 2:m_short){
    z_short[i, 1:n_county] <- z_shortR[i - 1, 1:n_county] +
      # z_shortR[i, 1:n_county] * sigma_short
      z_shortR[i, 1:n_county]
  }

  # for(t in 1:n_timestep){
  #   for(i in 1:n_county){
  #     eps_stR[t, i] ~ dnorm(0, sd = 1)
  #   }
  # }

  # create the spatio-temporal adjustment
  # eps_st[1, 1:n_county] <- eps_stR[1, 1:n_county] * sigma_st0
  # for(t in 2:n_timestep){
  #   # at subsequent time steps, there is autocorrelation (eta) between time steps
  #   eps_st[t, 1:n_county] <- eta * eps_st[t - 1, 1:n_county] +
  #     eps_stR[t, 1:n_county] * sigma_st
  # }



  # expected pig abundance
  for (i in 1:n_st) {
    log(lambda[i]) <- (X[i, 1:m_n] %*% beta[1:m_n] + # non time varying covariates
      X_short[i, 1:m_short] %*% z_short[1:m_short, county_idx[i]] + # B-spline vectors * basis vector coefficients
      eps_property[property[i]] + # property adjustment
      # eps_st[timestep[i], county_idx[i]] + # spatio-temporal adjustment
      log_area_km2[i])[1,1] # area offset
  }

})

simulate_saturating <- function(post, constants, data){

  within_data <- purrr::prepend(data, constants)

  with(within_data,{
    post <- as_tibble(post)
    n_mcmc <- nrow(post)
    sigma_property <- post |> pull("sigma_property") |> as.numeric()
    sigma_property_p <- post |> pull("sigma_property_p") |> as.numeric()
    sigma_short <- post |> pull("sigma_short") |> as.numeric()
    eps_propertyR <- post |> select(contains("eps_propertyR")) |> as.data.frame()
    eps_property_pR <- post |> select(contains("eps_property_pR")) |> as.data.frame()
    p_unique <- post |> select(contains("p_unique")) |> as.data.frame()
    rho <- post |> select(contains("rho")) |> as.data.frame()
    gamma <- post |> select(contains("gamma")) |> as.data.frame()
    beta_p <- post |> select(contains("beta_p")) |> as.data.frame()
    beta <- post |> select(contains("beta[")) |> as.data.frame()
    z_shortR <- post |> select(contains("z_shortR")) |> as.data.frame()

    # storage
    potential_area <- pr_area_sampled <- N <-  matrix(NA, n_mcmc, n_survey)
    p <- theta_star <- theta <- potential_area
    z_short <- array(NA, dim = c(m_short, n_county, n_mcmc))
    lambda <- matrix(NA, n_mcmc, n_st)

    eps_property <- eps_propertyR * sigma_property

    message("calculate potential area for each survey")
    pb <- txtProgressBar(max = n_survey, style = 3)
    for(i in 1:n_survey){

      potential_area[1:n_mcmc,i] <- exp(log(1 + p_unique[method[i]] * n_trap_m1[i]) + # all methods
                                          ((log_pi + 2 * (log(rho[method[i]]) + log_effort_per[i] - log(gamma[method[i]] + effort_per[i]))) * trap_snare_ind[i]) + # traps and snares only
                                          ((log(rho[method[i]]) + log_effort_per[i])) * shooting_ind[i])[,1] # shooting and aerial only

      # area_constraint[i] ~ dconstraint(potential_area[i] <= survey_area_km2[i])
      pr_area_sampled[,i] <- pmin(survey_area_km2[i], potential_area[,i])

      setTxtProgressBar(pb, i)
    }
    close(pb)

    message("Calculate lambda and p")
    pb <- txtProgressBar(max = n_mcmc, style = 3)
    for(m in 1:n_mcmc){
      for(i in 1:n_county){
        z_node1 <- paste0("z_shortR[1, ", i, "]")
        z_short[1, i, m] <- z_shortR[m, z_node1]
        for(j in 2:m_short){
          z_node2 <- paste0("z_shortR[", j-1, ", ", i, "]")
          z_node3 <- paste0("z_shortR[", j, ", ", i, "]")
          z_short[j, i, m] <- z_shortR[m, z_node2] +
            z_shortR[m, z_node3] * sigma_short[m]
        }
      }

      for(i in 1:n_st){
        lambda[m, i] <- exp(crossprod(X[i,], as.numeric(beta[m,])) + # non time varying covariates
                              X_short[i, ] %*% z_short[, county_idx[i], m] + # B-spline vectors * basis vector coefficients
                              eps_property[m, property[i]] + # property adjustment
                              # eps_st[timestep[i], county_idx[i]] + # spatio-temporal adjustment
                              log_area_km2[i])[1,1] # area offset
      }

      # probability of capture, given that an individual is in the surveyed area
      for(i in 1:n_survey){
        theta_star[m, i] <- boot::inv.logit((crossprod(X_p[i, ], as.numeric(beta_p[m, ])) + # beta_p[1:5] are random intercepts by method
                                               eps_property_pR[m, p_property_idx[i]] * sigma_property_p[m])[,1])

        theta[m,i] <- (pr_area_sampled[m,i] / survey_area_km2[i]) * theta_star[m,i]

        # the probability an individual is captured
        p[m, i] <- exp(log(theta[m, i]) +
                         sum(log(1 - theta[m ,start[i]:end[i]])) * not_first_survey[i]) # if > 1st survey
      }
      setTxtProgressBar(pb, m)
    }
    close(pb)

    # likelihood
    for(i in 1:n_survey){
      N[, i] <- rpois(n_mcmc, lambda[, survey_idx[i]] * p[, i])
    }
    return(list(
      N = N,
      lambda = lambda,
      p = p
    ))
  })
}
