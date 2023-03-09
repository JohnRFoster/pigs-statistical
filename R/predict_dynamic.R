simulate_nimble_dynamic <- function(mcmc, flags, constants, data, unit_lookup){

  ls <- list(mcmc = mcmc,
             unit_lookup = unit_lookup)
  ls <- purrr::prepend(ls, flags)
  ls <- purrr::prepend(ls, data)
  ls <- purrr::prepend(ls, constants)

  sim <- with(ls, {

    log_gamma <- mcmc |> select(contains("log_gamma")) |> as.matrix()
    log_rho <- mcmc |> select(contains("log_rho")) |> as.matrix()
    p_unique <- mcmc |> select(contains("p_unique")) |> as.matrix()
    beta_p <- mcmc |> select(contains("beta_p")) |> as.matrix()
    sigma_proc <- mcmc |> mutate(sigma_proc = 1/sqrt(tau_proc)) |> pull(sigma_proc)
    X_p <- as.matrix(X_p)
    n_mcmc <- nrow(beta_p)

    if(process_type %in% c("jamiesonBrooks", "dennisTaper")){
      beta_r <- mcmc |> pull(beta_r)
    }

    log_r_mu <- mcmc |> pull("log_r_mu")
    r_t <- matrix(NA, n_mcmc, n_pp)
    for(i in 1:n_pp){
      r_t[,i] <- log_r_mu
    }

    units <- unit_lookup |> mutate(index = 1:n())
    x <- mcmc |>
      select(contains("xn")) |>
      mutate(iter = 1:n()) |>
      pivot_longer(cols = -iter,
                   names_to = "node") |>
      mutate(index = as.numeric(str_extract(node, "\\d*(?=\\])"))) |>
      left_join(units) |>
      mutate(node = paste0("x[", property_idx, ", ", timestep, "]")) |>
      select(iter, value, node) |>
      pivot_wider(names_from = node,
                  values_from = value) |>
      select(-iter) |>
      as.matrix()

    if(spatial){
      s <- mcmc |> select(contains("s[")) |> as.matrix()
    } else {
      s <- matrix(1, n_mcmc, n_county)
      colnames(s) <- paste0("s[", 1:n_county, "]")
    }

    if(property_obs_effect){
      eps_property_pR <- mcmc |> select(contains("eps_property_pR")) |> as.matrix()
      sigma_property_p <- mcmc |> pull("sigma_property_p")
    } else {
      eps_property_pR <- matrix(1, n_mcmc, n_property)
      colnames(eps_property_pR) <- paste0("eps_property_pR[", 1:n_property, "]")
      sigma_property_p <- 0
    }

    log_potential_area <- matrix(NA, n_mcmc, n_survey)

    # potential area impacted by trap or snare removals
    for(i in 1:n_trap_snare){
      r_node <- paste0("log_rho[", method[ts_idx[i]], "]")
      g_node <- paste0("log_gamma[", method[ts_idx[i]], "]")
      p_node <- paste0("p_unique[", method[ts_idx[i]], "]")

      log_potential_area[, ts_idx[i]] <- log_pi +
        (2 * (log_rho[, r_node] + log_effort_per[ts_idx[i]] - log(exp(log_gamma[, g_node]) + effort_per[ts_idx[i]]))) +
        log(1 + (p_unique[, p_node] * n_trap_m1[ts_idx[i]]))
    }

    # potential area impacted by shooting removals
    for(i in 1:n_shooting){
      r_node <- paste0("log_rho[", method[shooting_idx[i]], "]")
      g_node <- paste0("log_gamma[", method[shooting_idx[i]], "]")
      p_node <- paste0("p_unique[", method[shooting_idx[i]], "]")

      log_potential_area[,shooting_idx[i]] <- log_rho[, r_node] +
        log_effort_per[shooting_idx[i]] -
        log(exp(log_gamma[, g_node]) + effort_per[shooting_idx[i]]) +
        log(1 + (p_unique[, p_node] * n_trap_m1[shooting_idx[i]]))
    }

    log_theta <- matrix(NA, n_mcmc, n_survey)
    for(i in 1:n_survey){

      log_pr_area_sampled <- pmin(log_survey_area_km2[i], log_potential_area[, i])

      b_nodes <- paste0("beta_p[", method[i])
      beta <- beta_p[,grep(b_nodes, colnames(beta_p), fixed = TRUE)]

      # eps_node <- paste0("eps_property_pR[", p_property_idx[i], "]")

      # probability of capture, given that an individual is in the surveyed area
      theta_star <- boot::inv.logit(
        beta_p[, 1] +
        X_p[i, 1] * beta_p[, 2] +
        X_p[i, 2] * beta_p[, 3] +
        X_p[i, 3] * beta_p[, 4])
        # eps_property_pR[, eps_node] * sigma_property_p)

      log_theta[, i] <- log_pr_area_sampled - log_survey_area_km2[i] + log(theta_star)

    }

    # the probability an individual is captured on the first survey
    p <- matrix(NA, n_mcmc, n_survey)
    for(i in 1:n_first_survey){
      p[, first_survey[i]] <- exp(log_theta[, first_survey[i]])
    }

    # the probability an individual is captured after the first survey
    for(i in 1:n_not_first_survey){
      j_seq <- start[not_first_survey[i]]:end[not_first_survey[i]]

      if(length(j_seq) == 1){
        cum_p <- log(1-exp(log_theta[, j_seq]))
      } else {
        cum_p <- apply(log_theta[, j_seq], 1, function(x) sum(log(1-exp(x))))
      }

      p[, not_first_survey[i]] <- exp(log_theta[, not_first_survey[i]] + cum_p)
    }

    N <- array(NA, dim = c(n_property, max(n_timesteps), n_mcmc))
    yp <- array(NA, dim = c(n_property, max(n_timesteps), max(n_reps, na.rm = T), n_mcmc))
    for(i in 1:n_property){

      x_node <- paste0("x[", i, ", 1]")

      # likelihood at first rep in first PP
      N[i, 1, ] <- rpois(n_mcmc, exp(x[,x_node]))
      yp[i, 1, 1, ] <- rbinom(n_mcmc, N[i, 1, ], p[, H[i, 1, 1]])
      z <- N[i, 1, ] - yp[i, 1, 1, ]

      # likelihood for reps after first pass in fist PP
      for(j in 2:n_reps[i, 1]){

        yp[i, 1, j, ] <- rbinom(n_mcmc, z, p[, H[i, 1, j]])
        removed <- apply(yp[i, 1, 1:j, ], 2, sum)
        z <- N[i, 1, ] - removed

      }

      for(t in 2:n_timesteps[i]){

        # population growth across time steps
        l_seq <- (timestep[i, t-1]+1):timestep[i, t]
        if(length(l_seq) == 1){
          log_r <- r_t[,l_seq]
        } else {
          log_r <- apply(r_t[,l_seq], 1, sum)
        }

        # sampling unit is the primary period in a property
        x_node <- paste0("x[", i, ", ", t-1, "]")
        s_node <- paste0("s[", county_idx[i], "]")

        if(process_type == "exponential"){
          mu <- x[, x_node] +
            log_r +
            log(s[, s_node])
        } else if(process_type == "ricker"){
          mu <-  x[, x_node] +
            exp(log_r + log(1 - exp(x[, x_node] - log_k[i]))) +
            log(s[, s_node])
        } else if(process_type == "gompertz"){
          mu <- x[, x_node] +
            exp(log_r +
                  log(1 - exp(log(log(exp(x[, x_node]) + 1)) - log(log(exp(log_k[i]) + 1))))) +
            log(s[, s_node])
        } else if(process_type == "jamiesonBrooks"){
          mu <- x[i, x_node] +
            log_r +
            beta_r*exp(x[i, x_node]) +
            log(s[, s_node])
        } else if(process_type == "dennisTaper"){
          mu <- log_r +
            (1 + beta_r)*x[i, x_node] +
            log(s[, s_node])
        }

        # mu <- pmin(25, mu)
        x_pred <- rnorm(n_mcmc, mu, sigma_proc)
        N[i, t, ] <- rpois(n_mcmc, exp(x_pred))

        # likelihood at first rep in > 1 PP
        yp[i, t, 1, ] <- rbinom(n_mcmc, N[i, t, ], p[, H[i, t, 1]])
        z <- N[i, t, ] - yp[i, t, 1, ]

        # likelihood for reps after first pass in > 1 PP
        for(j in 2:n_reps[i, t]){

          yp[i, t, j, ] <- rbinom(n_mcmc, z, p[, H[i, t, j]])
          removed <- apply(yp[i, t, 1:j, ], 2, sum)
          z <- N[i, t, ] - removed

        }
      }
    }
    list(p = p,
         yp = yp,
         N = N)
  })
  return(sim)
}
