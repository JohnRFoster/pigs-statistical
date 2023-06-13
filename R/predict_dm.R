calc_log_potential_area <- function(n_mcmc, n_survey, n_trap_snare, n_shooting, method, ts_idx, shooting_idx,
                                    log_rho, log_gamma, log_effort_per, effort_per, p_unique, n_trap_m1){
  log_pi <- log(pi)

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

  return(log_potential_area)
}

calc_log_theta <- function(n_survey, log_survey_area_km2, log_potential_area, beta_p, method, X_p){
  log_theta <- matrix(NA, nrow(beta_p), n_survey)
  for(i in 1:n_survey){

    log_pr_area_sampled <- pmin(log_survey_area_km2[i], log_potential_area[, i])

    b_nodes <- paste0("beta_p[", method[i])
    beta <- beta_p[,grep(b_nodes, colnames(beta_p), fixed = TRUE)]

    # eps_node <- paste0("eps_property_pR[", p_property_idx[i], "]")

    # probability of capture, given that an individual is in the surveyed area
    theta_star <- boot::inv.logit(
      beta[, 1] +
        X_p[i, 1] * beta[, 2] +
        X_p[i, 2] * beta[, 3] +
        X_p[i, 3] * beta[, 4])
    # eps_property_pR[, eps_node] * sigma_property_p)

    log_theta[, i] <- log_pr_area_sampled - log_survey_area_km2[i] + log(theta_star)
  }
  return(log_theta)
}

calc_p <- function(n_mcmc, n_first_survey, n_not_first_survey, first_survey, not_first_survey,
                   log_theta, start, end){
  n_survey <- n_first_survey + n_not_first_survey
  # the probability an individual is captured on the first survey
  p <- matrix(NA, nrow(log_theta), n_survey)
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
  return(p)
}

simulate_nimble_dm <- function(mcmc, flags, constants, data, unit_lookup){
  likelihood_nb <- TRUE
  ls <- list(mcmc = mcmc,
             unit_lookup = unit_lookup)
  ls <- append(ls, flags)
  ls <- append(ls, data)
  ls <- append(ls, constants)

  sim <- with(ls, {

    log_gamma <- mcmc |> select(contains("log_gamma")) |> as.matrix()
    log_rho <- mcmc |> select(contains("log_rho")) |> as.matrix()
    p_unique <- mcmc |> select(contains("p_unique")) |> as.matrix()
    beta_p <- mcmc |> select(contains("beta_p")) |> as.matrix()
    size <- mcmc |> select(contains("size")) |> as.matrix()
    X_p <- as.matrix(X_p)
    n_mcmc <- nrow(beta_p)

    zeta_pp <- if_else(grepl("zeta_pp", process_type), TRUE, FALSE)
    zeta_constant <- !zeta_pp

    phi_pp <- if_else(grepl("phi_pp", process_type), TRUE, FALSE)
    phi_constant <- !phi_pp

    if(phi_constant){
      phi <- mcmc |>
        pull(logit_phi) |>
        ilogit()
    } else if(phi_pp){
      phi <- mcmc |>
        select(contains("logit_phi[")) |>
        as.matrix() |>
        ilogit()
    }

    if(zeta_constant){
      zeta <- mcmc |>
        pull(log_zeta) |>
        exp()
    } else if(zeta_pp){
      zeta <- mcmc |>
        select(contains("log_zeta[")) |>
        as.matrix() |>
        exp()
    }

    # extend process types here

    units <- unit_lookup |> mutate(index = 1:n())
    N <- mcmc |>
      select(contains("xn")) |>
      mutate(iter = 1:n()) |>
      pivot_longer(cols = -iter,
                   names_to = "node") |>
      mutate(index = as.numeric(str_extract(node, "\\d*(?=\\])"))) |>
      left_join(units) |>
      mutate(node = paste0("N[", property_idx, ", ", timestep, "]")) |>
      select(iter, value, node) |>
      pivot_wider(names_from = node,
                  values_from = value) |>
      select(-iter) |>
      as.matrix()

    if(spatial){
      s <- mcmc |> select(contains("s[")) |> as.matrix()
    } else {
      s <- matrix(0, n_mcmc, n_county)
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

    ## -- Data model
    log_potential_area <- calc_log_potential_area(
      n_mcmc, n_survey, n_trap_snare, n_shooting, method, ts_idx, shooting_idx,
      log_rho, log_gamma, log_effort_per, effort_per, p_unique, n_trap_m1)

    log_theta <- calc_log_theta(n_survey, log_survey_area_km2, log_potential_area, beta_p, method, X_p)

    p <- calc_p(
      n_mcmc, n_first_survey, n_not_first_survey,
      first_survey, not_first_survey, log_theta, start, end)

    N_pred <- array(NA, dim = c(n_property, max(n_timesteps), n_mcmc))
    for(i in 1:n_property){

      n_node <- paste0("N[", i, ", 1]")
      s_node <- paste0("size[", i, "]")

      N_pred[i, 1, ] <- N[, n_node]
      dm <- N[, n_node]

      for(j in 2:n_pp_prop[i]){
        ts <- all_pp[i, j]

        if(phi_constant){
          S <- rbinom(n_mcmc, dm, phi)
        } else if (phi_pp){
          p_node <- paste0("logit_phi[", all_pp[i, j-1], "]")
          S <- rbinom(n_mcmc, dm, phi[, p_node])
        }

        if(zeta_constant){
          R <- rpois(n_mcmc, zeta * dm)
        } else if (zeta_pp){
          z_node <- paste0("log_zeta[", all_pp[i, j-1], "]")
          R <- rpois(n_mcmc, zeta[, z_node] * dm)
        }

        dm <- S + R

        if(ts %in% timestep[i, ]){
          ts_idx <- which(ts == timestep[i, ])
          N_pred[i, ts_idx, ] <- dm
          dm <- N[, paste0("N[", i, ", ", ts_idx, "]")]
        }

      }

    }

    # likelihood - predicted number removed given data model
    y_pred <- matrix(NA, n_mcmc, n_survey)
    for(i in 1:n_survey){
      z <- N_pred[p_property_idx[i], p_timestep_idx[i], ] - y_sum[i]
      z <- pmax(0, z)
      if(likelihood_nb){
        w <- size[,paste0("size[", p_county_idx[i], "]")]
        py <- w / ((p[,i] * z) + w)
        y_pred[,i] <- rnbinom(n_mcmc, w, py)
        if(any(is.na(y_pred[,i]))) stop("Negative binomial likelihood error")
      } else if(likelihood_binom){
        y_pred[,i] <- rbinom(n_mcmc, z, p[,i])
      }
    }

    list(p = p,
         yp = y_pred,
         N = N_pred)
  })
  return(sim)
}
