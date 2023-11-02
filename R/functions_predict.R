
#'@description Calculate the potential search area from posterior samples when traps or snares are used
#'@param log_rho vector of mcmc samples for rho (log scale)
#'@param log_gamma vector of mcmc samples for gamma (log scale)
#'@param p_unique vector of mcmc samples for p
#'@param effort_per the effort per trap/snare
#'@param n_trap_m1 the number of traps/snares used minus 1

trap_snare_lpa <- function(log_rho, log_gamma, p_unique, effort_per, n_trap_m1){
  log(pi) +
    (2 * (log_rho + log(effort_per) -
            log(exp(log_gamma) + effort_per))) +
    log(1 + (p_unique * n_trap_m1))
}

#'@description Calculate the potential search area from posterior samples when aerial methods (fixed wing, helicopter) are used
#'@param log_rho vector of mcmc samples for rho (log scale)
#'@param effort_per the effort per trap/snare

aerial_lpa <- function(log_rho, effort_per){
  log_rho + log(effort_per)
}

#'@description Calculate the potential search area from posterior samples when firearms are used
#'@param log_rho vector of mcmc samples for rho (log scale)
#'@param p_unique vector of mcmc samples for p
#'@param effort_per the effort per trap/snare
#'@param n_trap_m1 the number of traps/snares used minus 1

firearms_lpa <- function(log_rho, p_unique, effort_per, n_trap_m1){
  log_rho +
    log(effort_per) -
    log(1 + (p_unique * n_trap_m1))
}

#'@description Calculate the potential search area from posterior samples for a vector of methods (used in forecasting code)
#'@param log_rho vector of mcmc samples for rho (log scale), iterations (rows) by method (columns)
#'@param log_gamma vector of mcmc samples for gamma (log scale), iterations (rows) by method (columns)
#'@param p_unique vector of mcmc samples for p, iterations (rows) by method (columns)
#'@param effort_per the effort per unit
#'@param n_trap_m1 the number of units used minus 1

lpa <- function(method, log_rho, log_gamma, p_unique, effort_per, n_trap_m1){
  n_mcmc <- length(log_rho)
  n_reps <- length(method)
  log_potential_area <- matrix(NA, n_mcmc, n_reps)
  for(j in seq_along(method)){
    lr <- log_rho[,method[j]]
    ep <- effort_per[j]
    tcm1 <- n_trap_m1[j]
    if(method[j] == 1){
      p <- p_unique[,method[j]]
      log_potential_area[,j] <- firearms_lpa(lr, p, ep, tcm1)
    } else if(method[j] == 2 | method[j] == 3){
      log_potential_area[,j] <- aerial_lpa(lr, ep)
    } else if(method[j] == 4 | method[j] == 5){
      lg <- log_gamma[,method[j]-3]
      p <- p_unique[,method[j]-2]
      log_potential_area[,j] <- trap_snare_lpa(lr, lg, p, ep, tcm1)
    }
  }
  return(log_potential_area)
}

#'@description calculate the capture probability for a removal event
#'@param X matrix of county-level land cover covariates, first column is 1's for an intercept term
#'@param beta vector of coeffiencents
#'@param log_potential_area matrix of potential search area, iterations (rows) by replicate (columns)
#'@param area_property scalar value of the area (mk^2) of the property

calc_p <- function(X, beta1, beta, log_potential_area, area_property, method){
  n_mcmc <- nrow(beta)
  n_reps <- ncol(log_potential_area)
  p <- matrix(NA, n_mcmc, n_reps)
  for(mc in seq_len(n_mcmc)){
    log_theta <- numeric(n_reps)
    for(j in seq_len(n_reps)){

      # base probability of capture given an individual is the first survey
      log_theta[j] <- log(ilogit(beta1[mc] + inprod(X, beta[mc,]))) +
        pmin(0, log_potential_area[mc, j] - log(area_property))

      # the probability an individual is captured on the first survey
      if(j == 1){
        p[mc, j] <- exp(log_theta[j])
      } else {
        # the probability an individual is captured after the first survey
        for(k in 2:n_reps){
          p[mc, j] <- exp(log_theta[1] +
                            sum(log(1 - exp(log_theta[1:(j-1)]))))
        }
      }
    }
  }
  return(p)
}


data_posteriors <- function(samples, constants, data){

  require(dplyr)
  require(nimble)

  D <- append(constants, data)

  post <- with(D, {

    xH <- tibble(
      property = p_property_idx,
      pp = p_pp_idx
    ) |>
      group_by(property, pp) |>
      mutate(number = cur_group_id()) |>
      pull(number)

    log_potential_area <- matrix(NA, nrow(samples), n_survey)
    y <- matrix(NA, nrow(samples), n_survey)
    log_theta <- matrix(NA, nrow(samples), n_survey)

    pattern <- "(?<!beta_)p\\[\\d*\\]"
    p_detect <- str_detect(colnames(samples), pattern)
    if(any(p_detect)){
      calc_p <- FALSE
      p <- samples[,which(p_detect)]
    } else {
      calc_p <- TRUE
      p <- matrix(NA, nrow(samples), n_survey)
    }

    log_rho <- samples[, grep("log_rho", colnames(samples))]
    log_gamma <- samples[, grep("log_gamma", colnames(samples))]
    p_unique <- ilogit(samples[, grep("p_mu", colnames(samples))])
    beta1 <- samples[, grep("beta1", colnames(samples))]
    beta_p <- samples[, grep("beta_p", colnames(samples))]

    for(i in 1:n_survey){

      if(method[i] == 1){
        log_potential_area[, i] <- firearms_lpa(
          log_rho = log_rho[, method[i]],
          p_unique = p_unique[, method[i]],
          effort_per = effort_per[i],
          n_trap_m1 = n_trap_m1[i]
        )
      } else if(method[i] == 2 | method[i] == 3){
        log_potential_area[, i] <- aerial_lpa(
          log_rho = log_rho[, method[i]],
          effort_per = effort_per[i]
        )
      } else if(method[i] == 4 | method[i] == 5){
        log_potential_area[, i] <- trap_snare_lpa(
          log_rho = log_rho[, method[i]],
          log_gamma = log_gamma[, method[i] - 3],
          p_unique = p_unique[, method[i] - 2],
          effort_per = effort_per[i],
          n_trap_m1 = n_trap_m1[i]
        )
      }

      if(calc_p){

        pb <- txtProgressBar(max = nrow(samples), style = 3)
        for(m in 1:nrow(samples)){
          M <- samples[m,]

          # base probability of capture given an individual is the first survey
          # TODO fix to work with beta1 and beta_p by method if p is not saved from mcmc
          log_theta[m, ] <- log(ilogit(X_p %*% M[grep("beta_p", names(M))])) +
            pmin(0, log_potential_area[m, ] - log_survey_area_km2[i])

          # the probability an individual is captured on the first survey
          p[m, first_survey] <- exp(log_theta[m, first_survey])

          # the probability an individual is captured after the first survey
          for(i in 1:n_not_first_survey){
            p[m, not_first_survey[i]] <- exp(log_theta[m, start[not_first_survey[i]]] +
                                               sum(log(1 - exp(log_theta[m, start[not_first_survey[i]]:end[not_first_survey[i]]]))))
          }
          setTxtProgressBar(pb, m)
        }
        close(pb)
      }
      N <- as.numeric(samples[,paste0("xn[", xH[i], "]")])
      y[, i] <- rpois(length(N), N * p[,i])
    }

    list(
      y = y,
      p = p,
      potential_area = exp(log_potential_area),
      theta = exp(log_theta)
    )

  })
  return(post)
}

remove <- function(i, t, j, effort, effort_property, reps_property, methods,
                   log_rho, log_gamma, p_unique, X, beta1, beta_p, new_abundance){
  id <- j
  n_mcmc <- length(beta1)
  if(effort == 0){
    pass <- 0
    ep <- 0
    tc_m1 <- 0
    m <- 0
    predicted_take <- tibble(
      ens = 1:n_mcmc,
      take = 0,
      pass = pass)
  } else {
    E <- effort_property |>
      filter(q == effort)

    n <- reps_property |>
      filter(q == effort) |>
      pull(n_reps)

    # removal methods based on frequency of use in a given property
    m <- sample(methods, n, replace = TRUE)

    # extract effort attributes from property
    get_e <- function(e){
      lapply(m, function(i){
        E |> filter(method == i)
      }) |>
        bind_rows() |>
        pull(e)
    }

    ep <- get_e("effort_per")
    tc_m1 <- get_e("trap_count") - 1

    # the search area given method and effort
    log_potential_area <- lpa(
      method = m,
      log_rho = log_rho,
      log_gamma = log_gamma,
      p_unique = p_unique,
      effort_per = ep,
      n_trap_m1 = tc_m1
    )

    # probability of capture for each pass
    p <- calc_p(
      X = X,
      beta1 = beta1,
      beta = beta_p,
      log_potential_area = log_potential_area,
      area_property = E$area_property[1]
    )

    N <- new_abundance |>
      filter(id == j) |>
      pull(abundance)

    # remove
    y <- matrix(NA, n_mcmc, n)
    for(j in seq_len(n)){
      y[,j] <- rpois(n_mcmc, N * p[,j])
    }

    pass <- 1:n
    colnames(y) <- pass
    predicted_take <- y |>
      as_tibble() |>
      pivot_longer(cols = everything(),
                   names_to = "pass",
                   values_to = "take") |>
      mutate(pass = as.integer(pass)) |>
      group_by(pass) |>
      mutate(ens = 1:n()) |>
      ungroup()
  }

  y_rep <- tibble(
    pass = pass,
    property = i,
    start_pp = t, # the PP we started from
    PPNum = t+1, # the PP we are forecasting into (the PP pigs are removed)
    horizon = 1, # the number of PPs being forecasted across
    method = m,
    effort_per = ep,
    trap_count = tc_m1,
    effort_scenario = effort
  )

  left_join(predicted_take, y_rep) |>
    suppressMessages()
}

get_IC_posterior <- function(t, observations, n_samples){

    n_id <- observations |>
      filter(PPNum == t) |>
      pull(n_id)
    n_node <- paste0("xn[", n_id, "]")
    tibble(
      N = n_samples[,grep(n_node, colnames(n_samples), fixed = TRUE)],
      effort_scenario = NA,
      id = 1,
      ens = 1:nrow(n_samples),
      last_observation = t,
    )

}

get_IC_forecast <- function(t, observations, forecast_ens, prev_effort){

    N <- forecast_ens |>
      mutate(N = abundance_after_take) |>
      filter(PPNum == t,
             horizon == min(horizon))

    if(prev_effort){
      N <- N |>
        select(prev_effort_scenario, effort_scenario, N, ens) |>
        group_by(effort_scenario, prev_effort_scenario)
    } else {
      N <- N |>
        select(effort_scenario, abundance_after_take, N, ens) |>
        group_by(effort_scenario)
    }

    last_observation <- observations$PPNum[max(which(observations$PPNum < t))]

    N |>
      mutate(id = cur_group_id(),
             last_observation = last_observation) |>
      ungroup()
}

sim_N <- function(IC, start_date, fx_date, phi, zeta){

  new_abundance <- IC |>
    group_by(id) |>
    mutate(
      start_date = start_date,
      fx_date = fx_date,
      survived = rbinom(length(N), N, phi),
      recruited = rpois(length(N), N * zeta),
      abundance = survived + recruited) |>
    ungroup() |>
    select(-N)

}

remove_extinct <- function(df){
  not_extinct <- df |>
    group_by(prev_effort_scenario, effort_scenario) |>
    summarise(sum_ens = sum(abundance_after_take)) |>
    ungroup() |>
    suppressMessages()

  if(any(not_extinct$sum_ens == 0)){
    not_extinct <- not_extinct |> filter(sum_ens > 0)
    pes <- pull(not_extinct, prev_effort_scenario)
    es <- pull(not_extinct, effort_scenario)
    df_keep <- tibble()
    for(e in 1:nrow(not_extinct)){
      dfe <- df |>
        filter(prev_effort_scenario == pes[e]) |>
        filter(effort_scenario == es[e])
      df_keep <- bind_rows(df_keep, dfe)
    }
    return(df_keep)
  } else {
    return(df)
  }
}

abundance_m_take <- function(t, obs, new_abundance, y_ens_store){
  # prep new_abundance for joining with y_ens_store
  if(is.na(new_abundance$effort_scenario[1])) {
    new_abundance <- new_abundance |>
      select(-effort_scenario, -id)
    relationship <- "many-to-one"
  } else {
    if("prev_effort_scenario" %in% colnames(new_abundance)){
      new_abundance <- new_abundance |>
        select(-id, -prev_effort_scenario) |>
        rename(prev_effort_scenario = effort_scenario)
    } else {
      new_abundance <- new_abundance |>
        select(-id) |>
        rename(prev_effort_scenario = effort_scenario)
    }
    relationship <- "many-to-many"
  }

  abundance_take_ens <- y_ens_store |>
    filter(start_pp == t) |>
    group_by(effort_scenario, ens, horizon, start_pp, PPNum, property) |>
    summarise(total_take = sum(take)) |>
    ungroup() |>
    left_join(new_abundance, relationship = relationship) |>
    mutate(abundance_after_take = pmax(0, abundance - total_take)) |>
    suppressMessages()

  if(obs) abundance_take_ens <- abundance_take_ens |> mutate(prev_effort_scenario = NA)
  abundance_take_ens |>
    select(c("property", "start_date", "fx_date",
             "horizon", "start_pp", "PPNum", "last_observation",
             "survived", "recruited", "abundance",
             "total_take", "abundance_after_take",
             "prev_effort_scenario", "effort_scenario", "ens"))
}
