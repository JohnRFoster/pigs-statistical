
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

#'@description Calculate the potential search area from posterior samples when shooting methods (fixed wing, helicopter, hunting) are used
#'@param log_rho vector of mcmc samples for rho (log scale)
#'@param p_unique vector of mcmc samples for p
#'@param effort_per the effort per trap/snare
#'@param n_trap_m1 the number of traps/snares used minus 1

shooting_lpa <- function(log_rho, p_unique, effort_per, n_trap_m1){
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

calc_log_potential_area <- function(method, log_rho, log_gamma, p_unique, effort_per, n_trap_m1){
  n_mcmc <- length(log_rho)
  n_reps <- length(method)
  log_potential_area <- matrix(NA, n_mcmc, n_reps)
  for(j in seq_along(method)){
    lr <- log_rho[,method[j]]
    p <- p_unique[,method[j]]
    ep <- effort_per[j]
    tcm1 <- n_trap_m1[j]
    if(method[j] == 4 | method[j] == 5){
      lg <- log_gamma[,method[j]-3]
      log_potential_area[,j] <- trap_snare_lpa(lr, lg, p, ep, tcm1)
    } else {
      log_potential_area[,j] <- shooting_lpa(lr, p, ep, tcm1)
    }
  }
  return(log_potential_area)
}

#'@description calculate the capture probability for a removal event
#'@param X matrix of county-level land cover covariates, first column is 1's for an intercept term
#'@param beta vector of coeffiencents
#'@param log_potential_area matrix of potential search area, iterations (rows) by replicate (columns)
#'@param area_property scalar value of the area (mk^2) of the property

calc_p <- function(X, beta, log_potential_area, area_property){
  n_mcmc <- nrow(beta_p)
  n_reps <- ncol(log_potential_area)
  p <- matrix(NA, n_mcmc, n_reps)
  for(mc in seq_len(n_mcmc)){
    log_theta <- numeric(n_reps)
    for(j in seq_len(n_reps)){

      # base probability of capture given an individual is the first survey
      log_theta[j] <- log(ilogit(inprod(X, beta[mc,]))) +
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
    beta_p <- samples[, grep("beta_p", colnames(samples))]

    for(i in 1:n_survey){
      if(method[i] == 4 | method[i] == 5){
        log_potential_area[,i] <- trap_snare_lpa(
          log_rho = log_rho[,method[i]],
          log_gamma = log_gamma[,method[i]-3],
          p_unique = p_unique[,method[i]],
          effort_per = effort_per[i],
          n_trap_m1 = n_trap_m1[i]
        )
      } else {
        log_potential_area[,i] <- shooting_lpa(
          log_rho = log_rho[,method[i]],
          p_unique = p_unique[,method[i]],
          effort_per = effort_per[i],
          n_trap_m1 = n_trap_m1[i]
        )
      }

      if(calc_p){

        pb <- txtProgressBar(max = nrow(samples), style = 3)
        for(m in 1:nrow(samples)){
          M <- samples[m,]

          # base probability of capture given an individual is the first survey
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











