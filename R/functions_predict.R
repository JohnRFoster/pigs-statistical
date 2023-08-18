



calc_log_potential_area_trap_snare <- function(log_rho, log_gamma, p_unique, effort_per, n_trap_m1){
  log(pi) +
    (2 * (log_rho + log(effort_per) -
            log(exp(log_gamma) + effort_per))) +
    log(1 + (p_unique * n_trap_m1))
}
calc_log_potential_area_shooting <- function(log_rho, p_unique, effort_per, n_trap_m1){
  log_rho +
    log(effort_per) -
    log(1 + (p_unique * n_trap_m1))
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
        log_potential_area[,i] <- calc_log_potential_area_trap_snare(
          log_rho = log_rho[,method[i]],
          log_gamma = log_gamma[,method[i]-3],
          p_unique = p_unique[,method[i]],
          effort_per = effort_per[i],
          n_trap_m1 = n_trap_m1[i]
        )
      } else {
        log_potential_area[,i] <- calc_log_potential_area_shooting(
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











