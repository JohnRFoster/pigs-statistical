## functions to simulate pig population abundance

simulate_dm <- function(
    property_data,
    likelihood,
    phi_mu,
    psi_phi,
    sigma_dem,
    n_pp,
    c_road_den,
    c_rugged,
    c_canopy,
    beta_p,
    demographic_stochasticity,
    plot = FALSE
){

  n_county <- max(property_data$county)
  n_property <- max(property_data$property)

  # mean litters per year from VerCauteren et al. 2019 pg 64
  data_litters_per_year <- c(1, 2, 0.86, 1, 2.28, 2.9, 0.49, 0.85, 1.57)

  # mean litter size year from VerCauteren et al. 2019 pg 63
  data_litter_size <- c(5.6, 6.1, 5.6, 6.1, 4.2, 5.0, 5.0, 6.5, 5.5, 6.8,
                        5.6, 5.9, 4.9, 5.1, 4.5, 4.7, 5.3, 5.7, 7.4, 8.4,
                        4.7, 4.9, 3.0, 3.0, 4.8, 4.8, 4.2, 5.4, 4.7, 5.2, 5.4)

  data_litter_size <- round(data_litter_size)

  # ecological process
  zeta <- mean(data_litter_size)*28*1/365 # mean pigs produced per female per 28 days

  county <- property_data$county
  n_units <- n_property * n_pp
  survey_area <- property_data$area_property
  log_survey_area <- log(survey_area)

  unit <- matrix(1:n_units, n_property, n_pp, byrow = TRUE)

  # methods
  method_lookup <- tibble(
    idx = 1:5,
    method = c("FIREARMS", "FIXED WING", "HELICOPTER", "SNARE", "TRAPS"),
    freq = c(0.1, 0.05, 0.05, 0.3, 0.5),
    p_unique = runif(5),
    rho = c(runif(1, 0.01, 5), # firearms; p_mu[1]
            runif(1, 5, 20),   # fixed wing
            runif(1, 5, 20),   # helicopter
            runif(1, 0.5, 2),  # snare; gamma[1], p_mu[2]
            runif(1, 0.5, 2)), # traps; gamma[2], p_mu[3]
    gamma = c(0.5, 0.01, 0.01, runif(1, 0.5, 4), runif(1, 0.5, 4))
  )

  log_rho <- log(method_lookup$rho)
  log_gamma <- log(method_lookup$gamma)
  p_unique <- method_lookup$p_unique

  a_phi <- phi_mu * psi_phi
  b_phi <- (1 - phi_mu) * psi_phi

  # method <- sample.int(5, 1, prob = method_lookup$freq)
  method <- sample.int(5, 1)

  generate_take_data <- function(method){
    if(method == 1){
      effort_per <- rpois(1, 5) + 1
      trap_count <- rpois(1, 2) + 1
    } else if(method == 2){
      effort_per <- max(0.1, rnorm(1, 1.65, 1.1))
      trap_count <- 1
    } else if(method == 3){
      effort_per <- max(0.01, rexp(1, 0.75))
      trap_count <- sample.int(4, 1, prob = c(0.935, 0.049, 0.014, 0.002))
    } else if(method == 4){
      effort_per <- max(1, round(rexp(1, 1)))
      trap_count <- max(1, round(rnorm(1, 16, 10)))
    } else if(method == 5){
      effort_per <- max(1, round(rnorm(1, 1, 8)))
      trap_count <- max(1, round(rnorm(1, 1, 8)))
    }
    tibble(effort_per = effort_per, trap_count = trap_count)
  }

  # storage
  N <- matrix(NA, n_property, n_pp)
  phi <- matrix(NA, n_property, n_pp-1)
  pop_growth <- matrix(0, n_property, n_pp-1)
  E <- tibble()

  # initial abundance for each property
  N[, 1] <- property_data$initial_abundnace

  # for tracking calculations
  log_area_check <- p_check <- n_passes_track <- 0

  for(i in seq_len(n_property)){

    # determine which primary periods are sampled for each property
    sample_occ <- sample(n_pp, property_data$n_sample_occasions[i])

    for(t in seq_len(n_pp)){

      if(t %in% sample_occ){

        n_passes <- round(runif(1, 1.5, 5.4))
        n_passes_track <- c(n_passes_track, n_passes)
        log_theta <- C <- rep(NA, n_passes)
        for(j in seq_len(n_passes)){

          # determine take method and effort
          m <- sample.int(5, 1)
          e <- generate_take_data(m)
          effort_per <- pull(e, effort_per)
          log_effort_per <- log(effort_per)
          trap_count <- pull(e, trap_count)
          n_trap_m1 <- trap_count - 1

          if(m == 1){ # firearms
            log_potential_area <- log_rho[m] +
              log_effort_per -
              log(1 + (p_unique[m] * n_trap_m1))
          } else if(m == 2 | m ==3){ # fixed wing and helicopter
            log_potential_area <- log_rho[m] + log_effort_per
          } else if(m == 4 | m == 5){
            log_potential_area <- log(pi) +
              (2 * (log_rho[m] + log_effort_per -
                      log(exp(log_gamma[m-3]) + effort_per))) +
              log(1 + (p_unique[m-2] * n_trap_m1))
          }

          log_pr_area_sampled <- min(log_survey_area[i], log_potential_area)
          log_area_check <- c(log_area_check, log_survey_area[i] - log_potential_area)

          # probability of capture, given that an individual is in the surveyed area
          log_theta[j] <- log(
            boot::inv.logit(
              beta_p[m, 1] +
              beta_p[m, 2] * c_road_den[county[i]] +
              beta_p[m, 3] * c_rugged[county[i]] +
              beta_p[m, 4] * c_canopy[county[i]])) +
          min(0, log_potential_area - log_survey_area[i])

          if(j == 1){
            p <- exp(log_theta[1])
            n_avail <- N[i, t]
          } else {
            p <- exp(log_theta[1] + sum(log(1 - exp(log_theta[1:(j-1)]))))
            n_avail <- N[i, t] - sum(C[1:(j-1)])
          }
          p_check <- c(p_check, p)

          if(likelihood == "binomial"){
            C[j] <- min(rbinom(1, n_avail, p), n_avail)
          } else if(likelihood == "nb"){
            C[j] <- min(rnbinom(1, mu = n_avail * p, size = size[county[i]]), n_avail)
          } else if(likelihood == "poisson"){
            C[j] <- min(rpois(1, n_avail * p), n_avail)
          }

          et <- tibble(
            property = i,
            PPNum = t,
            take = C[j],
            method = m,
            effort_per = effort_per,
            trap_count = trap_count
          )
          E <- bind_rows(E, et)

        }
        z <- N[i, t] #- sum(C)
      } else {
        z <- N[i, t]
      }

      if(t < n_pp){
        # phi[i, t] <- boot::inv.logit(rnorm(1, logit_mean_phi, sigma_phi))
        phi[i, t] <- rbeta(1, a_phi, b_phi)
        if(demographic_stochasticity){
          # S <- rbinom(1, z, phi[i, t])  # survival
          # R <- rpois(1, zeta/2*z) # recruitment

          # z <- max(0.00001, z)
          # log_lam <- log(z) + log(zeta/2 + phi[i, t])
          lam <- z * zeta / 2 + z * phi[i, t]
          # Ndem <- max(0, rnorm(1, lam, sigma_dem))

        } else {
          S <- round(z * phi[i, t])
          R <- round(zeta * z / 2)
        }
        # N[i, t+1] <- S + R
        N[i, t+1] <- rpois(1, lam)
        if(N[i, t] != 0){
          pop_growth[i, t] <- N[i, t+1] / N[i, t]
        }
      }
    }
  }

  if(plot){
    par(mfrow = c(1,2))
    areas <- property_data$area_property
    plot(N[1,]/areas[1], type = "l", ylim = c(0, 10))
    for(i in 2:n_property) lines(N[i,]/areas[i], col = i)
    plot(N[1,], type = "l", ylim = c(0, max(N)))
    for(i in 2:n_property) lines(N[i,], col = i)
  }

  p_check <- p_check[-1]
  colnames(N) <- 1:n_pp

  nH <- property_data |>
    select(-initial_abundnace, -n_sample_occasions)

  N_long <- N |>
    as_tibble() |>
    mutate(property = 1:n(),
           county = county) |>
    pivot_longer(cols = c(-property, -county),
                 names_to = "PPNum",
                 values_to = "abundance") |>
    left_join(nH) |>
    mutate(PPNum = as.numeric(PPNum),
           density = abundance / area_property)

  M_long <- N_long |>
    group_by(county, PPNum, area_county, area_county_surveyed, proportion_county_surveyed, area_not_sampled) |>
    summarise(total_abundance = sum(abundance)) |>
    ungroup() |>
    mutate(scaled_county_abundance = total_abundance / area_county_surveyed * area_county,
           county_abundance_low = NA,
           county_abundance_high = NA)

  for(i in 1:nrow(M_long)){
    county_abundance <- pmax(
      rpois(10000, M_long$scaled_county_abundance[i]),
      M_long$total_abundance[i])
    M_long$county_abundance_low[i] <- quantile(county_abundance, 0.025)
    M_long$county_abundance_high[i] <- quantile(county_abundance, 0.975)
  }

  m_start <- N_long |>
    filter(PPNum == 1) |>
    group_by(county, area_not_sampled) |>
    summarise(avg_density = mean(density)) |>
    ungroup() |>
    mutate(abundance = round(area_not_sampled * avg_density))

  colnames(phi) <- 1:ncol(phi)
  phi_long <- phi |>
    as_tibble() |>
    mutate(property = 1:n(),
           county = county) |>
    pivot_longer(cols = c(-property, -county),
                 names_to = "PPNum",
                 values_to = "phi")

  N_unsampled <- matrix(NA, n_county, n_pp)
  N_unsampled[,1] <- m_start$abundance
  for(i in 1:n_county){
    for(t in 2:n_pp){
      phi_mu <- phi_long |>
        filter(county == i,
               PPNum == t-1) |>
        pull(phi) |>
        mean()
      N_unsampled[i, t] <- rbinom(1, N_unsampled[i, t-1], phi_mu) +
        rpois(1, N_unsampled[i, t-1] / 2 * zeta)
    }
  }

  colnames(N_unsampled) <- 1:n_pp
  Nu_long <- N_unsampled |>
    as_tibble() |>
    mutate(property = (1:n())+100,
           county = 1:n()) |>
    pivot_longer(cols = c(-property, -county),
                 names_to = "PPNum",
                 values_to = "abundance") |>
    mutate(PPNum = as.numeric(PPNum))

  M <- N_long |>
    select(property, county, PPNum, abundance) |>
    bind_rows(Nu_long) |>
    group_by(county, PPNum) |>
    summarise(latent_county_abundance = sum(abundance)) |>
    right_join(M_long)

  cH <- property_data |>
    select(county, property) |>
    group_by(county) |>
    mutate(Property_county = 1:n()) |>
    ungroup()

  areas <- property_data |>
    select(property, area_property)

  take <- left_join(E, areas) |>
    left_join(cH) |>
    group_by(property, PPNum) |>
    mutate(order = 1:n()) |>
    ungroup()

  timestep_df <- take |>
    select(property, PPNum) |>
    distinct() |>
    group_by(property) |>
    mutate(timestep = 1:n()) |>
    ungroup()

  take <- left_join(take, timestep_df)

  sum_area_surveyed <- take |>
    group_by(county, PPNum) |>
    distinct() |>
    summarise(sum_area_surveyed = sum(area_property))

  sampled_units <- take |>
    select(county, property, PPNum, timestep, Property_county) |>
    distinct() |>
    mutate(n_id = 1:n())

  n_timesteps <- property_data$n_sample_occasions


  county_areas <- property_data |>
    select(county, area_county) |>
    distinct()

  county_sampled_units <- sampled_units |>
    select(county, PPNum) |>
    distinct() |>
    mutate(m_id = 1:n()) |>
    left_join(county_areas)

  sampled_units <- left_join(sampled_units, county_sampled_units)

  timestep <- sampled_units |>
    select(-n_id, -m_id, -Property_county, -area_county) |>
    pivot_wider(names_from = timestep,
                values_from = PPNum) |>
    select(-county, -property)

  all_pp <- tibble()
  for(i in 1:n_property){
    sub <- filter(sampled_units, property == i)
    reps <- tibble(
      property = i,
      PPNum = min(sub$PPNum):max(sub$PPNum),
      timestep = PPNum
    ) |>
      mutate(timestep = timestep - min(timestep) + 1)
    all_pp <- bind_rows(all_pp, reps)
  }

  all_pp_wide_prop <- all_pp |>
    pivot_wider(names_from = timestep,
                values_from = PPNum)

  pop_growth_lookup <- all_pp_wide_prop |>
    pivot_longer(cols = -property,
                 names_to = "timestep",
                 values_to = "PPNum") |>
    filter(!is.na(PPNum)) |>
    group_by(property) |>
    filter(PPNum < max(PPNum)) |>
    ungroup() |>
    select(-PPNum) |>
    mutate(H = 1:n()) |>
    pivot_wider(values_from = H,
                names_from = timestep) |>
    select(-property)

  all_pp_wide <- all_pp_wide_prop |>
    select(-property)

  n_pp_include <- apply(all_pp_wide, 1, function(x) max(which(!is.na(x))))

  X_lookup <- tibble(
    county = 1:n_county,
    c_road_den = c_road_den,
    c_rugged = c_rugged,
    c_canopy = c_canopy
  )

  take <- left_join(take, X_lookup)

  X_p <- take |>
    select(starts_with("c_")) |>
    as.matrix()
  # X_p <- cbind(rep(1, nrow(X_p)), X_p)

  # Generate start and end indices for previous surveys ---------------------
  take$start <- 0
  take$end <- 0

  pb <- txtProgressBar(max = nrow(take), style = 3)
  for (i in 1:nrow(take)) {
    if (take$order[i] > 1) {
      idx <- which(take$county == take$county[i] &
                     take$property == take$property[i] &
                     take$timestep == take$timestep[i] &
                     take$order < take$order[i])
      take$start[i] <- idx[1]
      take$end[i] <- idx[length(idx)]
      tar_assert_identical(idx, take$start[i]:take$end[i])
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  y_sum <- take |>
    group_by(property, PPNum) |>
    mutate(ysum = cumsum(take) - take) |>
    ungroup() |>
    select(property, PPNum, ysum, order)

  y_rem <- take |>
    group_by(property, PPNum) |>
    summarise(ysum = sum(take)) |>
    ungroup()

  N_init <- take |>
    select(property, timestep, take) |>
    mutate(p = p_check,
           N = round((take+1) / p)) |>
    filter(timestep == 1) |>
    group_by(property) |>
    summarise(N = sum(N)) |>
    ungroup() |>
    pull(N)

  y_sum_wide <- left_join(all_pp, y_rem) |>
    mutate(ysum = if_else(is.na(ysum), 0, ysum)) |>
    select(-PPNum) |>
    pivot_wider(names_from = timestep,
                values_from = ysum) |>
    select(-property)

  county_units <- county_sampled_units

  M_join <- sampled_units |>
    select(county, PPNum, n_id, Property_county) |>
    left_join(county_units)

  M_wide <- M_join |>
    pivot_wider(id_cols = m_id, # long index for a county x pp unit
                values_from = n_id, # long index for a property x pp unit
                names_from = Property_county) |>
    select(-m_id) |>
    as.matrix()

  M_lookup <- matrix(NA, nrow(M_wide), ncol(M_wide))
  for(i in 1:nrow(M_wide)){
    vec <- which(!is.na(M_wide[i,]))
    M_lookup[i,1:length(vec)] <- M_wide[i, vec]
  }
  n_prop <- apply(M_lookup, 1, function(x) max(which(!is.na(x))))
  county_to_property <- left_join(sampled_units, county_units)


  constants <- list(
    n_timesteps = n_timesteps,
    n_county = n_county,
    n_survey = nrow(take),
    n_lpy = length(data_litters_per_year),
    n_ls = length(data_litter_size),
    n_property = n_property,
    n_county_units = nrow(county_units),
    n_prop = n_prop,
    n_first_survey = length(which(take$order == 1)),
    n_not_first_survey = length(which(take$order != 1)),
    n_trap_snare = length(which(take$method %in% 4:5)),
    n_shooting = length(which(take$method %in% 1:3)),
    n_units = nrow(sampled_units),
    n_method = 5,
    n_pp = max(n_pp_include),
    n_pp_prop = n_pp_include,
    all_pp = as.matrix(all_pp_wide),
    p_county_idx = take$county,
    n_county_idx = sampled_units$county,
    log_pi = log(pi),
    m_p = ncol(X_p),
    first_survey = which(take$order == 1),
    not_first_survey = which(take$order != 1),
    p_property_idx = take$property,
    p_pp_idx = take$PPNum,
    start = take$start,
    end = take$end,
    PPNum = as.matrix(timestep),
    method = take$method,
    ts_idx = which(take$method %in% 4:5),
    shooting_idx = which(take$method %in% 1:3),
    property_x = sampled_units$property,
    pp_x = sampled_units$PPNum,
    pp_len = 28,
    phi_prior_mean = 3.3,
    phi_prior_tau = 1,
    M_lookup = M_lookup,
    county = M_long$county,
    N_init = N_init,
    pH = as.matrix(pop_growth_lookup)
  )

  data <- list(
    y = take$take,
    K = data_litters_per_year,
    J = data_litter_size,
    y_sum = y_sum$ysum,
    rem = as.matrix(y_sum_wide),
    sum_prop_area = sum_area_surveyed$sum_area_surveyed,
    county_area = county_sampled_units$area_county,
    property_area = property_data$area_property,
    X_p = X_p,
    effort_per = take$effort_per,
    log_effort_per = log(take$effort_per),
    n_trap_m1 = take$trap_count - 1,
    log_survey_area_km2 = log(take$area_property)
  )

  return(
    list(
      constants = constants,
      data = data,
      take = take,
      N = N_long,
      M = M,
      beta_p = beta_p,
      method_lookup = method_lookup,
      log_pr_area_sampled = log_pr_area_sampled,
      county_to_property = county_to_property,
      pop_growth = pop_growth,
      phi = phi,
      p_check = p_check
    )
  )
}






