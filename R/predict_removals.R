library(targets)
library(tidyverse)
library(lubridate)
library(config)
library(spdep)
library(nimble)
setwd("C:/Users/John.Foster/OneDrive - USDA/Desktop/fosteR/pigs-statistical")
source("R/functions_predict.R")
source("R/functions_nimble.R")

n_fx <- 1
n_reps <- 3

removal_effort <- seq(0, 1, by = 0.25)

rep_num <- 2
out_dir <- "out/simulation"
model_dir <- "DM_recruitData_varyingEffort"
likelihood <- "poisson"
mcmc_config <- "customMCMC_conjugate"
rep <- paste0("simulation_", rep_num)

dest <- file.path(out_dir, model_dir, likelihood, mcmc_config, paste0("simulation_", rep_num))
rds <- read_rds(file.path(dest, "posterior.rds"))
samples <- rds$samples
sim_data <- read_rds(file.path(dest, "sim_data.rds"))
attach(sim_data$constants)
attach(sim_data$data)

sample_occasions <- tibble(
  property = property_x,
  PPNum = pp_x
) |>
  group_by(property, PPNum) |>
  mutate(timestep = cur_group_id()) |>
  ungroup() |>
  mutate(n_id = 1:n())

property_df <- tibble(
  y = sim_data$data$y,
  property = p_property_idx,
  sampled_pp = p_pp_idx,
  start = start,
  end = end,
  method = method,
  county = p_county_idx
) |>
  group_by(property, sampled_pp) |>
  mutate(timestep = cur_group_id()) |>
  ungroup() |>
  mutate(n_id = 1:n())

take <- sim_data$take

draws <- sample.int(nrow(samples), 500, replace = TRUE)

n_samples <- samples[draws, grep("xn[", colnames(samples), fixed = TRUE)]
n_mcmc <- nrow(n_samples)

pattern <- "(?<!beta_)p\\[\\d*\\]"
p_detect <- str_detect(colnames(samples), pattern)
p <- samples[,which(p_detect)]

mu_phi <- samples[,"logit_mean_phi"]
sigma_phi <- sigma_2_tau(samples[,"tau_phi"])
mean_ls <- exp(samples[,"log_mean_ls"])
zeta <- mean_ls*28*1/365

log_rho <- samples[, grep("log_rho", colnames(samples))]
log_gamma <- samples[, grep("log_gamma", colnames(samples))]
p_unique <- ilogit(samples[, grep("p_mu", colnames(samples))])
beta_p <- samples[, grep("beta_p", colnames(samples))]
beta1 <- samples[, grep("beta1", colnames(samples))]

X <- unique(X_p)




# determine the initial condition
# if the starting primary period has an observation, use the posterior for that PP
# otherwise, use the latest forecast


# get a specific date
get_date <- function(pp_filter){
  time_lookup |>
    filter(pp == pp_filter) |>
    pull(date) |>
    as.character()
}

effort_quants <- take |>
  group_by(property, method, area_property) |>
  reframe(
    effort_per = quantile(effort_per, removal_effort), q = removal_effort,
    trap_count = quantile(trap_count, removal_effort), q = removal_effort) |>
  ungroup() |>
  suppressMessages()

n_reps <- take |>
  select(property, PPNum) |>
  group_by(property, PPNum) |>
  tally() |>
  group_by(property) |>
  reframe(n_reps = quantile(n, removal_effort), q = removal_effort) |>
  mutate(n_reps = round(n_reps))

#


i <- 1
for(i in 1:n_property){

  message("Forecasting property ", i, "/", n_property)

  observations <- sample_occasions |>
    filter(property == i)

  start_pp <- min(observations$PPNum)
  end_pp <- max(observations$PPNum)
  fx_pp_seq <- start_pp:(end_pp + n_fx)

  fx_date_seq <- seq.Date(
    from = ymd("2023-09-01"),
    length.out = length(fx_pp_seq),
    by = "4 week")

  # need a timestep to PP/date lookup table for saving
  time_lookup <- tibble(
    timestep = seq_along(fx_pp_seq),
    pp = fx_pp_seq,
    date = fx_date_seq
  )

  effort_property <- effort_quants |>
    filter(property == i)

  reps_property <- n_reps |>
    filter(property == i)

  c <- take |>
    filter(property == i) |>
    slice(1) |>
    pull(county)

  methods <- take |>
    filter(property == i) |>
    pull(method)

  pb <- txtProgressBar(min = start_pp, max = end_pp, style = 3)

  for(t in start_pp:end_pp){

    # if the forecast issue date has an observation, use the abundance posterior
    # otherwise, use the latest forecast
    # need to keep track of effort scenarios to propagate how switching between
    # effort levels affects abundance
    obs <- t %in% observations$PPNum
    if(obs){
      IC <- get_IC_posterior(t, observations, n_samples)
    } else {
      prev_effort <- "prev_effort_scenario" %in% colnames(new_abundance)
      IC <- get_IC_forecast(t, observations, forecast_ens, prev_effort)
    }

    phi <- ilogit(rnorm(n_mcmc, mu_phi, sigma_phi)) # draw survival
    start_date <- get_date(t)
    fx_date <- get_date(t+1)

    # forecast new abundance
    new_abundance <- sim_N(IC, start_date, fx_date, phi, zeta)

    # removal effort scenarios
    y_ens_store <- tibble()
    groups <- unique(new_abundance$id)
    for(j in seq_along(groups)){
      for(k in seq_along(removal_effort)){

        y_attr <- remove(i, t, j, removal_effort[k], effort_property, reps_property, methods,
                         log_rho, log_gamma, p_unique, X[c,], beta1, beta_p, new_abundance)
        y_ens_store <- bind_rows(y_ens_store, y_attr)

      } # k
    } # j

    forecast_ens <- abundance_m_take(t, obs, new_abundance, y_ens_store)

    fx_dir <- file.path(dest, "forecast", paste0("property_", i))
    if(!dir.exists(fx_dir)) dir.create(fx_dir, recursive = TRUE)

    fx_file <- file.path(fx_dir, paste0(start_date, "_abundance.csv"))
    write_csv(forecast_ens, fx_file)

    y_file <- file.path(fx_dir, paste0(start_date, "_takeByMethod.csv"))
    write_csv(y_ens_store, y_file)

    setTxtProgressBar(pb, t)

  } # t
  close(pb)
}


