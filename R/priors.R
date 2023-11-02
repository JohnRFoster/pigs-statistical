## -------------------
# create informative priors for survival and recruitment of pigs
## -------------------


library(tidyverse)
library(lubridate)
library(nimble)

create_surv_prior <- function(logit = TRUE, sd_inflate = 1){
  data <- read_csv("data/insitu/Vital_Rate_Data.csv")

  data_usa <- data |>
    filter(country == "USA",
           time.period.end != "null",
           time.period.start != "null",
           !paper.ID %in% c(128, 1007, 130, 136)) |> # these papers don't have specified date ranges or are meta-analysis
    mutate(time.period.end = mdy(time.period.end),
           time.period.start = mdy(time.period.start))

  surv_data <- data_usa |>
    filter(!is.na(survival.prop)) |>
    select(
      unique.ID,
      paper.ID,
      N.hogs.in.study,
      contains("survival"),
      contains("hunting"),
      state,
      contains("time"),
      method.for.data)

  surv_mu <- surv_data |>
    mutate(weeks = as.numeric(time.period.end - time.period.start)/7,
           weeks4 = weeks / 4,
           survival.per.4week = survival.prop ^ (1/weeks4),
           logit.survival.per.4week = logit(survival.per.4week)) |>
    filter(survival.per.4week > 0) |>
    mutate(scale_factor = survival.per.4week/survival.prop)

  surv_mu_summary <- surv_mu |>
    summarise(mu = mean(survival.per.4week),
              mu.logit = mean(logit.survival.per.4week))

  surv_var <- surv_data |>
    filter(survival.var.type %in% c("SD", "95% CI"))

  surv_sd <- surv_var |>
    filter(survival.var.type == "SD") |>
    mutate(sd = as.numeric(survival.var))

  surv_sd_calc <- surv_var |>
    filter(survival.var.type == "95% CI") |>
    mutate(low.CI = as.numeric(str_extract(survival.var, "[[:graph:]]*(?=\\-)")),
           high.CI = as.numeric(str_extract(survival.var, "(?<=\\-)[[:graph:]]*")),
           sd_low = (low.CI - survival.prop) / -1.96,
           sd_high = (high.CI - survival.prop) / 1.96) |>
    group_by(unique.ID) |>
    summarise(sd = max(sd_high, sd_low))

  surv_var_join <- left_join(surv_var, surv_sd_calc) |>
    filter(survival.var.type != "SD")

  scale_ids <- surv_mu |>
    select(unique.ID, scale_factor)

  surv_variance <- bind_rows(surv_var_join, surv_sd) |>
    left_join(scale_ids) |>
    mutate(variance = sd^2,
           variance.4week = variance * scale_factor^2,
           sd.4week = sqrt(variance.4week))

  surv_sd_summary <- surv_variance |>
    pull(sd.4week) |>
    mean()

  mu <- surv_mu_summary$mu
  psi <- 1 / mean(surv_variance$variance.4week)
  alpha <- mu * psi
  beta <- (1 - mu) * psi

  # mu_phi <- mean(surv_mu$survival.per.4week)
  # var_phi <- var(surv_mu$survival.per.4week)
  #
  # k_phi <- mu_phi / (var_phi / mu_phi)
  # theta_phi <- var_phi / mu_phi
  #
  # k_theta <- function(mu, var){
  #   k <- mu / (var / mu)
  #   theta <- var / mu
  #   c(k, theta)
  # }
  #
  # g_mu <- k_theta(mean(surv_mu$survival.per.4week), var(surv_mu$survival.per.4week))
  # g_var <- k_theta(mean(surv_variance$variance.4week), var(surv_variance$variance.4week))
  #
  # n <- 100000
  # alpha <- rgamma(n, g_mu[1], 1/g_mu[2])
  # beta <- rgamma(n, g_var[1], 1/g_var[2])
  # phi <- numeric(n)
  # for(i in 1:n){
  #   phi[i] <- rbeta(1, alpha[i], beta[i])
  # }
  # hist(phi)

  return(list(
    surv_mu = ifelse(logit, surv_mu_summary$mu.logit, surv_mu_summary$mu),
    surv_sd = surv_sd_summary * sd_inflate,
    alpha = alpha,
    beta = beta
  ))
}


# repro_data <- data_usa |>
#   filter(!is.na(survival.prop)) |>
#   select(
#     unique.ID,
#     paper.ID,
#     N.hogs.in.study,
#     contains("pregnancy"),
#     contains("hunting"),
#     contains("litter"),
#     contains("fecundity"),
#     contains("age"),
#     demographics,
#     gender,
#     state,
#     contains("time"),
#     method.for.data)

