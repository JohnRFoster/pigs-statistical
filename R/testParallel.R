library(tidyverse)
library(targets)
library(nimble)
library(parallel)

source("R/nimble_saturating.R")
source("R/functions_nimble.R")
# st_d <- tar_read("st_d")
# group_d_dev <- tar_read("group_d_dev")
# group_d_dev$group_d |>
#   select(method, method_factor) |>
#   mutate(
#          method = as.numeric(method_factor),
#          trap_snare_indicator = as.numeric(method_factor %in% c("trap", "snare")),
#          trap_snare_idx = pmax(method - 3, 1)) |>
#   distinct()

nimls <- tar_read("nimble_lists_train")
str(nimls)

constants <- nimls$constants
data <- nimls$data

log_rho_priors_inf <- data.frame(
  muLog = c(-1.654013, 2.74274, 2.74274, 0.2777449, 0.6438721),
  sigmaLog = c(0.8092821, 0.5976205, 0.5976205, 0.2255966, 0.2225168)
)

log_rho_priors_unInf <- data.frame(
  mu = rep(0, 5),
  sigma = rep(1, 10)
)

gamma_priors_inf <- data.frame(
  shape = c(7.704547, 3.613148),
  rate = c(4.41925, 3.507449)
)

gamma_priors_unInf <- data.frame(
  shape = rep(0.1, 2),
  rate = rep(0.1, 2)
)

data$log_rho_prior <- log_rho_priors_unInf
data$gamma_prior <- gamma_priors_unInf
data$area_constraint <- rep(1, constants$n_survey)

inits <- function(){

  ls <- list(
    lambda = rpois(length(data$y), data$y) + 1,
    gamma = rgamma(5, 1, 1),
    rho = rgamma(5, 1, 1),
    beta = rnorm(constants$m_n, 0, 1),
    beta_p = rnorm(constants$m_p, 0, 1),
    p_unique = rbeta(constants$n_method, 1, 1),
    eps_propertyR = rnorm(constants$n_property, 0, 1),
    eps_property_pR = rnorm(constants$n_property, 0, 1),
    tau_car = runif(1, 10, 90),
    sigma_property = rexp_nimble(1, 1),
    sigma_property_p = rexp_nimble(1, 1),
    sigma_st0 = rexp_nimble(1, 1),
    sigma_st = rexp_nimble(1, 1),
    sigma_short = rexp_nimble(1, 1),
    eta = rbeta(1, 9, 1),
    eps_stR = matrix(rnorm(constants$n_timestep*constants$n_county, 0, 1), constants$n_timestep, constants$n_county),
    z_shortR = matrix(0, constants$m_short, constants$n_county)
  )

  return(ls)
}


params_check <- c(
  "beta",
  "gamma",
  "rho",
  "p_unique",
 # "tau",
  "sigma"
)

message("Build cluster")
cl <- makeCluster(3)
system.time(
  samples <- run_nimble_parallel(
    cl = cl,
    model_code = modelCode,
    model_data = data,
    model_constants = constants,
    model_inits = inits,
    n_iter = 5000,
    n_burnin = 0,
    params_check = params_check,
    max_iter = 200000,
    dest = "out/MO/dcar_precision",
    resetMV = TRUE,
    save_iter = TRUE
  )
)
stopCluster(cl)
message("Stop cluster")

# model_code = modelCode
# model_data = data
# model_constants = constants
# model_inits = inits
# init = model_inits()
# n_iter = 1000
# n_burnin = 0
# params_check = params_check
# max_iter = 5000
# resetMV = TRUE
# save_iter = TRUE
# effective_size = 5000
# max_psrf = 500
# calculate = TRUE
# use_conjugacy = FALSE


str(samples)

mcmc <- samples$samples
plot(mcmc[,"sigma_property"])


