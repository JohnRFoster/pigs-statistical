library(targets)
library(tidyverse)
library(lubridate)
library(config)
library(spdep)
library(readxl)
library(spatialreg)
library(parallel)
library(nimble)
library(nimbleHMC)
library(rgdal)

setwd("C:/Users/John.Foster/OneDrive - USDA/Desktop/fosteR/pigs-statistical/")

run_parallel <- TRUE
rep_num <- 1
# set.seed(rep_num)
out_dir <- "out/simulation"
model_dir <- "modifiedDM_betaSurvival_dataByMethod"
likelihood <- "poisson"
mcmc_config <- "customMCMC"
rep <- paste0("simulation_", rep_num)

dest <- file.path(out_dir, model_dir, likelihood, mcmc_config, rep)
if(!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

message("Writing samples to: ", dest)

phi_mu <- 0.75
psi_phi <- 3
sigma_dem <- 0.25

a_phi <- phi_mu * psi_phi
b_phi <- (1 - phi_mu) * psi_phi
par(mfrow = c(1, 1))
hist(rbeta(10000, a_phi, b_phi))

property_attributes <- expand_grid(
  property_density = round(runif(5, 0.3, 6), 2),
  proportion_county_surveyed = seq(0.1, 0.9, length.out = 5)
  # survival_rate = seq(0.6, 0.95, length.out = 4)
) |>
  mutate(
    county = as.numeric(as.factor(proportion_county_surveyed)),
    # survival_rate = round(ilogit(rnorm(n(), logit_mean_phi, sigma_phi)), 3),
    phi_mu = phi_mu,
    psi_phi = psi_phi,
    area_property = round(runif(n(), 5, 250), 2),
    initial_abundnace = round(property_density * area_property),
    n_sample_occasions = round(runif(n(), 2.5, 16.4))
  ) |>
  group_by(county) |>
  mutate(
    area_county_surveyed = sum(area_property),
    area_county = area_county_surveyed / proportion_county_surveyed,
    area_not_sampled = area_county - area_county_surveyed
  ) |>
  ungroup() |>
  arrange(county) |>
  mutate(property = 1:n())

# view(property_attributes)

demographic_stochasticity <- TRUE
n_county <- max(property_attributes$county)
n_property <- max(property_attributes$property)
n_pp <- 20
c_road_den <- rnorm(n_county)
c_rugged <- rnorm(n_county)
c_canopy <- rnorm(n_county)
beta_p <- matrix(rnorm(20), 5, 4)

source("R/functions_nimble.R")
source("R/functions_simulation.R")

sim_data <- simulate_dm(
  property_data = property_attributes,
  likelihood = likelihood,
  phi_mu = phi_mu,
  psi_phi = psi_phi,
  sigma_dem = sigma_dem,
  n_pp = n_pp,
  c_road_den = c_road_den,
  c_rugged = c_rugged,
  c_canopy = c_canopy,
  beta_p = beta_p,
  demographic_stochasticity = demographic_stochasticity,
  plot = TRUE
)

source("R/priors.R")
phi_prior <- create_surv_prior()

constants <- sim_data$constants
# constants$phi_mu_a <- 1
# constants$phi_mu_b <- 1
constants$phi_mu_a <- phi_prior$alpha
constants$phi_mu_b <- phi_prior$beta
data <- sim_data$data

likelihood_nb <- if_else(likelihood == "nb", TRUE, FALSE)
likelihood_binom <- ifelse(likelihood == "binomial", TRUE, FALSE)
likelihood_poisson <- ifelse(likelihood == "poisson", TRUE, FALSE)
spatial <- FALSE

post_dir <- file.path(out_dir, "modifiedDM_betaSurvival_uninformative", likelihood, mcmc_config, "simulation_1")

message("==== Test build ====")
source("R/dm_inits.R")
source("R/nimble_dm_2.R")
Rmodel <- nimbleModel(
  code = modelCode,
  constants = constants,
  data = data,
  inits = inits(data, constants),
  calculate = TRUE
)
print(warnings())
message("\n\n")


monitors_add <- c("xn", "p", "log_theta")
params_check <- c(
  "beta_p",
  "beta1",
  "log_gamma",
  "log_rho",
  "phi_mu",
  "psi_phi",
  "log_mean_ls",
  "p_mu"
)
model_flags <- list(likelihood = likelihood,
                    spatial = spatial,
                    demographic_stochasticity = demographic_stochasticity)

custom_samplers <- tribble(
  ~node,            ~type,
  "log_mean_ls",    "slice",
  "phi_mu",         "slice",
  "psi_phi",        "slice",
  "log_rho",        "AF_slice"
)
n_iter <- 10000
n_chains <- 1
max_iter <- 80000

sim_data$params_check <- params_check
sim_data$custom_samplers <- custom_samplers
write_rds(sim_data,
          file = file.path(dest, "sim_data.rds"))
write_rds(property_attributes,
          file = file.path(dest, "property_attributes.rds"))

if(run_parallel){


  message("==== Build cluster ===")
  cl <- makeCluster(3)
  source("R/run_nimble_parallel.R")
  samples <- run_nimble_parallel(
    cl = cl,
    model_code = modelCode,
    model_data = data,
    model_constants = constants,
    model_inits = inits,
    inits_dir = NULL,
    n_iter = n_iter,
    params_check = params_check,
    monitors_add = monitors_add,
    state.col = "xn",
    model_flags = model_flags,
    custom_samplers = custom_samplers,
    effective_size = 5000,
    max_iter = max_iter,
    max_psrf = 50,
    calculate = TRUE,
    use_conjugacy = TRUE,
    resetMV = TRUE,
    save_iter = TRUE,
    dest = dest
  )
  stopCluster(cl)
  message("Stop cluster")


} else {

  Rmodel <- nimbleModel(
    code = modelCode,
    constants = constants,
    data = data,
    inits = inits(data, constants),
    calculate = TRUE
  )
  # warnings()

  # Rmodel$calculate()
  # Rmodel$simulate(Rmodel$getDependencies('logit_phi'))

  Rmodel$initializeInfo()

  # default MCMC configuration
  mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)
  mcmcConf$addMonitors(monitors_add)

  if(!is.null(custom_samplers)){
    for(i in seq_len(nrow(custom_samplers))){
      node <- custom_samplers$node[i]
      type <- custom_samplers$type[i]
      mcmcConf$removeSampler(node)
      mcmcConf$addSampler(node, type)
    }
  }

  for(i in 1:5){
    node <- paste0("beta_p[", i, ", ", 1:constants$m_p, "]")
    node <- c(paste0("beta1[", i, "]"), node)
    mcmcConf$removeSampler(node)
    mcmcConf$addSampler(node, "AF_slice")
  }

  mcmcConf$printSamplers(byType = TRUE)

  Rmcmc <- buildMCMC(mcmcConf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc)

  n_iter <- 5
  n_chains <- 1

  samples <- runMCMC(
    Cmcmc,
    niter = n_iter,
    nchains = n_chains,
    # thin = 5,
    samplesAsCodaMCMC = TRUE
  )

  params <- split_out(samples, "xn")
  params2 <- split_out(params$params, "pop_growth")
  params3 <- split_out(params2$params, "phi")
  params4 <- split_out(params3$params, "log_lambda_1")
  plot(params4$params)

  write_rds(
    samples,
    file = file.path(dest, "allSamples.rds")
  )
}



