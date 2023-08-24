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
rep_num <- 70 # starting with five using inv. gamma for initial state prior
# set.seed(rep_num)
out_dir <- "out/simulation"
model_dir <- "DM_recruitData_varyingEffort"
likelihood <- "poisson"
mcmc_config <- "customMCMC_conjugate"
rep <- paste0("simulation_", rep_num)

dest <- file.path(out_dir, model_dir, likelihood, mcmc_config, rep)
if(!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

logit_mean_phi <- 1.8
sigma_phi <- 1.2

property_attributes <- expand_grid(
  property_density = round(runif(10, 0.3, 6), 2),
  proportion_county_surveyed = seq(0.1, 0.9, length.out = 5)
  # survival_rate = seq(0.6, 0.95, length.out = 4)
) |>
  mutate(
    county = as.numeric(as.factor(proportion_county_surveyed)),
    # survival_rate = round(ilogit(rnorm(n(), logit_mean_phi, sigma_phi)), 3),
    logit_mean_phi = logit_mean_phi,
    sigma_phi = sigma_phi,
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
beta_p <- rnorm(4)

source("R/functions_nimble.R")
source("R/functions_simulation.R")

sim_data <- simulate_dm(
  property_data = property_attributes,
  likelihood = likelihood,
  logit_mean_phi = logit_mean_phi,
  sigma_phi = sigma_phi,
  n_pp = n_pp,
  c_road_den = c_road_den,
  c_rugged = c_rugged,
  c_canopy = c_canopy,
  beta_p = beta_p,
  demographic_stochasticity = demographic_stochasticity
)

constants <- sim_data$constants
data <- sim_data$data

likelihood_nb <- if_else(likelihood == "nb", TRUE, FALSE)
likelihood_binom <- ifelse(likelihood == "binomial", TRUE, FALSE)
likelihood_poisson <- ifelse(likelihood == "poisson", TRUE, FALSE)
spatial <- FALSE

inits <- make_inits_function_dm(
  inits_dir = NULL,
  constants = constants,
  data = data,
  demographic_stochasticity = demographic_stochasticity
)

message("==== Test build ====")
source("R/nimble_dm_2.R")
Rmodel <- nimbleModel(
  code = modelCode,
  constants = constants,
  data = data,
  inits = inits(),
  calculate = TRUE
)
message("\n\n")


monitors_add <- c("xn", "lambda", "p", "log_theta")
params_check <- c(
  "beta_p",
  "log_gamma",
  "log_rho",
  "logit_mean_phi",
  "sigma_phi",
  "mean_ls",
  "p_mu"
)
model_flags <- list(likelihood = likelihood,
                    spatial = spatial,
                    demographic_stochasticity = demographic_stochasticity)

custom_samplers <- tibble(
  node = c("log_gamma", "log_rho", "beta_p"),
  type = c("AF_slice", "AF_slice", "AF_slice")
)

n_iter <- 50000
n_chains <- 1
max_iter <- 500000


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

  inits <- make_inits_function_dm(
    inits_dir = NULL,
    constants = constants,
    data = data,
    demographic_stochasticity = demographic_stochasticity
  )


  Rmodel <- nimbleModel(
    code = modelCode,
    constants = constants,
    data = data,
    inits = inits(),
    calculate = TRUE
  )
  # warnings()

  # Rmodel$calculate()
  # Rmodel$simulate(Rmodel$getDependencies('logit_phi'))

  Rmodel$initializeInfo()

  library(compareMCMCs)
  nimbleMCMCdefs <- list(

    default = function(model){
      mcmcConf <- configureMCMC(model)
      mcmcConf
    },

    RW = function(model){
      mcmcConf <- configureMCMC(model, onlyRW = TRUE)
      mcmcConf
    }

    # AFblock_data = function(model){
    #   mcmcConf <- configureMCMC(model)
    #   mcmcConf$removeSamplers(c('log_gamma', 'log_rho', 'p_mu'))
    #   mcmcConf$addSampler(target = c("log_rho[1]", "p_mu[1]"), type = "AF_slice")
    #   mcmcConf$addSampler(target = c("log_rho[2]", "p_mu[2]"), type = "AF_slice")
    #   mcmcConf$addSampler(target = c("log_rho[3]", "p_mu[3]"), type = "AF_slice")
    #   mcmcConf$addSampler(target = c("log_rho[4]", "p_mu[4]", "log_gamma[1]"), type = "AF_slice")
    #   mcmcConf$addSampler(target = c("log_rho[5]", "p_mu[5]", "log_gamma[2]"), type = "AF_slice")
    #   mcmcConf
    # },
    #
    # RWblock_data = function(model){
    #   mcmcConf <- configureMCMC(model)
    #   mcmcConf$removeSamplers(c('log_gamma', 'log_rho', 'beta_p', 'p_mu'))
    #   mcmcConf$addSampler(target = c("log_rho[1]", "p_mu[1]"), type = "RW_block")
    #   mcmcConf$addSampler(target = c("log_rho[2]", "p_mu[2]"), type = "RW_block")
    #   mcmcConf$addSampler(target = c("log_rho[3]", "p_mu[3]"), type = "RW_block")
    #   mcmcConf$addSampler(target = c("log_rho[4]", "p_mu[4]", "log_gamma[1]"), type = "RW_block")
    #   mcmcConf$addSampler(target = c("log_rho[5]", "p_mu[5]", "log_gamma[2]"), type = "RW_block")
    #   mcmcConf
    # }

  )

  init_vals <- inits()

  mcmcResults_nimble <- compareMCMCs(
    modelInfo = list(code = modelCode,
                     data = data,
                     constants = constants, # centered
                     inits = init_vals),
    # Omit monitors argument to
    # use default monitors: top-level parameters
    MCMCs = c('default',
              'RW'), # Ditto
    nimbleMCMCdefs = nimbleMCMCdefs,
    MCMCcontrol = list(niter = 4000, burnin = 100)
  )

  # MCMC efficiency for a parameter is defined as the effective sample size divided by computation time in seconds
  # It is the number of effectively independent samples generated per second.
  res <- combineMetrics(mcmcResults_nimble)
  res$byMCMC
  res$byParameter |>
    as_tibble() |>
    filter(Parameter == "log_mean_ls")

  mcmcResults_nimble[['default']]$metrics$byParameter
  # mcmcResults_nimble[['AFslice_data']]$metrics$byParameter
  mcmcResults_nimble[['AFblock_data']]$metrics$byParameter
  mcmcResults_nimble[['RWblock_data']]$metrics$byParameter
  make_MCMC_comparison_pages(mcmcResults_nimble, modelName = "mis_dm")

  # default MCMC configuration
  mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)
  mcmcConf$addMonitors(monitors_add)

  # custom_samplers <- tibble(
  #   node = c("beta_p", "log_gamma"),
  #   type = c("AF_slice", "AF_slice")
  # )

  if(!is.null(custom_samplers)){
    for(i in seq_len(nrow(custom_samplers))){
      node <- custom_samplers$node[i]
      type <- custom_samplers$type[i]
      mcmcConf$removeSampler(node)
      mcmcConf$addSampler(node, type)
    }
  }

  # cv <- runCrossValidate(mcmcConf, 20, lossFunction = "predictive")

  mcmcConf$printSamplers(byType = TRUE)

  Rmcmc <- buildMCMC(mcmcConf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc)

  n_iter <- 20000
  n_chains <- 3

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



