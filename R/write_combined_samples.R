library(targets)
library(tidyverse)
library(lubridate)
library(nimble)
library(coda)
library(parallel)
setwd("C:/Users/John.Foster/OneDrive - USDA/Desktop/fosteR/pigs-statistical/")
source("R/functions_nimble.R")
source("R/functions_predict.R")


overide_existing_samples <- FALSE

rep_num <- 7
out_dir <- "out/simulation"
model_dir <- "modifiedDM_betaSurvival_Informative"
likelihood <- "poisson"
mcmc_config <- "customMCMC"
rep <- paste0("simulation_", rep_num)

dest <- file.path(out_dir, model_dir, likelihood, mcmc_config, rep)
message("Checking MCMC for ", file.path(model_dir, likelihood, mcmc_config, rep))

# out_dir <- "out/data"
# model_dir <- "modifiedDM_recruitData_varyingEffort_captureByMethod"
# likelihood <- "poisson"
# mcmc_config <- "customMCMC"
# dest <- file.path(out_dir, model_dir, likelihood, mcmc_config)


mm <- get_mcmc_chunks(dest, start = 1)

all_nodes <- colnames(mm$params[[1]])

params_check <- c(
  "beta_p",
  "beta1",
  "log_gamma",
  "log_rho",
  "phi_mu",
  "psi_phi",
  # "tau_dem",
  "log_mean_ls",
  "p_mu"
)

j <- unlist(lapply(params_check, function(x) grep(x, all_nodes)))

params <- mm$params[,j]

message("Creating traceplots...")
png(filename = file.path(dest, "mcmcTimeseries%03d.png"))
plot(params)
dev.off()
message("  done")

message("Calculating PSRF...")
psrf <- gelman.diag(params, multivariate = FALSE)
print(psrf)

message("Calculating effective samples...")
effective_samples <- effectiveSize(params)
print(effective_samples)

message("Calculating burnin...")
ff <- tempfile()
png(filename = ff)
GBR <- gelman.plot(params)
dev.off()
unlink(ff)

burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
print(burnin)

if(is.na(burnin) | burnin >= 0.9*nrow(params[[1]])){
  burnin <- 60000
}

# burnin <- 100

mcmc_1 <- as.matrix(window(mm$params, start = burnin))
mcmc_2 <- as.matrix(window(mm$predict, start = burnin))

draws <- sample.int(nrow(mcmc_1), 5000, replace = TRUE)

samples <- cbind(mcmc_1[draws,], mcmc_2[draws,])

if(grepl("simulation", out_dir)){
  sim_data <- read_rds(file.path(dest, "sim_data.rds"))
} else {
  sim_data <- read_rds(file.path(dest, "nimbleList.rds"))
}

constants <- sim_data$constants
data <- sim_data$data

post <- data_posteriors(samples, constants, data)

write_rds(
  list(samples = samples,
       y_pred = post$y,
       potential_area = post$potential_area,
       # theta = post$theta,
       psrf = psrf,
       effective_samples = effective_samples,
       burnin = burnin),
  file = file.path(dest, "posterior.rds")
)

message("DONE!")
