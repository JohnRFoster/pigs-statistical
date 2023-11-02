library(parallel)
library(targets)
library(tidyverse)
library(lubridate)
library(nimble)


n_nodes <- 25
n_iter <- 80000
n_chains <- 3

n_pp <- 20
phi_mu <- 0.78
psi_phi <- 6
sigma_dem <- 0.25

inits_dir <- NULL
out_dir <- "out/simulation"
model_dir <- "modifiedDM_betaSurvival_dataByMethod"
sim_dir <- file.path(out_dir, model_dir)
past_reps <- list.files(sim_dir) |> as.numeric()
if(length(past_reps) == 0){
  rep_num <- 1
} else {
  rep_num <- max(past_reps) + 1
}

message("\nBegin simulations for replicate ", rep_num)

dest <- file.path(out_dir, model_dir, rep_num)

message("Writing simulations to: ", dest)

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

custom_samplers <- tribble(
  ~node,            ~type,
  "log_mean_ls",    "slice",
  "phi_mu",         "slice",
  "psi_phi",        "slice",
  "log_rho",        "AF_slice"
)

source("R/nimble_dm_2.R")
source("R/server_fit_simulations_parallel.R")

message("\nStart Cluster")
cl <- makeCluster(n_nodes)

message("=== Begin Simulation ===")
system.time(
  fit_simulations_parallel(
    cl = cl,
    rep_num = rep_num,
    phi_mu = phi_mu,
    psi_phi = psi_phi,
    sigma_dem = sigma_dem,
    n_pp = n_pp,
    modelCode = modelCode,
    n_iter = n_iter,
    n_chains = n_chains,
    params_check = params_check,
    monitors_add = monitors_add,
    custom_samplers = custom_samplers,
    inits_dir = inits_dir,
    dest = dest
  )
)

stopCluster(cl)

message("Simulations complete!")
message("Combing output and writing CSVs...")
source("R/write_simulation_results.R")
message("\nDONE!")
