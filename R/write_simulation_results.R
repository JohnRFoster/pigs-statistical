library(tidyverse)
library(lubridate)
library(nimble)

model_dir <- "modifiedDM_betaSurvival_dataByMethod"
sim_dir <- file.path("out/simulation", model_dir)

sim_runs <- list.files(sim_dir)

all_samples <- tibble()
all_take <- tibble()
all_N <- tibble()
all_M <- tibble()
all_beta_p <- tibble()
all_pop_growth <- tibble()
all_phi <- tibble()
all_properties <- tibble()
all_methods <- tibble()
all_y <- tibble()
all_area <- tibble()
property_lookup <- tibble()

first <- TRUE
# sim_runs <- sim_runs[1:2]

if(first){
  prev_tasks <- "none"
} else {
  prev_tasks <- read_csv("out/simulation/KnownMethodParameters.csv") |>
    pull(simulation) |>
    unique()
}

pb <- txtProgressBar(max = length(sim_runs), style = 3)

for(i in seq_along(sim_runs)){
  # print(i)
  run_experiments <- list.files(file.path(sim_dir, sim_runs[i]))
  for(j in seq_along(run_experiments)){

    path <- file.path(sim_dir, sim_runs[i], run_experiments[j])
    simulation_files <- list.files(path)

    rds <- read_rds(file.path(path, "posteriorEval.rds"))
    bad_mcmc <- rds$bad_mcmc
    task_id <- rds$task_id
    already_collated <- task_id %in% prev_tasks

    if(bad_mcmc | already_collated) next

    samples <- read_csv(file.path(path, "posteriorSamples.csv")) |>
      mutate(simulation = task_id) |>
      suppressMessages()

    all_samples <- bind_rows(all_samples, samples)

    y_pred <- read_csv(file.path(path, "posteriorPredictedTake.csv")) |>
      suppressMessages()
    colnames(y_pred) <- 1:ncol(y_pred)
    y_pred <- y_pred  |>
      mutate(simulation = task_id)
    all_y <- bind_rows(all_y, y_pred)

    pot_area <- read_csv(file.path(path, "posteriorPredictedPotentialArea.csv")) |>
      suppressMessages()
    colnames(pot_area) <- 1:ncol(pot_area)
    pot_area <- as_tibble(pot_area) |>
      mutate(simulation = task_id)
    all_area <- bind_rows(all_area, pot_area)

    sim_data <- read_rds(file.path(path, "sim_data.rds"))

    take <- sim_data$take |>
      mutate(simulation = task_id,
             p_id = 1:n())
    all_take <- bind_rows(all_take, take)

    N <- sim_data$N |>
      mutate(simulation = task_id)
    all_N <- bind_rows(all_N, N)

    M <- sim_data$M |>
      mutate(simulation = task_id)
    all_M <- bind_rows(all_M, M)

    county_to_property <- sim_data$county_to_property |>
      mutate(simulation = task_id)

    property_lookup <- bind_rows(property_lookup, county_to_property)

    # need a lookup table for known data model covariates
    bH <- tibble(
      method_idx = rep(1:nrow(sim_data$beta_p), ncol(sim_data$beta_p)),
      position = rep(1:ncol(sim_data$beta_p), each = nrow(sim_data$beta_p)),
      actual = as.numeric(sim_data$beta_p),
      simulation = task_id
    )

    all_beta_p <- bind_rows(all_beta_p, bH)

    popH <- as_tibble(sim_data$pop_growth) |>
      mutate(simulation = task_id)
    all_pop_growth <- bind_rows(all_pop_growth, popH)

    # table for which phi matches to each property and PP
    pH <- sim_data$constants$pH |>
      as_tibble() |>
      mutate(property = 1:n()) |>
      pivot_longer(cols = -property,
                   names_to = "timestep",
                   values_to = "idx") |>
      filter(!is.na(idx)) |>
      mutate(simulation = task_id)

    # table for which primary periods are observed at each property
    ph_1 <- sim_data$constants$PPNum |>
      as_tibble() |>
      mutate(property = 1:n()) |>
      pivot_longer(cols = -property,
                   names_to = "timestep",
                   values_to = "PPNum") |>
      filter(!is.na(PPNum)) |>
      group_by(property) |>
      filter(PPNum < max(PPNum)) |>
      ungroup() |>
      mutate(simulation = task_id,
             observation = 1) |>
      select(-timestep)

    # the known phi value for each PP in each property
    ph_2 <- sim_data$phi |>
      as_tibble() |>
      mutate(property = 1:n()) |>
      pivot_longer(cols = -property,
                   names_to = "PPNum",
                   values_to = "phi") |>
      mutate(PPNum = as.integer(PPNum)) |>
      filter(!is.na(phi)) |>
      mutate(simulation = task_id)

    # all pp between first and last sample for each property
    ph_3 <- sim_data$constants$all_pp |>
      as_tibble() |>
      mutate(property = 1:n()) |>
      pivot_longer(cols = -property,
                   names_to = "timestep",
                   values_to = "PPNum") |>
      filter(!is.na(PPNum)) |>
      mutate(simulation = task_id) |>
      group_by(property) |>
      filter(PPNum < max(PPNum)) |>
      ungroup()

    phj_1 <- left_join(ph_3, ph_2) |>
      suppressMessages()
    phj_2 <- left_join(phj_1, pH) |>
      suppressMessages()
    phj_3 <- left_join(phj_2, ph_1) |>
      mutate(observation = if_else(is.na(observation), 0, observation)) |>
      suppressMessages()
    all_phi <- bind_rows(all_phi, phj_3)

    method_lookup <- sim_data$method_lookup |>
      mutate(simulation = task_id)

    all_methods <- bind_rows(all_methods, method_lookup)


    property_attributes <- read_rds(file.path(path, "property_attributes.rds"))
    pa <- property_attributes |>
      select(-phi_mu, -psi_phi, -initial_abundnace) |>
      mutate(simulation = task_id)

    all_properties <- bind_rows(all_properties, pa)

  }
  setTxtProgressBar(pb, i)
}
close(pb)

analysis_dir <- "analysis/simulation"

path <- file.path(analysis_dir, model_dir)
if(!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

ammend_prev <- function(df, name, first, file_name){

  dest <- file.path(path, file_name)

  if(first){
    write_csv(df, dest)
  } else {
    prev <- read_csv(dest)

    out <- bind_rows(prev, df)
    assign(name, out, envir = .GlobalEnv)
    write_csv(out, dest)
  }
  message(name, " done\n")
}

ammend_prev(all_samples, "all_samples", first, "Samples.csv")
ammend_prev(all_y, "all_y", first, "PredictedTake.csv")
ammend_prev(all_take, "all_take", first, "KnownTake.csv")
ammend_prev(all_N, "all_N", first, "KnownPropertyAbundance.csv")
ammend_prev(all_M, "all_M", first, "KnownCountyAbundance.csv")
ammend_prev(all_beta_p, "all_beta_p", first, "KnownDataCovariates.csv")
ammend_prev(all_pop_growth, "all_pop_growth", first, "KnownPopulationGrowth.csv")
ammend_prev(all_methods, "all_methods", first, "KnownMethodParameters.csv")
ammend_prev(all_properties, "all_properties", first, "Properties.csv")
ammend_prev(property_lookup, "property_lookup", first, "propertyLookup.csv")

X <- all_take |>
  select(simulation, county, property, area_property, contains("c_")) |>
  distinct()

select_pivot_longer <- function(df, node){
  df |>
    select(contains(node), simulation) |>
    pivot_longer(cols = -c(simulation),
                 names_to = "node")
}

beta1_long <- all_samples |>
  select_pivot_longer("beta1") |>
  mutate(method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
         position = 1)

recovered <- function(df){
  df |> mutate(parameter_recovered = if_else(actual >= low & actual <= high, 1, 0))
}

my_summary <- function(df){
  df |>
    summarise(low = quantile(value, 0.025),
              med = quantile(value, 0.5),
              high = quantile(value, 0.975))
}

## capture probability intercepts ------
beta1_recovery <- beta1_long |>
  group_by(simulation, node, method_idx, position) |>
  my_summary() |>
  left_join(all_beta_p) |>
  ungroup() |>
  recovered()

out_dir <- "out/simulation"

table(beta1_recovery$parameter_recovered)
sum(beta1_recovery$parameter_recovered) / nrow(beta1_recovery)

beta1_residual <- beta1_long |>
  left_join(all_beta_p) |>
  mutate(value = value - actual) |>
  group_by(node, position, method_idx) |>
  my_summary() |>
  ungroup()

write_csv(beta1_recovery, file.path(out_dir, "summaryRecovered_captureIntercept.csv"))
write_csv(beta1_residual, file.path(out_dir, "summaryResidual_captureIntercept.csv"))
message("capture intercepts done")

## capture probability covariates ------
beta_p_long <- all_samples |>
  select_pivot_longer("beta_p") |>
  mutate(method_idx = as.numeric(str_extract(node, "(?<=\\[)\\d")),
         position = as.numeric(str_extract(node, "(?<=\\, )\\d")) + 1)

beta_p_recovery <- beta_p_long |>
  group_by(simulation, node, method_idx, position) |>
  my_summary() |>
  left_join(all_beta_p) |>
  ungroup() |>
  recovered()

beta_p_residual <- beta_p_long |>
  left_join(all_beta_p) |>
  mutate(value = value - actual) |>
  group_by(node, position, method_idx) |>
  my_summary() |>
  ungroup()

write_csv(beta_p_recovery, file.path(out_dir, "summaryRecovered_captureCovariates.csv"))
write_csv(beta_p_residual, file.path(out_dir, "summaryResidual_captureCovariates.csv"))
message("capture covariates done")

## gamma ------

gH <- all_methods |>
  select(idx, gamma, method, simulation) |>
  filter(method %in% c("SNARE", "TRAPS")) |>
  mutate(idx = idx - 3) |>
  rename(actual = gamma)

gamma_long <- all_samples |>
  select_pivot_longer("log_gamma[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d"))) |>
  mutate(value = exp(value))

gamma_recovery <- gamma_long |>
  group_by(simulation, node, idx) |>
  my_summary() |>
  left_join(gH) |>
  ungroup() |>
  recovered()

gamma_residual <- gamma_long |>
  left_join(gH)|>
  mutate(value = value - actual) |>
  group_by(node, idx) |>
  my_summary() |>
  ungroup()

write_csv(gamma_recovery, file.path(out_dir, "summaryRecovered_saturationConstant.csv"))
write_csv(gamma_residual, file.path(out_dir, "summaryResidual_saturationConstant.csv"))
message("saturation constant done")

## rho ------
rH <- all_methods |>
  select(idx, rho, method, simulation) |>
  rename(actual = rho)

rho_long<- all_samples |>
  select_pivot_longer("log_rho[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d"))) |>
  mutate(value = exp(value))

rho_recovery <- rho_long |>
  group_by(simulation, node, idx) |>
  my_summary() |>
  left_join(rH) |>
  ungroup() |>
  recovered()

rho_residual <- rho_long |>
  left_join(rH)|>
  mutate(value = value - actual) |>
  group_by(node, idx) |>
  my_summary() |>
  ungroup()

write_csv(rho_recovery, file.path(out_dir, "summaryRecovered_searchArea.csv"))
write_csv(rho_residual, file.path(out_dir, "summaryResidual_searchArea.csv"))
message("search area done")

## unique area ------
pH <- all_methods |>
  select(idx, p_unique, method, simulation) |>
  rename(actual = p_unique)

p_mu_long <- all_samples |>
  select_pivot_longer("p_mu[") |>
  mutate(idx = as.numeric(str_extract(node, "(?<=\\[)\\d"))) |>
  mutate(value = ilogit(value))

p_mu_recovery <- p_mu_long |>
  group_by(simulation, node, idx) |>
  my_summary() |>
  left_join(pH) |>
  ungroup() |>
  recovered()

p_mu_residual <- p_mu_long |>
  left_join(pH)|>
  mutate(value = value - actual) |>
  group_by(node, idx) |>
  my_summary() |>
  ungroup()

write_csv(p_mu_recovery, file.path(out_dir, "summaryRecovered_uniqueArea.csv"))
write_csv(p_mu_residual, file.path(out_dir, "summaryResidual_uniqueArea.csv"))
message("unique area done")

## litter size ------
actual <- 5.290323
ls_long <- all_samples |>
  select_pivot_longer("log_mean_ls") |>
  mutate(value = exp(value))

ls_recovery <- ls_long |>
  group_by(simulation, node) |>
  my_summary() |>
  mutate(actual = actual) |>
  ungroup() |>
  recovered()

ls_residual <- ls_long |>
  mutate(actual = actual) |>
  mutate(value = value - actual) |>
  group_by(node) |>
  my_summary() |>
  ungroup()

write_csv(ls_recovery, file.path(out_dir, "summaryRecovered_litterSize.csv"))
write_csv(ls_residual, file.path(out_dir, "summaryResidual_litterSize.csv"))
message("liter size done")

## survival ------
actual <- 0.75
phi_long <- all_samples |>
  select_pivot_longer("phi_mu")

phi_recovery <- phi_long |>
  group_by(simulation, node) |>
  my_summary() |>
  mutate(actual = actual) |>
  ungroup() |>
  recovered()

phi_residual <- phi_long |>
  mutate(actual = actual) |>
  mutate(value = value - actual) |>
  group_by(node) |>
  my_summary() |>
  ungroup()

actual <- 3
psi_phi_long <- all_samples |>
  select_pivot_longer("psi_phi")

psi_phi_recovery <- psi_phi_long |>
  group_by(simulation, node) |>
  my_summary() |>
  mutate(actual = actual) |>
  ungroup() |>
  recovered()

psi_phi_residual <- psi_phi_long |>
  mutate(actual = actual) |>
  mutate(value = value - actual) |>
  group_by(node) |>
  my_summary() |>
  ungroup()

write_csv(bind_rows(phi_recovery, psi_phi_recovery), file.path(out_dir, "summaryRecovered_survial.csv"))
write_csv(bind_rows(phi_residual, psi_phi_residual), file.path(out_dir, "summaryResidual_survial.csv"))
message("survival done")


## abundance ------
abundance <- left_join(property_lookup, all_N) |>
  select(n_id, county, property, PPNum, abundance, density, simulation, area_property)

xn <- all_samples |>
  select_pivot_longer("xn[") |>
  filter(!is.na(value)) |>
  mutate(n_id = as.numeric(str_extract(node, "(?<=\\[)\\d*"))) |>
  left_join(abundance) |>
  filter(!is.na(abundance))

xn_posterior <- xn |>
  mutate(estimated_density = value / area_property) |>
  group_by(node, n_id, county, property, PPNum, abundance, simulation, area_property, density) |>
  summarise(low_abundance = quantile(value, 0.025),
            med_abundance = quantile(value, 0.5),
            high_abundance = quantile(value, 0.975),
            low_density = quantile(estimated_density, 0.025),
            med_density = quantile(estimated_density, 0.5),
            high_density = quantile(estimated_density, 0.975)) |>
  ungroup()

write_csv(xn_posterior, file.path(out_dir, "summaryPropertyAbundance.csv"))
message("posterior abundance done")

xn_error <- xn |>
  mutate(estimated_density = value / area_property) |>
  group_by(node, n_id, county, property, PPNum, abundance, simulation, area_property, density) |>
  summarise(mae_abundance = mean(abs(value - abundance)),
            mae_density = mean(abs(estimated_density - density)),
            mpe_abundance = mean(abs((value+1) - (abundance+1))/(abundance+1))*100,
            mpe_density = mean(abs((estimated_density+0.1) - (density+0.1))/(density+0.1))*100,
            mbias_abundance = mean(value - abundance),
            mbias_density = mean(estimated_density - density),
            mse_abundance = mean((value - abundance)^2),
            mse_density = mean((estimated_density - density)^2),
            rmse_abundance = (sqrt(mse_abundance)),
            rmse_density = (sqrt(mse_density))) |>
  ungroup() |>
  arrange(simulation, property, PPNum) |>
  group_by(simulation, property) |>
  mutate(delta = PPNum - lag(PPNum)) |>
  ungroup()

write_csv(xn_error, file.path(out_dir, "PropertyAbundanceMetrics.csv"))
message("abundance metrics done")

posterior_take <- all_y |>
  pivot_longer(cols = -c(simulation),
               names_to = "p_id") |>
  filter(!is.na(value)) |>
  mutate(p_id = as.numeric(p_id)) |>
  group_by(simulation, p_id) |>
  my_summary() |>
  ungroup() |>
  left_join(all_take)

write_csv(posterior_take, file.path(out_dir, "summaryPredictedTake.csv"))
message("posterior take done")

