library(targets)
library(tidyverse)
library(lubridate)
library(config)
library(spdep)
library(readxl)
library(spatialreg)
library(parallel)
library(nimble)
library(rgdal)

setwd("C:/Users/John.Foster/OneDrive - USDA/Desktop/fosteR/pigs-statistical/")
source("R/functions_nimble.R")

run_parallel <- TRUE
# set.seed(rep_num)
out_dir <- "out/data"
model_dir <- "modifiedDM_recruitData_varyingEffort_captureByMethod"
likelihood <- "poisson"
mcmc_config <- "customMCMC"

dest <- file.path(out_dir, model_dir, likelihood, mcmc_config)
if(!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

data <- read_csv("data/insitu/sample_equal.csv")

# data <- data |>
#   arrange(agrp_prp_id, start.date, timestep_4) |>
#   group_by(agrp_prp_id, timestep_4) |>
#   mutate(cumsum_take = cumsum(take)) |>
#   filter(!all(cumsum_take == 0)) |>
#   ungroup() |>
#   group_by(agrp_prp_id, timestep_4) |>
#   mutate(two_plus_takes = n() >= 2) |>
#   filter(two_plus_takes) |>
#   group_by(agrp_prp_id) |>
#   mutate(n_timesteps = length(unique(timestep_4))) |>
#   filter(n_timesteps >= 2) |>
#   arrange(cnty_name, agrp_prp_id, timestep_4) |>
#   ungroup()

property_filter <- data |>
  # filter(method != "SNARE") |>
  group_by(agrp_prp_id) |>
  summarise(n = sum(take)) |>
  filter(n > 0) |>
  pull(agrp_prp_id)

tar_assert_equal_lengths(property_filter, unique(data$agrp_prp_id))

take_df <- data |>
  # filter(method != "SNARE") |>
  select(-alws_agrprop_id, -property.size, -two_plus_takes, -n_timesteps, -p, -`...1`) |>
  rename(Method = method,
         Property = agrp_prp_id,
         PPNum = timestep_4) |>
  filter(Property %in% property_filter) |>
  mutate(County = sprintf("%05d", fips)) |>
  mutate(Method = if_else(grepl("TRAPS", Method), "TRAPS", Method),
         County = as.character(County),
         method_idx = as.numeric(as.factor(Method)),
         property_idx = as.numeric(as.factor(Property)),
         county_idx = as.numeric(as.factor(County))) |>
  arrange(Property, start.date, PPNum)
  # group_by(Property, PPNum) |>
  # mutate(cumsum_take = cumsum(take)) |>
  # filter(!all(cumsum_take == 0)) |>
  # ungroup()

property_timestep <- take_df |>
  select(Property, PPNum) |>
  distinct() |>
  group_by(Property) |>
  mutate(timestep = 1:n()) |>
  ungroup()

take_df <- left_join(take_df, property_timestep)

# property_filter_2 <- take_df |>
#   select(Property, PPNum) |>
#   distinct() |>
#   group_by(Property) |>
#   count() |>
#   filter(n >= 5) |>
#   pull(Property)
#
# take_df <- take_df |>
#   filter(Property %in% property_filter_2) |>
#   mutate(property_idx = as.numeric(as.factor(Property)))

all_county_units <- take_df |>
  select(Property, County, county_idx, PPNum, property_area_km2) |>
  distinct() |>
  group_by(Property) |>
  mutate(timestep = 1:n()) |>
  ungroup() |>
  group_by(County, Property) |>
  mutate(cnty_property = cur_group_id()) |>
  ungroup()

sum_prop_area <- all_county_units |>
  group_by(County, PPNum) |>
  summarise(sum_area = sum(property_area_km2)) |>
  ungroup()

n_prop_county <- all_county_units |>
  mutate(n_idx = 1:n()) |>
  select(County, PPNum, cnty_property, n_idx) |>
  pivot_wider(names_from = cnty_property,
              values_from = n_idx,
              id_cols = c(County, PPNum)) |>
  arrange(County, PPNum) |>
  select(-County, -PPNum) |>
  as.matrix()

n_prop <- apply(n_prop_county, 1, function(x) length(which(!is.na(x))))

n_prop_county1 <- matrix(NA, nrow(n_prop_county), max(n_prop))
for(i in 1:nrow(n_prop_county)){
  vec <- n_prop_county[i, which(!is.na(n_prop_county[i,]))]
  n_prop_county1[i, 1:n_prop[i]] <- vec
}

pp <- all_county_units |>
  group_by(Property) |>
  mutate(timestep_idx = 1:n()) |>
  ungroup() |>
  select(Property, timestep, timestep_idx) |>
  pivot_wider(names_from = timestep_idx,
              values_from = timestep,
              id_cols = Property)

n_property <- max(take_df$property_idx)
all_pp <- tibble()
for(i in 1:n_property){
  sub <- filter(take_df, property_idx == i)
  reps <- tibble(
    property_idx = i,
    PPNum = min(sub$PPNum):max(sub$PPNum),
    timestep = PPNum
  ) |>
    mutate(timestep = timestep - min(timestep) + 1)
  all_pp <- bind_rows(all_pp, reps)
}

all_pp_wide_prop <- all_pp |>
  pivot_wider(names_from = timestep,
              values_from = PPNum)

pop_growth_lookup <- all_pp |>
  group_by(property_idx) |>
  filter(PPNum < max(PPNum)) |>
  ungroup() |>
  select(-PPNum) |>
  mutate(H = 1:n()) |>
  pivot_wider(values_from = H,
              names_from = timestep) |>
  select(-property_idx)

all_pp_wide <- all_pp_wide_prop |>
  select(-property_idx)
n_pp_include <- apply(all_pp_wide, 1, function(x) max(which(!is.na(x))))

# create data frame with all timesteps for each property
st_all <- take_df |>
  # map_dfr(~ create_all_timesteps(sample_units, .)) |>
  mutate(state_cd = if_else(st_gsa_state_cd < 10,
                            paste0("0", st_gsa_state_cd),
                            as.character(st_gsa_state_cd)),
         fips = paste0(state_cd, countyfp)) |>
  rowwise() |>
  mutate(midpoint = ifelse(start.date == end.date,
                           as.numeric(start.date),
                           as.numeric(start.date) +
                             (as.numeric(end.date) - as.numeric(start.date)) / 2)) |>
  suppressMessages()

# impose stochastic ordering of events by adding jitter
# we are assuming the following order of events when the events have the same midpoint
# e.g., are on the same day:
# 1. (trap or snare), with order random
# 2. (heli or plane), with order random
# 3. hunting
st_all$jittered_midpoint <- NA
for (i in 1:nrow(st_all)) {
  if(is.na(st_all$midpoint[i])){
    st_all$jittered_midpoint[i] <- 1
  } else if (st_all$Method[i] %in% c('TRAPS', 'SNARE')) {
    st_all$jittered_midpoint[i] <- st_all$midpoint[i] + runif(1, min = 0, max = .01)
  } else if (st_all$Method[i] %in% c('FIXED WING', 'HELICOPTER')) {
    st_all$jittered_midpoint[i] <- st_all$midpoint[i] + runif(1, min = .02, max = .03)
  } else if (st_all$Method[i] == "FIREARMS"){
    st_all$jittered_midpoint[i] <- st_all$midpoint[i] + runif(1, min = .04, max = .05)
  }
}

# now compute orders of events based on jittered midpoints
st_order <- st_all |>
  ungroup() |>
  group_by(st_gsa_state_cd, County, Property, PPNum) |>
  mutate(order = order(jittered_midpoint)) |>
  ungroup() |>
  arrange(Property, PPNum, order) |>
  mutate(p = 1:n())

# Generate start and end indices for previous surveys ---------------------
st_order$start <- 0
st_order$end <- 0

pb <- txtProgressBar(max = nrow(st_order), style = 3)
for (i in 1:nrow(st_order)) {
  if (st_order$order[i] > 1) {
    idx <- which(st_order$County == st_order$County[i] &
                   st_order$Property == st_order$Property[i] &
                   st_order$PPNum == st_order$PPNum[i] &
                   st_order$order < st_order$order[i])
    st_order$start[i] <- idx[1]
    st_order$end[i] <- idx[length(idx)]
    tar_assert_identical(idx, st_order$start[i]:st_order$end[i])
  }
  setTxtProgressBar(pb, i)
}
close(pb)

sampled_units <- st_order |>
  select(County, Property, property_idx, PPNum, timestep) |>
  distinct() |>
  mutate(n_id = 1:n())

county_areas_xl <- readxl::read_xls("data/counties/area.xls")
county_areas <- county_areas_xl |>
  rename(County = STCOU) |>
  select(Areaname, County, LND110210D) |>
  mutate(county_area_km2 = LND110210D * 2.589) |>
  select(County, county_area_km2)

county_sampled_units <- sampled_units |>
  select(County, PPNum) |>
  distinct() |>
  mutate(m_id = 1:n()) |>
  left_join(county_areas)

sampled_units <- left_join(sampled_units, county_sampled_units)

n_timesteps <- sampled_units |>
  group_by(Property) |>
  count() |>
  pull(n)

timestep <- sampled_units |>
  select(-n_id, -m_id, -county_area_km2, -property_idx) |>
  pivot_wider(names_from = timestep,
              values_from = PPNum) |>
  select(-County, -Property)

y_rem <- st_order |>
  group_by(property_idx, PPNum) |>
  summarise(ysum = sum(take)) |>
  ungroup()

N_init <- st_order |>
  group_by(property_idx, timestep) |>
  summarise(N = round((sum(take) + 1)/0.15)) |>
  ungroup() |>
  pivot_wider(names_from = timestep,
              values_from = N) |>
  select(-property_idx) |>
  as.matrix()

y_sum_wide <- left_join(all_pp, y_rem) |>
  mutate(ysum = if_else(is.na(ysum), 0, ysum)) |>
  select(-PPNum) |>
  pivot_wider(names_from = timestep,
              values_from = ysum) |>
  select(-property_idx)

sum_area_surveyed <- st_order |>
  group_by(County, PPNum) |>
  distinct() |>
  summarise(sum_area_surveyed = sum(property_area_km2))

# for spatial process, will add later
# fips <- unique(st_order$fips)
# states <- unique(st_order$state)
#
# counties <- st_order |>
#   select(state, cnty_name) |>
#   distinct() |>
#   pull(cnty_name)
#
#
# # need to figure out which counties have been removed
# shp <- readOGR(file.path("data", "counties"), "dtl.cnty.lower48.meters")
# counties_missing <- setdiff(as.character(fips), shp@data$FIPS)
#
# shp <- shp[shp@data$FIPS %in% as.character(fips),]
#
# st_order |> filter(fips %in% counties_missing)
#
# sf_use_s2(FALSE)
# nb <- poly2nb(shp, row.names = shp$FIPS)
# D <- nb2WB(nb)
# tar_assert_equal_lengths(D$num, fips)
#
# # remove islands from neighbors list
# islands <- which(D$num == 0)
# non_islands <- which(D$num != 0)
# fips_non_islands <- attributes(nb)$region.id[non_islands]
# shp_non_islands <- shp[shp@data$FIPS %in% fips_non_islands,]
# nb_non_islands <- poly2nb(shp_non_islands, row.names = shp_non_islands$FIPS)
# D_non_islands <- nb2WB(nb_non_islands)

# mean litter size year from VerCauteren et al. 2019 pg 63
data_litter_size <- round(c(5.6, 6.1, 5.6, 6.1, 4.2, 5.0, 5.0, 6.5, 5.5, 6.8,
                      5.6, 5.9, 4.9, 5.1, 4.5, 4.7, 5.3, 5.7, 7.4, 8.4,
                      4.7, 4.9, 3.0, 3.0, 4.8, 4.8, 4.2, 5.4, 4.7, 5.2, 5.4))

Xp <- tar_read("group_d_train")$X_p

survey_obs <- tar_read("survey_obs")

# generate centered and scaled versions of these numeric variables
center_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

group_d <- survey_obs |>
  rename(County = FIPS) |>
  select(County, rural.road.density, mean.ruggedness, mean.canopy.density) |>
  distinct()

st_x <- left_join(st_order, group_d) |>
  mutate(method_factor = as.factor(Method),
         c_road_den = center_scale(rural.road.density),
         c_rugged = center_scale(mean.ruggedness),
         c_canopy = center_scale(mean.canopy.density))
nrow(st_x)

X_p <- st_x |>
  select(c_road_den, c_rugged, c_canopy) |>
  as.matrix()

tar_assert_identical(nrow(X_p), nrow(st_order))

areas <- all_county_units |>
  select(Property, property_area_km2) |>
  distinct() |>
  pull(property_area_km2)

st_order <- st_order |>
  group_by(Property, PPNum) |>
  mutate(n_id = cur_group_id()) |>
  ungroup()

constants <- list(
  n_timesteps = n_timesteps,
  n_county = length(unique(st_order$County)),
  n_survey = nrow(st_order),
  n_ls = length(data_litter_size),
  n_property = n_property,
  n_county_units = nrow(county_sampled_units),
  n_prop = n_prop,
  n_first_survey = length(which(st_order$order == 1)),
  n_not_first_survey = length(which(st_order$order != 1)),
  n_units = nrow(sampled_units),
  n_id = st_order$n_id,
  n_method = max(st_order$method_idx),
  n_pp = max(n_pp_include),
  n_pp_prop = n_pp_include,
  all_pp = as.matrix(all_pp_wide),
  m_p = ncol(X_p),
  first_survey = which(st_order$order == 1),
  not_first_survey = which(st_order$order != 1),
  p_property_idx = st_order$property_idx,
  p_pp_idx = st_order$PPNum,
  start = st_order$start,
  end = st_order$end,
  PPNum = as.matrix(timestep),
  method = st_order$method_idx,
  property_x = sampled_units$property_idx,
  pp_x = sampled_units$PPNum,
  pp_len = 28,
  # phi_prior_mean = phi_prior$surv_mu,
  # phi_prior_tau = sigma_2_tau(phi_prior$surv_sd),
  M_lookup = n_prop_county1,
  county = as.numeric(as.factor(st_order$County)),
  pH = as.matrix(pop_growth_lookup),
  property_area = areas,
  log_pi = log(pi)
)

source("R/priors.R")
phi_prior <- create_surv_prior()
constants$phi_mu_a <- 1
constants$phi_mu_b <- 1
# constants$phi_mu_a <- phi_prior$alpha
# constants$phi_mu_b <- phi_prior$beta

data <- list(
  y = st_order$take,
  J = data_litter_size,
  rem = as.matrix(y_sum_wide),
  sum_prop_area = sum_area_surveyed$sum_area_surveyed,
  county_area = county_sampled_units$county_area_km2,
  X_p = X_p,
  log_effort_per = log(st_order$effort_per),
  effort_per = st_order$effort_per,
  n_trap_m1 = st_order$trap_count - 1,
  log_survey_area_km2 = log(st_order$property_area_km2)
)

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

property_lookup <- st_order |>
  select(Property, property_idx, property_area_km2, state, cnty_name, fips) |>
  distinct()

method_lookup <- st_order |>
  select(method_idx, Method) |>
  distinct() |>
  arrange(method_idx)

write_rds(
  list(
    data = data,
    constants = constants,
    property_lookup = property_lookup,
    method_lookup = method_lookup,
    unit_lookup = sampled_units,
    params_check = params_check,
    all_pp = all_pp
  ),
  file = file.path(dest, "nimbleList.rds")
)


source("R/dm_inits.R")
post_dir <- "out/simulation/modifiedDM_betaSurvival_uninformative/poisson/customMCMC/simulation_1"

# custom_samplers <- tribble(
#   ~node,            ~type,
#   "tau_phi",        "slice",
#   "tau_dem",        "slice",
#   "log_mean_ls",    "slice",
#   "logit_mean_phi", "slice",
#   "log_rho",        "AF_slice"
# )

custom_samplers <- NULL
monitors_add <- c("xn", "p", "log_theta")
model_flags <- list(likelihood = likelihood,
                    spatial = FALSE,
                    demographic_stochasticity = TRUE)
n_iter <- 25000
max_iter <- 1000000

if(run_parallel){
  source("R/nimble_dm_2.R")
  source("R/run_nimble_parallel.R")

  message("==== Build cluster ===")
  cl <- makeCluster(3)
  system.time(
    samples <- run_nimble_parallel(
      cl = cl,
      model_code = modelCode,
      model_data = data,
      model_constants = constants,
      model_inits = inits,
      inits_dir = post_dir,
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
  )
  stopCluster(cl)
  message("Stop cluster")

  # cl <- makeCluster(3)
  # model_code = modelCode
  # model_data = data
  # model_constants = constants
  # model_inits = inits
  # init = model_inits()
  # n_iter = 250
  # n_burnin = 0
  # params_check = params_check
  # custom_samplers = custom_samplers
  # model_flags = model_flags
  # state.col = "xn"
  # max_iter = 5000
  # resetMV = TRUE
  # save_iter = TRUE
  # effective_size = 5000
  # max_psrf = 500
  # calculate = FALSE
  # use_conjugacy = TRUE
  # dest = dest
  # stopCluster(cl)

} else {

  source("R/dm_inits.R")
  i_test <- inits(data, constants)
  # source("R/nimble_ZIP.R")
  source("R/nimble_dm_2.R")
  spatial <- FALSE
  Rmodel <- nimbleModel(
    code = modelCode,
    constants = constants,
    data = data,
    inits = i_test,
    calculate = TRUE
  )

  # sim_nodes <- Rmodel$getDependencies(c("log_mean_ls", "logit_mean_phi"),
  #                                     self = FALSE,
  #                                     downstream = TRUE)
  # Rmodel$simulate(sim_nodes)
  # warnings()

  # Rmodel$calculate()

  # check initialization
  Rmodel$initializeInfo()

  # default MCMC configuration
  mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)
  mcmcConf$addMonitors(monitors_add)

  if(!is.null(custom_samplers)){
    for(i in seq_len(nrow(custom_samplers))){
      node <- custom_samplers$node[i][[1]]
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
  samples <- runMCMC(Cmcmc, nchains = 1, niter = 3)

  params_check <- c(
    "beta_p",
    "beta1",
    "log_gamma",
    "log_rho",
    "logit_mean_phi",
    "tau_phi",
    "log_mean_ls",
    "p_mu"
  )
  all_nodes <- colnames(samples[[1]])
  j <- unlist(lapply(params_check, function(x) grep(x, all_nodes)))

  params <- mcmc.list()
  for(i in 1:length(samples)){
    params[[i]] <- as.mcmc(samples[[i]][,j])
  }
  plot(params)

}





for(i in 1:constants$n_property){
  for(j in 2:constants$n_pp_prop[i]){
    xx <- Rmodel$calculate(paste0("R[", i, ", ", j-1, "]"))
    if(is.nan(xx)) stop()
  }
}







