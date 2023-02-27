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

config_name <- "ricker"

config_ls <- get(config = config_name)
list2env(config_ls, .GlobalEnv)


data <- read_csv("data/insitu/sample_dynamic.csv")
sample_units <- data |>
  select(-alws_agrprop_id, -property.size, -two_plus_takes, -n_timesteps, -p, -`...1`) |>
  arrange(agrp_prp_id, timestep_4) |>
  rename(timestep = timestep_4) |>
  mutate(method = if_else(grepl("TRAPS", method), "TRAPS", method))

### need to fill in missing timesteps for each property
# each spacio-temporal unit
space_time_units <- sample_units |>
  select(agrp_prp_id, timestep) |>
  distinct()

# all properties
properties <- space_time_units |>
  pull(agrp_prp_id) |>
  unique() |>
  sort()

property_info <- sample_units |>
  select(agrp_prp_id, state, cnty_name, st_gsa_state_cd, cnty_gsa_cnty_cd, fips,
         property_area_km2, state_abr, countyfp) |>
  distinct()

# function to fill in missing timesteps and information for a given property
# create_all_timesteps <- function(df, prop){
#   # take data for property
#   take <- df |>
#     filter(agrp_prp_id == prop) |>
#     arrange(timestep)
#
#   timesteps <- unique(take$timestep)
#   all_time <- min(timesteps):max(timesteps)
#   missing_timesteps <- setdiff(all_time, timesteps)
#
#   if(length(missing_timesteps) > 0){
#     # property info
#     info <- property_info |>
#       filter(agrp_prp_id == prop)
#
#     st <- tibble(agrp_prp_id = prop,
#                  timestep = rep(missing_timesteps, each=2)) |>
#       left_join(info) |>
#       bind_rows(take) |>
#       arrange(timestep)
#     return(st)
#   } else {
#     return(take)
#   }
# }

# create data frame with all timesteps for each property
st_all <- sample_units |>
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
  } else if (st_all$method[i] %in% c('TRAPS, CAGE', 'SNARE')) {
    st_all$jittered_midpoint[i] <- st_all$midpoint[i] + runif(1, min = 0, max = .01)
  } else if (st_all$method[i] %in% c('FIXED WING', 'HELICOPTER')) {
    st_all$jittered_midpoint[i] <- st_all$midpoint[i] + runif(1, min = .02, max = .03)
  } else if (st_all$method[i] == "FIREARMS"){
    st_all$jittered_midpoint[i] <- st_all$midpoint[i] + runif(1, min = .04, max = .05)
  }
}

# now compute orders of events based on jittered midpoints
st_order <- st_all |>
  ungroup() |>
  group_by(st_gsa_state_cd, cnty_name, countyfp, agrp_prp_id, timestep) |>
  mutate(order = order(jittered_midpoint)) |>
  ungroup() |>
  arrange(agrp_prp_id, timestep, order) |>
  mutate(p = 1:n())

# Generate start and end indices for previous surveys ---------------------
st_order$start <- 0
st_order$end <- 0

pb <- txtProgressBar(max = nrow(st_order), style = 3)
for (i in 1:nrow(st_order)) {
  if (st_order$order[i] > 1) {
    idx <- which(st_order$fips == st_order$fips[i] &
                   st_order$agrp_prp_id == st_order$agrp_prp_id[i] &
                   st_order$timestep == st_order$timestep[i] &
                   st_order$order < st_order$order[i])
    st_order$start[i] <- idx[1]
    st_order$end[i] <- idx[length(idx)]
    tar_assert_identical(idx, st_order$start[i]:st_order$end[i])
  }
  setTxtProgressBar(pb, i)
}
close(pb)



fips <- unique(st_order$fips)
states <- unique(st_order$state)

counties <- st_order |>
  select(state, cnty_name) |>
  distinct() |>
  pull(cnty_name)


# st_d <- tar_read("st_d")
# train_d <- tar_read("group_d_train")
# group_d_train <- train_d$group_d
# short_basis <- tar_read("short_basis")

# need to figure out which counties have been removed
shp <- readOGR(file.path("data", "counties"), "dtl.cnty.lower48.meters")
counties_missing <- setdiff(as.character(fips), shp@data$FIPS)

shp <- shp[shp@data$FIPS %in% as.character(fips),]

st_order |> filter(fips %in% counties_missing)

sf_use_s2(FALSE)
nb <- poly2nb(shp, row.names = shp$FIPS)
D <- nb2WB(nb)
tar_assert_equal_lengths(D$num, fips)

# remove islands from neighbors list
islands <- which(D$num == 0)
non_islands <- which(D$num != 0)
fips_non_islands <- attributes(nb)$region.id[non_islands]
shp_non_islands <- shp[shp@data$FIPS %in% fips_non_islands,]
nb_non_islands <- poly2nb(shp_non_islands, row.names = shp_non_islands$FIPS)
D_non_islands <- nb2WB(nb_non_islands)


properties_cnty <- st_order |>
  select(cnty_name, agrp_prp_id) |>
  arrange(agrp_prp_id) |>
  distinct()

n_property <- nrow(properties_cnty)

# the number timesteps in each property
timesteps_df <- st_order |>
  select(agrp_prp_id, timestep) |>
  distinct() |>
  group_by(agrp_prp_id)

n_timesteps_prop <- timesteps_df |>
  tally() |>
  ungroup() |>
  left_join(properties_cnty) |>
  arrange(agrp_prp_id)

# the number of repeated removals in each space time unit
reps_df <- st_order |>
  group_by(agrp_prp_id, timestep) |>
  tally() |>
  ungroup() |>
  left_join(properties_cnty) |>
  arrange(agrp_prp_id)

max_timesteps <- max(n_timesteps_prop$n)
max_reps <- max(reps_df$n)

y <- H <- array(NA, dim = c(n_property, max_timesteps, max_reps))
n_reps <- timestep <-  matrix(NA, n_property, max_timesteps)
n_timesteps <- rep(NA, nrow(timesteps_df))
for(i in 1:n_property){
  prop <- properties[i]
  n_timesteps[i] <- n_timesteps_prop |>
    filter(agrp_prp_id == prop) |>
    pull(n)

  prop_sub <- st_order |>
    filter(agrp_prp_id == prop) |>
    select(take, timestep, order, p)

  y_sub <- prop_sub |>
    select(-p) |>
    pivot_wider(names_from = order, values_from = take) |>
    select(-timestep) |>
    as.matrix()
  y[i, 1:n_timesteps[i], 1:ncol(y_sub)] <- y_sub

  p_sub <- prop_sub |>
    select(-take) |>
    pivot_wider(names_from = order, values_from = p) |>
    select(-timestep) |>
    as.matrix()
  H[i, 1:n_timesteps[i], 1:ncol(p_sub)] <- p_sub

  t_sub <- timesteps_df |>
    filter(agrp_prp_id == prop) |>
    pull(timestep)
  timestep[i, 1:length(t_sub)] <- t_sub

  reps <- reps_df |>
    filter(agrp_prp_id == prop) |>
    pull(n)
  n_reps[i, 1:n_timesteps[i]] <- reps
}

for(i in 1:n_property){
  for(t in 2:n_timesteps[i]){
    cc <- (timestep[i, t-1]+1):timestep[i, t]
    if(any(is.na(cc))) print(i)
  }
}

# y[i, 1:n_timesteps[i],]
# n_reps[i,1:n_timesteps[i]]
# delta_t[i, 1:n_timesteps[i]]

st_order$method_factor <- as.factor(st_order$method)

method.vec <- as.numeric(st_order$method_factor)
# method.init <- method.vec
# method.init[is.na(method.vec)] <- 5
# method.init[!is.na(method.vec)] <- NA

Xp <- tar_read("group_d_train")$X_p

survey_obs <- tar_read("survey_obs")
group_d <- survey_obs |>
  # filter(FIPS %in% as.double(fips)) |>
  select(#FIPS,
         state, countyname,
         starts_with("c_"), rural.road.density, prop.pub.land,
         mean.ruggedness, mean.canopy.density)  |>
  # mutate(fips = as.character(FIPS)) |>
  distinct() |>
  rename(cnty_name = countyname,
         state_abr = state) |>
  mutate(cnty_name = if_else(grepl("ST ", cnty_name), gsub("ST ", "ST. ", cnty_name), cnty_name))

st_x <- left_join(st_order, group_d) |>
  mutate(method_factor = as.factor(method))
nrow(st_x)

X_p <- st_x |>
  select(c_road_den, c_rugged, c_canopy)

tar_assert_identical(nrow(X_p), nrow(st_order))

property_lookup <- st_order |>
  mutate(cnty_idx = as.numeric(as.factor(fips)),
         property_idx = as.numeric(as.factor(agrp_prp_id)),
         island = if_else(fips %in% fips_non_islands, 0, 1)) |>
  select(cnty_name, state, cnty_idx, property_idx, island, fips, property_area_km2) |>
  distinct()

unit_lookup <- st_order |>
  rename(pp = timestep) |>
  select(agrp_prp_id, pp) |>
  distinct() |>
  group_by(agrp_prp_id) |>
  mutate(timestep = 1:n()) |>
  ungroup() |>
  mutate(property_idx = as.numeric(as.factor(agrp_prp_id))) |>
  select(agrp_prp_id, timestep, pp, property_idx)

model_flags <- list(
  process_type = process_type,
  spatial = spatial,
  property_obs_effect = property_obs_effect
)

data <- list(
  y = y,
  X_p = X_p,
  effort_per = st_order$effort_per,
  log_effort_per = log(st_order$effort_per),
  n_trap_m1 = st_order$trap_count - 1,
  log_survey_area_km2 = log(st_order$property_area_km2)
)

constants <- list(
  H = H,
  n_reps = n_reps,
  n_pp = max(timestep, na.rm = TRUE),
  n_timesteps = n_timesteps,
  n_reps = n_reps,
  county_idx = property_lookup$cnty_idx,
  log_pi = log(pi),
  adj = D$adj,
  weights = D$weights,
  num = D$num,
  n_edges = length(D$adj),
  n_county = length(counties),
  n_survey = nrow(st_order),
  n_property = nrow(properties_cnty),
  m_p = ncol(X_p),
  n_beta_p = ncol(X_p) + 1,
  first_survey = which(st_order$order == 1),
  n_first_survey = length(which(st_order$order == 1)),
  not_first_survey = which(st_order$order != 1),
  n_not_first_survey = length(which(st_order$order != 1)),
  n_method = length(levels(st_order$method_factor)),
  p_property_idx = st_order$agrp_prp_id |> as.factor() |> as.numeric(),
  start = st_order$start,
  end = st_order$end,
  timestep = timestep,
  method = method.vec,
  ts_idx = which(method.vec %in% 4:5),
  n_trap_snare = length(which(method.vec %in% 4:5)),
  shooting_idx = which(method.vec %in% 1:3),
  n_shooting = length(which(method.vec %in% 1:3)),
  n_units = nrow(unit_lookup),
  property_x = unit_lookup$property_idx,
  timestep_x = unit_lookup$timestep,
  log_k = log(property_lookup$property_area_km2) + log(25)
  # log_property_area = log(property_lookup$property_area_km2),
  # log_max_density = log(15)
)

inits <- function(){
  N_init <- (rowSums(data$y, dims = 2, na.rm = TRUE) + 1)
  ls <- list(
    N = N_init,
    x = log(N_init),
    log_gamma = rnorm(5),
    log_rho = rnorm(5),
    beta_p = matrix(rnorm(constants$n_method*constants$n_beta_p, 0, 0.1), constants$n_method, constants$n_beta_p),
    p_mu = rnorm(constants$n_method),
    log_r_mu = runif(1, 0, 1),
    tau_proc = runif(1, 1, 10)
  )
  return(ls)
}



dest <- file.path(out_dir, model_dir)
if(!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

write_rds(
  list(
    data = data,
    constants = constants,
    property_lookup = property_lookup,
    unit_lookup = unit_lookup,
    model_flags = model_flags,
    params_check = params_check
  ),
  file = file.path(dest, "nimbleList.rds")
)



source(model_file)
source("R/functions_nimble.R")

custom_samplers <- data.frame(
  node = c("beta_p",   "log_rho",  "log_gamma", "p_mu",     "log_r_mu"),
  type = c("AF_slice", "AF_slice", "AF_slice",  "AF_slice", "ess")
)


if(run_parallel){
  message("Build cluster")
  cl <- makeCluster(3)
  system.time(
    samples <- run_nimble_parallel(
      cl = cl,
      model_code = modelCode,
      model_data = data,
      model_constants = constants,
      model_inits = inits,
      n_iter = n_iter,
      n_burnin = 0,
      params_check = params_check,
      custom_samplers = custom_samplers,
      model_flags = model_flags,
      state.col = "xn",
      max_iter = 200000,
      dest = dest,
      resetMV = TRUE,
      save_iter = TRUE
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
  # calculate = TRUE
  # use_conjugacy = FALSE
  # dest = dest
  # stopCluster(cl)

} else {

  exponential <- if_else(process_type == "exponential", TRUE, FALSE)
  ricker <- if_else(process_type == "ricker", TRUE, FALSE)
  gompertz <- if_else(process_type == "gompertz", TRUE, FALSE)

  Rmodel <- nimbleModel(
    code = modelCode,
    constants = constants,
    inits = inits(),
    data = data
  )

  # check initialization
  Rmodel$initializeInfo()

  # default MCMC configuration
  mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE, monitors = c("xn", params_check))

  if(!is.null(custom_samplers)){
    for(i in seq_len(nrow(custom_samplers))){
      node <- custom_samplers$node[i]
      type <- custom_samplers$type[i]
      mcmcConf$removeSampler(node)
      mcmcConf$addSampler(node, type)
    }
  }

  # these print statements will not display when running in parallel
  mcmcConf$printMonitors()
  mcmcConf$printSamplers(byType = TRUE)

  Rmcmc <- buildMCMC(mcmcConf)
  Cmodel <- compileNimble(Rmodel)
  Cmcmc <- compileNimble(Rmcmc)
  samples <- runMCMC(Cmcmc, nchains = 1, niter = 200)

  ss  <- as.matrix(samples)
  # plot(exp(ss[, "log_r_mu"]), type="l")
  # plot(exp(ss[, "xn[211]"]), type="l")
}





