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

config_name <- "dm"
config_ls <- get(config = config_name)
list2env(config_ls, .GlobalEnv)

data <- read_csv("data/insitu/sample_dynamic.csv")

property_filter <- data |>
  group_by(agrp_prp_id) |>
  summarise(n = sum(take)) |>
  filter(n > 0) |>
  pull(agrp_prp_id)

sample_units <- data |>
  filter(agrp_prp_id %in% property_filter) |>
  select(-alws_agrprop_id, -property.size, -two_plus_takes, -n_timesteps, -p, -`...1`) |>
  arrange(agrp_prp_id, timestep_4) |>
  rename(timestep = timestep_4) |>
  mutate(method = if_else(grepl("TRAPS", method), "TRAPS", method))

properties <- sample_units |>
  select(agrp_prp_id, timestep) |>
  distinct() |>
  group_by(agrp_prp_id) |>
  count() |>
  filter(n >= 3) |>
  pull(agrp_prp_id) |>
  unique() |>
  sort()

sample_units <- sample_units |>
  filter(agrp_prp_id %in% properties)

### need to fill in missing timesteps for each property
# each spacio-temporal unit
space_time_units <- sample_units |>
  select(agrp_prp_id, timestep) |>
  distinct()




property_info <- sample_units |>
  select(agrp_prp_id, state, cnty_name, st_gsa_state_cd, cnty_gsa_cnty_cd, fips,
         property_area_km2, state_abr, countyfp) |>
  distinct()

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

st_order <- st_order |>
  group_by(agrp_prp_id, timestep) |>
  mutate(t = 1:n(),
         ysum = cumsum(take) - take) |>
  ungroup()

ts_idx <- st_order |>
  select(agrp_prp_id, timestep) |>
  distinct() |>
  group_by(agrp_prp_id) |>
  mutate(ts_idx = 1:n())

st_order <- left_join(st_order, ts_idx)

n_county_idx <- st_order |>
  select(agrp_prp_id, fips) |>
  distinct() |>
  pull(fips) |>
  as.factor() |>
  as.numeric()

p_timestep_idx <- st_order$ts_idx
y_sum <- st_order$ysum

county_areas_xl <- readxl::read_xls("data/counties/area.xls")
county_areas <- county_areas_xl |>
  select(Areaname, STCOU, LND110210D) |>
  mutate(area_km2 = LND110210D * 2.589) |>
  rename(fips = STCOU)

# n_county_units
# sum_prop_area
# county_area
# n_prop_cnty
# n_prop
#
# all_county_units <- st_order |>
#   select(fips, timestep) |>
#   left_join(county_areas) |>
#   select(-LND110210D) |>
#   mutate(m_county_idx = as.numeric(as.factor(fips))) |>
#   distinct() |>
#   group_by(fips, timestep) |>
#   mutate(t = 1:n()) |>
#   ungroup()
#
# n_county_units <- nrow(all_county_units)
#
# st_order |>
#   select(fips, timestep, agrp_prp_id, t)
#
# for(i in 1:n_county_units){
#   sub <- st_order |>
#     filter(fips == all_county_units$fips[i],
#            timestep == all_county_units)
# }


max_timesteps <- max(n_timesteps_prop$n)
max_reps <- max(reps_df$n)

y_wide <- H <- array(NA, dim = c(n_property, max_timesteps, max_reps))
n_reps <- timestep <-  matrix(NA, n_property, max_timesteps)
n_timesteps <- rep(NA, n_property)
all_pp <- tibble()
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
  y_wide[i, 1:n_timesteps[i], 1:ncol(y_sub)] <- y_sub

  # p_sub <- prop_sub |>
  #   select(-take) |>
  #   pivot_wider(names_from = order, values_from = p) |>
  #   select(-timestep) |>
  #   as.matrix()
  # H[i, 1:n_timesteps[i], 1:ncol(p_sub)] <- p_sub

  t_sub <- timesteps_df |>
    filter(agrp_prp_id == prop) |>
    pull(timestep)
  timestep[i, 1:length(t_sub)] <- t_sub

  reps <- tibble(
    agrp_prp_id = prop,
    timestep = min(t_sub):max(t_sub),
    pp = timestep,
  ) |>
    mutate(timestep = timestep - min(timestep) + 1)
  all_pp <- bind_rows(all_pp, reps)

  # reps <- reps_df |>
  #   filter(agrp_prp_id == prop) |>
  #   pull(n)
  # n_reps[i, 1:n_timesteps[i]] <- reps
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

all_pp_wide <- all_pp |>
  pivot_wider(names_from = timestep,
              values_from = pp) |>
  select(-agrp_prp_id)

n_pp <- apply(all_pp_wide, 1, function(x) max(which(!is.na(x))))

model_flags <- list(
  process_type = process_type,
  spatial = spatial,
  property_obs_effect = property_obs_effect,
  likelihood = likelihood
)

data <- list(
  y_wide = y_wide,
  y_sum = st_order$ysum,
  y = st_order$take,
  X_p = X_p,
  effort_per = st_order$effort_per,
  log_effort_per = log(st_order$effort_per),
  n_trap_m1 = st_order$trap_count - 1,
  log_survey_area_km2 = log(st_order$property_area_km2)
)

source("R/priors.R")
phi_prior <- create_surv_prior(logit = TRUE, sd_inflate = 1.2) |>
  suppressMessages()

constants <- list(
  # n_reps = n_reps,
  n_timesteps = n_timesteps,
  # n_pp = max(timestep, na.rm = TRUE),
  n_edges = length(D$adj),
  n_county = length(counties),
  n_survey = nrow(st_order),
  n_property = nrow(properties_cnty),
  n_beta_p = ncol(X_p) + 1,
  n_first_survey = length(which(st_order$order == 1)),
  n_not_first_survey = length(which(st_order$order != 1)),
  n_trap_snare = length(which(method.vec %in% 4:5)),
  n_shooting = length(which(method.vec %in% 1:3)),
  n_units = nrow(unit_lookup),
  n_method = length(levels(st_order$method_factor)),
  n_pp = max(n_pp),
  n_pp_prop = n_pp,
  all_pp = as.matrix(all_pp_wide),
  # H = H,
  # reps = st_order$order,
  p_county_idx = st_order$fips |> as.factor() |> as.numeric(),
  n_county_idx = n_county_idx,
  log_pi = log(pi),
  adj = D$adj,
  weights = D$weights,
  num = D$num,
  m_p = ncol(X_p),
  first_survey = which(st_order$order == 1),
  not_first_survey = which(st_order$order != 1),
  p_property_idx = st_order$agrp_prp_id |> as.factor() |> as.numeric(),
  p_timestep_idx = p_timestep_idx,
  start = st_order$start,
  end = st_order$end,
  timestep = timestep,
  method = method.vec,
  ts_idx = which(method.vec %in% 4:5),
  shooting_idx = which(method.vec %in% 1:3),
  property_x = unit_lookup$property_idx,
  timestep_x = unit_lookup$timestep,
  phi_prior_mean = phi_prior$surv_mu,
  phi_prior_sd = phi_prior$surv_sd
)

dest <- file.path(out_dir, model_dir, process_type)
if(!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

write_rds(
  list(
    data = data,
    constants = constants,
    property_lookup = property_lookup,
    unit_lookup = unit_lookup,
    model_flags = model_flags,
    params_check = params_check,
    all_pp = all_pp
  ),
  file = file.path(dest, "nimbleList.rds")
)


source(model_file)
source("R/functions_nimble.R")

# inits <- make_inits_function(NULL, constants = constants, data = data)
# inits <- make_inits_function(file.path(out_dir, inits_dir))

zeta_pp <- if_else(grepl("zeta_pp", process_type), TRUE, FALSE)
zeta_constant <- !zeta_pp

phi_pp <- if_else(grepl("phi_pp", process_type), TRUE, FALSE)
phi_constant <- !phi_pp

inits <- make_inits_function_dm(NULL, process_type, constants, data,
                                phi_constant, phi_pp, zeta_constant, zeta_pp)

# custom_samplers <- tibble(
#   node = c("p_mu"),
#   type = c("AF_slice")
# )

custom_samplers <- NULL


# for(i in seq_len(constants$n_method)){
#   beta_node <- paste0("beta_p[", i, ", ", 1:constants$n_beta_p, "]")
#   beta_type <- "RW_block"
#
#   # rho_node <- paste0("log_rho[", i, "]")
#   # rho_type <- "slice"
#   #
#   # gamma_node <- paste0("log_gamma[", i, "]")
#   # gamma_type <- "slice"
#
#   p_node <- paste0("p_mu[", i, "]")
#   p_type <- "slice"
#
#   cs_i <- tibble(
#     node = c(list(beta_node), p_node),
#     type = c(beta_type, p_type)
#   )
#
#   custom_samplers <- bind_rows(custom_samplers, cs_i)
#
# }


if(run_parallel){
  source(model_file)
  # source("R/functions_nimble.R")

  message("==== Build cluster ===")
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
      custom_samplers = NULL,
      model_flags = model_flags,
      calculate = FALSE,
      state.col = "xn",
      max_iter = 500000,
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
  # calculate = FALSE
  # use_conjugacy = TRUE
  # dest = dest
  # stopCluster(cl)

} else {

  likelihood_nb <- if_else(likelihood == "nb", TRUE, FALSE)
  likelihood_binom <- ifelse(likelihood_nb, FALSE, TRUE)

  inits <- make_inits_function_dm(NULL, process_type, constants, data,
                                  phi_constant, phi_pp, zeta_constant, zeta_pp)
  i_test <- inits()

  source(model_file)
  Rmodel <- nimbleModel(
    code = modelCode,
    constants = constants,
    data = data,
    inits = i_test,
    calculate = FALSE
  )
  Cmodel <- compileNimble(Rmodel)
  Cmodel$calculate()

  # check initialization
  Cmodel$initializeInfo()

  # default MCMC configuration
  mcmcConf <- configureMCMC(Cmodel, useConjugacy = TRUE)#, monitors = c("xn", params_check))

  if(!is.null(custom_samplers)){
    for(i in seq_len(nrow(custom_samplers))){
      node <- custom_samplers$node[i][[1]]
      type <- custom_samplers$type[i]
      mcmcConf$removeSampler(node)
      mcmcConf$addSampler(node, type)
    }
  }

  # these print statements will not display when running in parallel
  # mcmcConf$printMonitors()
  # mcmcConf$printSamplers(byType = TRUE)
  Rmcmc <- buildMCMC(mcmcConf)
  Cmcmc <- compileNimble(Rmcmc)
  samples <- runMCMC(Cmcmc, nchains = 1, niter = 5)

  ss  <- as.matrix(samples)
  # plot(exp(ss[, "log_r_mu"]), type="l")
  # plot(exp(ss[, "xn[211]"]), type="l")
  # plot(exp(ss[, "log_rho[4]"]), type="l")
  # plot(exp(ss[, "log_zeta_global"]), type="l")
  # plot(ilogit(ss[, "logit_phi"]), type="l")
  # plot(exp(ss[, "beta_r"]), type="l")
}





