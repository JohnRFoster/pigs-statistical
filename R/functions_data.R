## functions called throughout the target pipeline


## ---------------------------- Data ingest -----------------------------------

get_timestep_df <- function(file) readRDS(file)

get_configs <- function(file) config::get()

get_states <- function(states_filter){
  stcsv <- read_csv("data/counties/statePostalCodes.csv")
  df <- stcsv
  if(states_filter != "all"){
    df <- stcsv |>
      filter(Postal %in% states_filter)
  }
  return(df)
}

get_shp <- function(){
  shp <- readOGR(file.path("data", "counties"), "dtl.cnty.lower48.meters")
  # shp <- shp[shp@data$STATE_NAME == "Missouri",]
  return(shp)
}

get_traps <- function(file){
  read_csv(file) %>%
    mutate(method = "trap",
           effort = trap.nights,
           y = Take,
           #ST_GSA_STATE_CD = ST_FIPS,
           ST_GSA_STATE_CD = ST_GSA_STATE_CD,
           #CNTY_GSA_CNTY_CD = CNTY_FIPS,
           CNTY_GSA_CNTY_CD = CNTY_GSA_CNTY_CD,
           trap_count = trap.count) %>%
    dplyr::select(CNTY_NAME, ST_GSA_STATE_CD, CNTY_GSA_CNTY_CD, start.date,
                  end.date, AGRP_PRP_ID, method, effort, y, TOTAL.LAND,
                  trap_count) %>%
    rename_with(tolower)
}

get_aerial <- function(file){
  read_csv(file) %>%
    mutate(method = CMP_NAME,
           effort = Flight.Hours,
           y = Take,
           #ST_GSA_STATE_CD = ST_FIPS,
           #CNTY_GSA_CNTY_CD = CNTY_FIPS,
           trap_count = VEHICLES) %>%
    dplyr::select(CNTY_NAME, ST_GSA_STATE_CD, CNTY_GSA_CNTY_CD,
                  Start.Date, End.Date, trap_count,
                  AGRP_PRP_ID, method, effort, y, TOTAL.LAND) %>%
    rename_with(tolower) #%>%
  #separate(start.date, into = c('m', 'd', 'yr'), sep = '/') %>%
  # separate(start.date, into = c('yr', 'm', 'd'), sep = '-') %>%
  # mutate(m = sprintf('%02d', as.numeric(m)),
  #        d = sprintf('%02d', as.numeric(d)),
  #        start.date = paste(yr, m, d, sep = '-'),
  #        start.date = as.Date(start.date)) %>%
  # dplyr::select(-m, -d, -yr) %>%
  # separate(end.date, into = c('m', 'd', 'yr'), sep = '/') %>%
  # mutate(m = sprintf('%02d', as.numeric(m)),
  #        d = sprintf('%02d', as.numeric(d)),
  #        end.date = paste(yr, m, d, sep = '-'),
  #        end.date = as.Date(end.date)) %>%
  # dplyr::select(-m, -d, -yr)
}

get_snares <- function(file){
  read_csv(file) %>%
    filter(!is.na(ST_NAME),
           !is.na(CNTY_NAME)) %>%
    mutate(method = 'snare',
           effort = trap.nights,
           y = Take,
           # ST_GSA_STATE_CD = ST_FIPS,
           # CNTY_GSA_CNTY_CD = CNTY_FIPS,
           trap_count = trap.nights / event.length) %>%
    dplyr::select(CNTY_NAME, ST_GSA_STATE_CD, CNTY_GSA_CNTY_CD, start.date,
                  end.date, AGRP_PRP_ID, method, effort, y, TOTAL.LAND,
                  trap_count) %>%
    rename_with(tolower)
}

get_firearms <- function(file){
  read_csv(file) %>%
    filter(!is.na(ST_NAME),
           !is.na(CNTY_NAME)) %>%
    mutate(method = 'firearms',
           effort = Hunt.Hours,
           y = Take,
           # ST_GSA_STATE_CD = ST_FIPS,
           # CNTY_GSA_CNTY_CD = CNTY_FIPS,
           Start.date = WT_WORK_DATE,
           End.date = WT_WORK_DATE,
           trap_count = FIREARMS) %>%
    dplyr::select(CNTY_NAME, ST_GSA_STATE_CD, CNTY_GSA_CNTY_CD, Start.date,
                  End.date, AGRP_PRP_ID, method, effort, y, TOTAL.LAND,
                  trap_count) %>%
    rename_with(tolower)
}



  # Merge all of the take data into one data frame
insitu_all <- function(traps.data, firearms.data, snares.data, aerial.data){
  ins_df <- full_join(traps.data, firearms.data) %>%
    full_join(snares.data) %>%
    full_join(aerial.data) %>%
    distinct %>%
    mutate(insitu_id = 1:n(),
           method = tolower(method)) %>%
    filter(!is.na(y))
  return(ins_df)
}

  # insitu_df <- insitu_all(traps, firearms, snares, aerial)


# resolve duplicate property values - when there are multiple values, take max
# why are we taking the max property value when there are duplicate records?
resolve_duplicate <- function(insitu.data){
  insitu.data |>
    distinct(agrp_prp_id, total.land) %>%
    group_by(agrp_prp_id) %>%
    summarize(n_areas = length(unique(total.land)),
              property_area_acres = max(total.land, na.rm = TRUE),
              # properties with all NA areas get -Inf
              # but the following line changes -Inf to NA
              property_area_acres = ifelse(is.infinite(property_area_acres),
                                           NA,
                                           property_area_acres),
              property_area_km2 = 0.00404686 * property_area_acres) %>%
    ungroup %>%
    arrange(property_area_acres)
  #property_area_df <- property_area_df[property_area_df$property_area_km2 > 1, ]
}

  # property_area_df <- resolve_duplicate(insitu_df)

# merge property areas back into insitu data, filter on area
join_area_filter <- function(ins, prop){
  insitu_property <- ins %>%
    left_join(prop) %>%
    filter(!is.na(property_area_acres),
           property_area_km2 >= 1.8,
           effort > 0)
  return(insitu_property)
}

  # insitu_property <- join_area_filter(insitu_df, property_area_df)

get_fips <- function(file){
  fips <- read_csv(file,
                   col_names = c("state", "statefp", "countyfp",
                                 "countyname", "classfp"),
                   comment = '#')
}


# fetch all fips data from census
merge_fips <- function(ins.p, fips){
  usa_fips <- fips %>%
    mutate(st_gsa_state_cd = parse_number(statefp),
           cnty_gsa_cnty_cd = parse_number(countyfp))

  # merge fips data to insitu pig data
  merged_join <- left_join(ins.p, usa_fips) %>%
    arrange(countyname, start.date, end.date) %>%
    filter(st_gsa_state_cd != 15) # exclude Hawaii

  return(merged_join)
}

# with_fips <- merge_fips(insitu_property, fips)

clean_merge <- function(merged_d){
  min_date <- min(merged_d$start.date)
  max_date <- max(merged_d$end.date)

  interval <- 4 # number of weeks that comprise one 'primary period'

  start_dates <- seq(min_date, max_date, by = paste(interval, "week"))
  end_dates <- c(start_dates[-1] - 1, max_date)
  tar_assert_identical(length(start_dates), length(end_dates))
  tar_assert_true(min(merged_d$start.date) >= min_date)
  tar_assert_true(max(merged_d$start.date) <= max_date)

  timestep_df <- tibble(start_dates, end_dates) %>%
    mutate(timestep = 1:n())
  timestep_df$month <- month(timestep_df$end_dates)
  timestep_df$year <- year(timestep_df$end_dates)

  # for each row in the merged data, insert the integer primary period timestep
  merged_d$timestep <- NA
  pb <- txtProgressBar(max = nrow(merged_d), style = 3)
  for (i in 1:nrow(merged_d)) {
    after_start <- which(timestep_df$start_dates <= merged_d$start.date[i]) %>% max
    before_end <- which(timestep_df$end_dates >= merged_d$end.date[i]) %>% min
    if (after_start == before_end) {
      # then the start and end date is contained within a primary period
      merged_d$timestep[i] <- timestep_df$timestep[before_end]
    } # otherwise, timestep[i] will be left as NA and filtered out later
    setTxtProgressBar(pb, i)
  }

  # restrict take data to efforts that are contained in one primary period
  cleaned_d <- merged_d %>%
    filter(!is.na(timestep)) %>%
    arrange(start.date)

  return(cleaned_d)
}

  # clean_d <- clean_merge(with_fips)

# compute ordering based on time interval midpoints
order_interval <- function(cleaned, prop.area){
  order_df <- cleaned %>%
    left_join(prop.area) %>%
    distinct %>%
    rowwise %>%
    mutate(midpoint = ifelse(start.date == end.date,
                             as.numeric(start.date),
                             as.numeric(start.date) +
                               (as.numeric(end.date) - as.numeric(start.date)) / 2)
    ) %>%
    ungroup

  return(order_df)
}

# order_df <- order_interval(clean_d, property_area_df)

# impose stochastic ordering of events by adding jitter
# we are assuming the following order of events when the events have the same midpoint
# e.g., are on the same day:
# 1. (trap or snare), with order random
# 2. (heli or plane), with order random
# 3. hunting
order_stochastic <- function(order.df){
  order_df <- order.df
  order_df$jittered_midpoint <- NA
  pb <- txtProgressBar(max = nrow(order_df), style = 3)
  for (i in 1:nrow(order_df)) {
    if (order_df$method[i] %in% c('Trap', 'Snare')) {
      order_df$jittered_midpoint[i] <- order_df$midpoint[i] + runif(1, min = 0, max = .01)
    } else if (order_df$method[i] %in% c('Fixed Wing', 'Helicopter')) {
      order_df$jittered_midpoint[i] <- order_df$midpoint[i] + runif(1, min = .02, max = .03)
    } else {
      order_df$jittered_midpoint[i] <- order_df$midpoint[i] + runif(1, min = .04, max = .05)
    }
    setTxtProgressBar(pb, i)
  }
  return(order_df)
}

# order_mid <- order_stochastic(order_df)

# now compute orders of events based on jittered midpoints
order_of_events <- function(order_df, clean_df){
  order_event <- order_df %>%
    ungroup %>%
    group_by(st_gsa_state_cd, cnty_name, countyfp, agrp_prp_id, timestep) %>%
    mutate(order = order(jittered_midpoint),
           has_multi = any(order > 1),
           any_ties = any(duplicated(jittered_midpoint)),
           n_survey = n()) %>%
    arrange(st_gsa_state_cd, cnty_name, countyfp, agrp_prp_id, timestep, order) %>%
    ungroup()

  # merge survey orders back into original data
  clean_order <- clean_df %>%
    left_join(order_event)

  return(order_event)
}

# order_event <- order_of_events(order_mid, clean_d)



# next, link the counties in the shapefile to the counties in the insitu data
get_county_shp <- function(fips, shp){

  # Create a neighborhood adjacency matrix
  IDs <- row.names(as(shp, "data.frame"))

  county_df <- shp@data %>%
    as_tibble %>%
    distinct(STATE_FIPS, CNTY_FIPS) %>%
    rename(statefp = STATE_FIPS, countyfp = CNTY_FIPS) %>%
    right_join(fips) %>%
    filter(state != 'AK') %>%
    mutate(county_fips_factor = factor(countyfp),
           countyname = gsub("County", "", x = countyname),
           countyname = tolower(countyname),
           countyname = trimws(countyname))

  return(county_df)
}


## Load covariates and merge into data --------------------

get_center_nontime_covar <- function(file.spatial, file.crop){
  spatial_covs <- read_csv(file.spatial)

  crop_covs <- read_csv(file.crop) %>%
    dplyr::select(State, FIPS, Group_Name, mean.Prop.Crop) %>%
    spread(Group_Name, mean.Prop.Crop)
  colnames(crop_covs) <- c("State", "FIPS", "Cereals",
                           "Fruit.Nuts", "Other" ,"Root.Tuber", "Vegetables.Melons")

  process_covs <- full_join(spatial_covs, crop_covs) %>%
    mutate(c_hydroden = c(scale(log(.01 + total.hydro.density))),
           c_hetero = c(scale(mean.habitat.hetero)),
           c_lewis = c(scale(mean.lewis.pig.density)),
           c_carnrich = c(scale(mean.carnivore.richness)),
           c_crop = c(scale(mean.crop.cover)),
           c_pasture = c(scale(mean.pasture.cover)),
           c_evergreen = c(scale(evergreen)),
           c_deciduous = c(scale(deciduous)),
           c_cerealCrop = c(scale(Cereals)),
           c_fruitNut = c(scale(Fruit.Nuts)),
           c_rootTuber = c(scale(Root.Tuber)),
           c_vegetablesMelons = c(scale(Vegetables.Melons)),
           c_tree = c(scale(mean.tree.cover)))

  return(process_covs[!is.na(process_covs$c_hydroden),])
}

# timeVar1 <- 'data/covariates/TimeVaryingCompletedPredictors/FINAL.Process.Model.time.varying.crop.acerage.22Dec2017.csv' %>%
#   read_csv %>%
#   spread(Group_Name, Prop.Crop)
# colnames(timeVar1) <- c("State", "FIPS","year", "Crop.Acres", "Total.Acres", "Cereals", "Fruit.Nuts", "Other")
# timeVar1 <- timeVar1[, -c(9,10)]
# # add in a row for each month in each year
# reptimes <- 12
# idx_rep <- rep(1:nrow(timeVar1), reptimes)
# rep_df <- timeVar1[idx_rep, ]
# #timeVar1[1:10,] == rep_df[(nrow(timeVar1)+1):(nrow(timeVar1)+10), ] # these should be all true or NA
# rep_df$month <- #rep(1:12, nrow(timeVar1))
#   c(
#     rep(1, nrow(timeVar1)), # I can't think of how to do this faster
#     rep(2, nrow(timeVar1)),
#     rep(3, nrow(timeVar1)),
#     rep(4, nrow(timeVar1)),
#     rep(5, nrow(timeVar1)),
#     rep(6, nrow(timeVar1)),
#     rep(7, nrow(timeVar1)),
#     rep(8, nrow(timeVar1)),
#     rep(9, nrow(timeVar1)),
#     rep(10, nrow(timeVar1)),
#     rep(11, nrow(timeVar1)),
#     rep(12, nrow(timeVar1))
#   )
#
# timeVar1 <- rep_df
# #as.data.frame(rep_df[1:100,2:9])
# #as.data.frame(rep_df[nrow(timeVar1):(nrow(timeVar1) +100),2:9])
#
# timeVar2 <- 'data/covariates/TimeVaryingCompletedPredictors/FINAL.Process.Model.Predictors.timevarying.tmin.tmax.ppt.14Apr2021.csv' %>%
#   read_csv #%>%
# #dplyr::select(state_name, fips, year, month, ndvi, precip, tmin, tmax)
# #dplyr::select(state_name, fips, year, month, precip, tmin, tmax)
# timeVar2$FIPS <- as.character(timeVar2$fips)
# timeVar2$month <- as.integer(timeVar2$month)
# #-- process time varying covariates
# timeVar <- full_join(timeVar1, timeVar2) %>%
#   mutate(c_cerealTV = c(scale(Cereals)),
#          c_fruitNutTV = c(scale(Fruit.Nuts)),
#          #c_ndviTV = c(scale(ndvi)),
#          c_precipTV = c(scale(prcp)),
#          c_tminTV = c(scale(tmin)),
#          c_tmaxTV = c(scale(tmax)))
# timeVar$endYear <- timeVar$year
# timeVar$endMonth <- as.numeric(timeVar$month)
# timeVar3 <- full_join(x=timeVar, y=timestep_df, by=c("year", "month"))
# timeVar3$STATE_NAME = timeVar3$state_name

# process_covs2 <- left_join(y=process_covs, x=timeVar3, by=c("FIPS", "STATE_NAME")) # this might be the answer!
#as.data.frame(process_covs2[1:13,])
#mean(is.na(process_covs2$c_hydroden))
# apply(process_covs2, 2, function(x) mean(is.na(x)))
# process_covs3 <- process_covs2[!is.na(process_covs2$c_hydroden),]
# apply(process_covs3, 2, function(x) mean(is.na(x)))
#process_covs4 <- process_covs3[!is.na(process_covs3$c_ndviTV), ]
# process_covs2 <- process_covs3

# for missing values, impute with state mean
mean_impute <- function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
}

# generate centered and scaled versions of these numeric variables
center_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# Deal with observation covariates, a small percentage of which are missing
get_obs_covars <- function(file){
  obs_covs <- file %>%
    read_csv %>%
    mutate(FIPS = as.character(FIPS)) %>%
    dplyr::select(-starts_with('sd'), -NAME)

  obs_covs <- obs_covs %>%
    group_by(STATE_NAME) %>%
    mutate(rural.road.density = mean_impute(rural.road.density),
           prop.pub.land = mean_impute(prop.pub.land),
           mean.ruggedness = mean_impute(mean.ruggedness),
           mean.canopy.density = mean_impute(mean.canopy.density))
  return(obs_covs)
}

# create a data set that merges surveys with county info
merge_survey_county <- function(clean_d, county_df, process_covs){
  d <- clean_d %>%
    mutate(countyname = trimws(cnty_name)) %>%
    left_join(county_df) %>%
    group_by(method) %>%
    mutate(scaled_effort = c(scale(log(1 + effort))),
           FIPS = paste0(sprintf('%02d', st_gsa_state_cd),
                         sprintf('%03d', cnty_gsa_cnty_cd)),
           Year = as.numeric(substr(start.date, 1, 4))) %>%
    ungroup |>
    mutate(fips_in_covs = FIPS %in% process_covs$FIPS) %>%
    filter(fips_in_covs)  %>%
    left_join(process_covs) |>
    filter(!is.na(c_pasture))

  # find the area of each county
  area_df <- read_excel(file.path('data', 'counties', 'area.xls')) %>%
    mutate(FIPS = STCOU)

  merged_d <- d %>%
    left_join(area_df)

  # ensure that we have area for every county
  tar_assert_true(!any(is.na(merged_d$LND010190D)))

  return(merged_d)
}

# compute start and end indices for previous surveys
start_end <- function(merged_d){
  survey_d <- merged_d %>%
    # mutate(idx = 1:n()) %>%
    mutate(idx=1:nrow(merged_d)) %>%
    dplyr::select(idx, state, countyname, FIPS, start.date, end.date, timestep, order,
                  agrp_prp_id, trap_count, method, effort, y, scaled_effort, property_area_acres,
                  order, Year,
                  starts_with('c_'), LND010190D) %>%
    arrange(FIPS, agrp_prp_id, timestep, order) %>%
    mutate(method_factor = factor(method),
           FIPS_timestep = paste(FIPS, timestep, sep = '_'),
           property_factor = as.factor(agrp_prp_id))
  # mean(is.na(survey_d$c_tree))
  return(survey_d)
}

# Create a smaller data frame with spatiotemporal covs to avoid  ----------
# redundant computations
train_test <- function(survey_d, property_area_df, timestep_df, shp, states_filter){

  ecoregions <- readOGR(file.path('data', 'ecoregions'),
                        'Cnty.lower48.EcoRegions.Level2') %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    as_tibble

  set.seed(123)

  st_d <- survey_d %>%
    filter(state %in% states_filter)

  # if("MO" %in% states_filter){
  #   st_d <- st_d |> filter(countyname != "ST CLAIR")
  # }

  st_d <- st_d |>
    distinct(FIPS_timestep, FIPS, timestep, LND010190D, agrp_prp_id,
             property_factor,
             c_hydroden, c_carnrich, c_crop, c_pasture, c_tree,
             c_evergreen, c_deciduous, c_cerealCrop, c_fruitNut,
             c_rootTuber, c_vegetablesMelons,
             c_lewis, c_hetero,
             #c_cerealTV, c_fruitNutTV, c_ndviTV, c_precipTV, c_tminTV, c_tmaxTV,
             .keep_all=TRUE
    ) %>%
    #select(c_cerealTV, c_fruitNutTV, c_ndviTV, c_precipTV, c_tminTV, c_tmaxTV) %>%
    # split into training, dev, & validation sets
    # mutate(group = sample(c('train', 'test', 'dev'),
    #                       size = n(),
    #                       #size=nrow(st_d),
    #                       replace = TRUE,
    #                       prob = c(0.7, 0.3, 0.05))) %>% # TODO: change these later
    left_join(property_area_df) %>%
    filter(property_area_km2 > 1) %>%
    left_join(ecoregions) %>%
    left_join(timestep_df) %>%
    mutate(county_index = match(FIPS, as.character(shp@data$FIPS)))

  return(st_d)
}


# Generate spatial neighbors ----------------------------------------------
spatial_neighbors <- function(st_d, timestep_df, shp, state_long){

  shp <- shp[shp@data$STATE_NAME %in% state_long,]
  shp <- shp[toupper(shp@data$NAME) %in% unique(st_d$countyname),]


  sf_use_s2(FALSE)
  nb <- poly2nb(shp, row.names = shp$FIPS)

  # generate neighborhood data for spatial prior
  listw <- nb2listw(nb, style = 'B', zero.policy = TRUE)
  B <- as(listw, 'symmetricMatrix')
  # B is suitable for building N, N_edges, node1, and node2
  # following http://mc-stan.org/users/documentation/case-studies/icar_stan.html

  # for use in dcar_normal implemented in nimble
  D <- nb2WB(nb)

  # B-splines for abundance over time
  n_year <- length(unique(substr(timestep_df$start_dates,
                                 1, 4)))
  short_basis <- splines::bs(st_d$timestep,
                             df = n_year * 2, intercept = TRUE)

  return(list(short_basis = short_basis,
              nb = nb,
              B = B,
              D = D))
}


# Generating a design matrix
design_matrix <- function(st_d){
  current.na.action <- options('na.action')
  options(na.action='na.pass')
  X <- model.matrix(~ 0 +
                      c_hydroden +
                      c_lewis +
                      c_hetero +
                      c_carnrich +
                      c_crop +
                      c_pasture +
                      #c_evergreen + # exclude evergreen, deciduous, cereal crop, pasture, tmax, tmin, precip to reduce parameters
                      #c_deciduous +
                      #c_cerealCrop +
                      #c_fruitNut+ # exclude fruitNut, rootTuber, vegetablesMelons, CerealTV, fruitNutTV because too many NAs
                      #c_rootTuber +
                      #c_vegetablesMelons +
                      #c_cerealTV +
                      #c_fruitNutTV +
                      #c_ndviTV +
                      # c_precipTV+ TODO: add back in when I get TV data
                      # c_tminTV + TODO: add back in when I get TV data
                      # c_tmaxTV + TODO: add back in when I get TV data
                      c_tree
                    ,
                    data = st_d, na.action="na.pass")

  # find the number of NAs in each column
  # apply(X, 2, function(x){mean(is.na(x))})
}


# Add observation covariates to the survey data
merge_obs_covs <- function(survey.d, obs_covs, shp){
  survey_d <- survey.d %>%
    left_join(obs_covs) %>%
    mutate(rural.road.density = mean_impute(rural.road.density),
           prop.pub.land = mean_impute(prop.pub.land),
           mean.ruggedness = mean_impute(mean.ruggedness),
           mean.canopy.density = mean_impute(mean.canopy.density),
           c_road_den = center_scale(rural.road.density),
           c_rugged = center_scale(mean.ruggedness),
           c_canopy = center_scale(mean.canopy.density),
           county_index = match(FIPS, as.character(shp@data$FIPS)))

  tar_assert_true(!any(is.na(survey_d$c_road_den)))
  tar_assert_true(!any(is.na(survey_d$county_index)))

  return(survey_d)
}

# need a function to create the following outputs from model:
# X_p
# start and end indices
# survey_idx to match space-time units
get_survey_outputs <- function(group, survey_d, st_d, property_area_df) {
  tar_assert_true(group %in% c('train', 'test', 'dev'))

  if(group == 'train') group <- c(group, "dev")

  group_d <- survey_d %>%
    filter(FIPS_timestep %in% unique(st_d$FIPS_timestep)) %>%
    left_join(property_area_df)

  get_training <- function(v, n){
    group_d |>
      select(idx, .data[[v]], y) |>
      group_by(.data[[v]]) |>
      mutate(test_county = sample(c("train", "test"),
                                  length(unique(.data[[v]])),
                                  replace = TRUE,
                                  prob = c(0.7, 0.3))) |>
      ungroup() |>
      mutate(y = ifelse(test_county == "train", y, NA)) |>
      select(idx, y)
  }

  fips_train <- get_training("FIPS") |> rename("y_fips_train" = y)
  timestep_train <- get_training("timestep") |> rename("y_timestep_train" = y)
  fips_timestep_train <- get_training("FIPS_timestep") |> rename("y_fips_timestep_train" = y)
  property_train <- get_training("agrp_prp_id") |> rename("y_property_train" = y)

  group_d <- group_d |>
    left_join(fips_train) |>
    left_join(timestep_train) |>
    left_join(fips_timestep_train) |>
    left_join(property_train) |>
    mutate(test_county = sample(c("train", "test"),
                                n(),
                                replace = TRUE,
                                prob = c(0.7, 0.3))) |>
    mutate(y_train = ifelse(test_county == "train", y, NA)) |>
    select(-test_county)

  tar_assert_true(all(group_d$property_area_km2 > 1))
  tar_assert_true(nrow(group_d) < nrow(survey_d))

  # Generate start and end indices for previous surveys ---------------------
  group_d$start <- 1
  group_d$end <- 1

  pb <- txtProgressBar(max = nrow(group_d), style = 3)
  for (i in 1:nrow(group_d)) {
    if (group_d$order[i] > 1) {
      idx <- which(group_d$FIPS == group_d$FIPS[i] &
                     group_d$agrp_prp_id == group_d$agrp_prp_id[i] &
                     group_d$timestep == group_d$timestep[i] &
                     group_d$order < group_d$order[i])
      group_d$start[i] <- idx[1]
      group_d$end[i] <- idx[length(idx)]
      tar_assert_true(identical(idx, group_d$start[i]:group_d$end[i]))
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)

  # indices to match surveys to rows in the spatiotemporal data
  group_d <- group_d %>%
    mutate(survey_idx = match(FIPS_timestep, st_d$FIPS_timestep),
           effort_per = effort / trap_count)
  tar_assert_true(!any(is.na(group_d$survey_idx)))
  tar_assert_true(!any(st_d$property_area_acres <= 0))
  tar_assert_true(!any(survey_d$effort == 0))

  # make a design matrix for detection probs
  X_p <- model.matrix(~ 0 +
                        method * c_road_den +
                        method * c_rugged +
                        method * c_canopy,
                      data = group_d)
  tar_assert_identical(nrow(X_p), nrow(group_d))

  return(list(X_p = X_p,
              group_d = group_d))
}

make_nimble_lists <- function(st_d, X, short_basis, timestep_df, train_d, survey_d, shp){

  data <- list(

    # spatiotemporal info
    area_km2 = st_d$property_area_km2,
    log_area_km2 = log(st_d$property_area_km2),
    X = X,
    X_short = short_basis$short_basis,

    # take
    y = train_d$group_d$y,
    scaled_effort = train_d$group_d$scaled_effort,
    trap_count = train_d$group_d$trap_count,
    n_trap_m1 = train_d$group_d$trap_count - 1,
    survey_area_km2 = train_d$group_d$property_area_km2,
    log_survey_area_km2 = log(train_d$group_d$property_area_km2),
    effort = train_d$group_d$effort,
    effort_per = train_d$group_d$effort_per,
    log_effort_per= log(train_d$group_d$effort_per)
  )

  constants <- list(

    # spatiotemporal info
    n_st = nrow(st_d),
    m_n = ncol(X),
    n_property = length(unique(st_d$agrp_prp_id)),
    property = st_d$agrp_prp_id |> as.factor() |> as.numeric(),
    m_short = ncol(short_basis$short_basis),
    n_timestep = nrow(timestep_df),
    timestep = st_d$timestep,
    log_pi = log(pi),

    # data for spatial indexing and ICAR priors
    n_county = st_d$countyname |> unique() |> length(),
    county_idx = st_d$countyname |> as.factor() |> as.numeric(),
    n_edges = length(short_basis$B@i),
    node1 = short_basis$B@i + 1, # add one to offset zero-based index
    node2 = short_basis$B@j + 1,
    adj = short_basis$D$adj,
    weights = short_basis$D$weights,
    num = short_basis$D$num,
    n_island =  which(card(short_basis$nb) == 0),
    island_idx = length(which(card(short_basis$nb) == 0)),
    n_survey = nrow(train_d$group_d),
    survey_idx = train_d$group_d$survey_idx,

    # take
    m_p = ncol(train_d$X_p),
    X_p = train_d$X_p,
    order = train_d$group_d$order,
    not_first_survey = as.numeric(train_d$group_d$order != 1),
    p_county_idx = train_d$group_d$countyname |> as.factor() |> as.numeric(),
    p_property_idx = train_d$group_d$agrp_prp_id |> as.factor() |> as.numeric(),
    start = train_d$group_d$start,
    end = train_d$group_d$end,
    n_method = length(levels(survey_d$method_factor)),
    method = as.numeric(train_d$group_d$method_factor),
    trap_snare_ind = as.numeric(train_d$group_d$method_factor %in% c("trap", "snare")),
    shooting_ind = as.numeric(!train_d$group_d$method_factor %in% c("trap", "snare")),
    trap_snare_idx = pmax(as.numeric(train_d$group_d$method_factor) - 3, 1)
  )


  return(list(data = data, constants = constants))

}
