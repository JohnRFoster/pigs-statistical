# srun --pty --account="iwctml" -t 0-02:00 --mem=32000 /bin/bash
# module load swset/2018.05  gcc/7.3.0 r/3.5.2-py27 r-rgdal/1.2-16-py27 r-sf/0.5-5-py27 r-rcpp/1.0.0-py27


library(plyr)
library(rgdal)
library(maptools)
#if (!require(gpclib)) install.packages("gpclib", type="source")
gpclibPermit()
library(spdep)
library(readxl)
library(assertthat)
library(tidyverse)
library(splines)
library(lubridate)


setwd("/project/iwctml/mtabak/APHIS/abundance/wild-pigs/")
outDir <- "./data/"

setwd("/Users/mikeytabak/Desktop/APHIS/abundanceModeling/wild-pigs/")

# Cleaning the MIS data ---------------------------------------------------

# read in the ecoregion data
ecoregions <- readOGR(file.path('data', 'ecoregions'), 
                      'Cnty.lower48.EcoRegions.Level2') %>%
  as.data.frame(stringsAsFactors = FALSE) %>%
  as_tibble


# Load and process in situ data -------------------------------------------
traps <- read_csv("data/insitu/feral.swine.effort.take.trap.ALL.daily.events.2021-03-25.csv") %>%
  mutate(method = "trap", 
         effort = trap.nights, 
         y = Take, 
         #ST_GSA_STATE_CD = ST_FIPS,
         ST_GSA_STATE_CD = ST_GSA_STATE_CD,
         #CNTY_GSA_CNTY_CD = CNTY_FIPS, 
         CNTY_GSA_CNTY_CD = CNTY_GSA_CNTY_CD,
         trap_count = trap.count) %>%
  dplyr::select(CNTY_NAME, ST_GSA_STATE_CD, CNTY_GSA_CNTY_CD, start.date, end.date, 
                AGRP_PRP_ID, method, effort, y, TOTAL.LAND, trap_count)

aerial <- read_csv("data/insitu/feral.swine.effort.take.aerial.ALL.daily.events.2021-03-23.csv") %>%
  mutate(method = CMP_NAME,
         effort = Flight.Hours,
         y = Take,
         #ST_GSA_STATE_CD = ST_FIPS,
         #CNTY_GSA_CNTY_CD = CNTY_FIPS, 
         trap_count = VEHICLES) %>%
  dplyr::select(CNTY_NAME, ST_GSA_STATE_CD, CNTY_GSA_CNTY_CD, 
                Start.Date, End.Date, trap_count,
                AGRP_PRP_ID, method, effort, y, TOTAL.LAND) %>%
  rename(start.date = Start.Date,
         end.date = End.Date) #%>%
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

snares <- read_csv('data/insitu/feral.swine.effort.take.snare.ALL2018-02-23.csv') %>%
  filter(!is.na(ST_NAME),
         !is.na(CNTY_NAME)) %>%
  mutate(method = 'snare',
         effort = trap.nights,
         y = Take, 
         ST_GSA_STATE_CD = ST_FIPS, 
         CNTY_GSA_CNTY_CD = CNTY_FIPS, 
         trap_count = trap.nights / event.length) %>%
  dplyr::select(CNTY_NAME, ST_GSA_STATE_CD, CNTY_GSA_CNTY_CD, start.date, end.date, 
                AGRP_PRP_ID, method, effort, y, 
                TOTAL.LAND, trap_count)

firearms <- read_csv('data/insitu/feral.swine.effort.take.firearms.ALL2018-02-24.csv') %>%
  filter(!is.na(ST_NAME), 
         !is.na(CNTY_NAME)) %>%
  mutate(method = 'firearms', 
         effort = Hunt.Hours, 
         y = Take, 
         ST_GSA_STATE_CD = ST_FIPS, 
         CNTY_GSA_CNTY_CD = CNTY_FIPS, 
         Start.date = WT_WORK_DATE, 
         End.date = WT_WORK_DATE, 
         trap_count = FIREARMS) %>%
  dplyr::select(CNTY_NAME, ST_GSA_STATE_CD, CNTY_GSA_CNTY_CD, 
                Start.date, End.date, 
                AGRP_PRP_ID, method, effort, y, TOTAL.LAND, 
                trap_count)

# set all column names to lower case
names(traps) <- tolower(names(traps))
names(firearms) <- tolower(names(firearms))
names(snares) <- tolower(names(snares))
names(aerial) <- tolower(names(aerial))


# Merge all of the take data into one data frame
insitu <- full_join(traps, firearms) %>%
  full_join(snares) %>% 
  full_join(aerial) %>%
  distinct %>%
  mutate(insitu_id = 1:n(), 
         method = tolower(method)) %>%
  filter(!is.na(y))
write_rds(insitu, paste0(outDir, 'insitu.rds'))



# resolve duplicate property values - when there are multiple values, take max
property_area_df <- distinct(insitu, agrp_prp_id, total.land) %>%
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

# merge property areas back into insitu data, filter on area
insitu <- insitu %>%
  left_join(property_area_df) %>%
  filter(!is.na(property_area_acres), 
         property_area_km2 >= 1.8, 
         effort > 0)

# fetch all fips data from census
usa_fips <- read_csv("data/fips/national_county.txt", 
                     col_names = c("state", "statefp", "countyfp", 
                                   "countyname", "classfp"), 
                     comment = '#') %>%
  mutate(st_gsa_state_cd = parse_number(statefp), 
         cnty_gsa_cnty_cd = parse_number(countyfp))

# merge fips data to insitu pig data
merged_d <- left_join(insitu, usa_fips) %>%
  arrange(countyname, start.date, end.date) %>%
  filter(st_gsa_state_cd != 15) # exclude Hawaii


# Generate time interval breaks for take data ---------------------------------------------------
min_date <- min(merged_d$start.date)
max_date <- max(merged_d$end.date)

interval <- 4 # number of weeks that comprise one 'primary period'

start_dates <- seq(min_date, max_date, by = paste(interval, "week"))
end_dates <- c(start_dates[-1] - 1, max_date)
assert_that(length(start_dates) == length(end_dates))
assert_that(min(merged_d$start.date) >= min_date)
assert_that(max(merged_d$start.date) <= max_date)

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
clean_d <- merged_d %>%
  filter(!is.na(timestep)) %>%
  arrange(start.date)

# fraction of records excluded by this rule
1 - nrow(clean_d) / nrow(merged_d)


# compute ordering based on time interval midpoints
order_df <- clean_d %>%
  left_join(property_area_df) %>%
  distinct %>%
  rowwise %>%
  mutate(midpoint = ifelse(start.date == end.date, 
                           as.numeric(start.date), 
                           as.numeric(start.date) + 
                             (as.numeric(end.date) - as.numeric(start.date)) / 2)
         ) %>%
  ungroup 

# impose stochastic ordering of events by adding jitter
# we are assuming the following order of events when the events have the same midpoint
# e.g., are on the same day:
# 1. (trap or snare), with order random
# 2. (heli or plane), with order random
# 3. hunting
order_df$jittered_midpoint <- NA
for (i in 1:nrow(order_df)) {
  if (order_df$method[i] %in% c('Trap', 'Snare')) {
    order_df$jittered_midpoint[i] <- order_df$midpoint[i] + runif(1, min = 0, max = .01)
  } else if (order_df$method[i] %in% c('Fixed Wing', 'Helicopter')) {
    order_df$jittered_midpoint[i] <- order_df$midpoint[i] + runif(1, min = .02, max = .03)
  } else {
    order_df$jittered_midpoint[i] <- order_df$midpoint[i] + runif(1, min = .04, max = .05)
  }
}

# now compute orders of events based on jittered midpoints
order_df <- order_df %>%
  ungroup %>%
  group_by(st_gsa_state_cd, cnty_name, countyfp, agrp_prp_id, timestep) %>%
  mutate(order = order(jittered_midpoint), 
         has_multi = any(order > 1), 
         any_ties = any(duplicated(jittered_midpoint)), 
         n_survey = n()) %>%
  arrange(st_gsa_state_cd, cnty_name, countyfp, agrp_prp_id, timestep, order) %>%
  ungroup()

# merge survey orders back into original data
clean_d <- clean_d %>%
  left_join(order_df)


# next, link the counties in the shapefile to the counties in the insitu data
shp <- readOGR(file.path("data", "counties"), "dtl.cnty.lower48.meters")

# Create a neighborhood adjacency matrix 
IDs <- row.names(as(shp, "data.frame"))

county_df <- shp@data %>%
  as_tibble %>%
  distinct(STATE_FIPS, CNTY_FIPS) %>%
  rename(statefp = STATE_FIPS, countyfp = CNTY_FIPS) %>%
  right_join(usa_fips) %>%
  filter(state != 'AK') %>%
  mutate(county_fips_factor = factor(countyfp), 
         countyname = gsub("County", "", x = countyname), 
         countyname = tolower(countyname), 
         countyname = trimws(countyname)) 

# create a data set that merges surveys with county info
d <- clean_d %>%
  mutate(countyname = trimws(cnty_name)) %>%
  left_join(county_df) %>%
  group_by(method) %>%
  mutate(scaled_effort = c(scale(log(1 + effort))), 
         FIPS = paste0(sprintf('%02d', st_gsa_state_cd), 
                           sprintf('%03d', cnty_gsa_cnty_cd)), 
         Year = as.numeric(substr(start.date, 1, 4))) %>%
  ungroup


## Load covariates and merge into data --------------------
spatial_covs <- "data/covariates/FINAL.Process.Model.nontime.varying.14Jan2018.csv" %>%
  read_csv

crop_covs <- 'data/covariates/FINAL.Process.Model.crop.acerage.08Jan2018.csv' %>%
  read_csv %>%
  dplyr::select(State, FIPS, Group_Name, mean.Prop.Crop) %>%
  spread(Group_Name, mean.Prop.Crop)
colnames(crop_covs) <- c("State", "FIPS", "Cereals", 
                         "Fruit.Nuts", "Other" ,"Root.Tuber", "Vegetables.Melons")

timeVar1 <- 'data/covariates/TimeVaryingCompletedPredictors/FINAL.Process.Model.time.varying.crop.acerage.22Dec2017.csv' %>%
  read_csv %>%
  spread(Group_Name, Prop.Crop)
colnames(timeVar1) <- c("State", "FIPS","year", "Crop.Acres", "Total.Acres", "Cereals", "Fruit.Nuts", "Other")
timeVar1 <- timeVar1[, -c(9,10)]
# add in a row for each month in each year
reptimes <- 12
idx_rep <- rep(1:nrow(timeVar1), reptimes)
rep_df <- timeVar1[idx_rep, ]
#timeVar1[1:10,] == rep_df[(nrow(timeVar1)+1):(nrow(timeVar1)+10), ] # these should be all true or NA
rep_df$month <- #rep(1:12, nrow(timeVar1))
  c(
    rep(1, nrow(timeVar1)), # I can't think of how to do this faster
    rep(2, nrow(timeVar1)),
    rep(3, nrow(timeVar1)),
    rep(4, nrow(timeVar1)),
    rep(5, nrow(timeVar1)),
    rep(6, nrow(timeVar1)),
    rep(7, nrow(timeVar1)),
    rep(8, nrow(timeVar1)),
    rep(9, nrow(timeVar1)),
    rep(10, nrow(timeVar1)),
    rep(11, nrow(timeVar1)),
    rep(12, nrow(timeVar1))
  )

timeVar1 <- rep_df
#as.data.frame(rep_df[1:100,2:9])
#as.data.frame(rep_df[nrow(timeVar1):(nrow(timeVar1) +100),2:9])

timeVar2 <- 'data/covariates/TimeVaryingCompletedPredictors/FINAL.Process.Model.Predictors.timevarying.tmin.tmax.ppt.14Apr2021.csv' %>%
  read_csv #%>%
  #dplyr::select(state_name, fips, year, month, ndvi, precip, tmin, tmax)
  #dplyr::select(state_name, fips, year, month, precip, tmin, tmax)
timeVar2$FIPS <- as.character(timeVar2$fips)
timeVar2$month <- as.integer(timeVar2$month)
#---figure out where NAs are coming from -- NAs are because some counties don't have 
#--- data for these covariates in csv
# crop_covs1 <- read.csv('data/covariates/FINAL.Process.Model.crop.acerage.08Jan2018.csv')[,1:6]
# crop_covs <- crop_covs1 %>%
#   dplyr::select(State, FIPS, Group_Name, mean.Prop.Crop) %>%
#   spread(Group_Name, mean.Prop.Crop)
# crops <- read.csv("data/covariates/FINAL.Process.Model.crop.acerage.08Jan2018.csv")
# colnames(crop_covs)[ apply(crop_covs, 2, anyNA) ]
# crops2 <- crop_covs[apply(crop_covs, 1, anyNA),]
# dim(crops2)

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
assert_that(!any(is.na(process_covs$Cereals)))

#-- process time varying covariates
timeVar <- full_join(timeVar1, timeVar2) %>%
  mutate(c_cerealTV = c(scale(Cereals)),
         c_fruitNutTV = c(scale(Fruit.Nuts)), 
         #c_ndviTV = c(scale(ndvi)),
         c_precipTV = c(scale(prcp)),
         c_tminTV = c(scale(tmin)),
         c_tmaxTV = c(scale(tmax)))
timeVar$endYear <- timeVar$year
timeVar$endMonth <- as.numeric(timeVar$month)

# find the timestep for each datapoint in the timeVar data
# timeVar3 <- merge(x=timeVar, y=timestep_df,
#                   by.x=c("year", "month"),
#                   by.y=c("year", "month")#, all.x=TRUE, all.y=FALSE # by removing the all=TRUE statement, I'm getting rid of covariate data that is outside of when we have MIS data
#                   )
timeVar3 <- full_join(x=timeVar, y=timestep_df, by=c("year", "month"))
timeVar3$STATE_NAME = timeVar3$state_name

# join timeVar with other covs
# I need to make all of the process_covs have the same value for each timestep
#process_covs2 <- merge(x=process_covs, y=timeVar3, by="FIPS", all.x=TRUE, all.y=FALSE)
#process_covs2 <- full_join(x=process_covs, y=timeVar3, by="FIPS")
#process_covs2 <- merge(x=process_covs, y=timeVar3, all=TRUE)
#process_covs2 <- left_join(x=process_covs, y=timeVar3, by=c("FIPS", "STATE_NAME"))
process_covs2 <- left_join(y=process_covs, x=timeVar3, by=c("FIPS", "STATE_NAME")) # this might be the answer!
#as.data.frame(process_covs2[1:13,])
#mean(is.na(process_covs2$c_hydroden))
apply(process_covs2, 2, function(x) mean(is.na(x)))
process_covs3 <- process_covs2[!is.na(process_covs2$c_hydroden),]
apply(process_covs3, 2, function(x) mean(is.na(x)))
#process_covs4 <- process_covs3[!is.na(process_covs3$c_ndviTV), ]
process_covs2 <- process_covs3
#process_covs2[process_covs2$timestep == 100, ]

# Deal with observation covariates, a small percentage of which are missing
obs_covs <- 'data/covariates/FINAL.Process.Model.Observation.Covariates.12Jan2018.csv' %>%
  read_csv %>%
  mutate(FIPS = as.character(FIPS)) %>%
  dplyr::select(-starts_with('sd'), -NAME)

# for missing values, impute with state mean
mean_impute <- function(x) {
  ifelse(is.na(x), mean(x, na.rm = TRUE), x)
}

obs_covs <- obs_covs %>%
  group_by(STATE_NAME) %>%
  mutate(rural.road.density = mean_impute(rural.road.density), 
         prop.pub.land = mean_impute(prop.pub.land), 
         mean.ruggedness = mean_impute(mean.ruggedness), 
         mean.canopy.density = mean_impute(mean.canopy.density))

# generate centered and scaled versions of these numeric variables
center_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

# verify that fips codes in the pig data match the covariate fips codes
# d <- d %>%
#   mutate(fips_in_covs = FIPS %in% process_covs$FIPS)
# do this including timevarrying
d <- d %>%
  mutate(fips_in_covs = FIPS %in% process_covs2$FIPS)

sum(!d$fips_in_covs)
# ^this many surveys occurred in counties that are not in the covariate data
# print them:
d %>%
  filter(!fips_in_covs)

# let's filter these out
d <- d %>%
  filter(fips_in_covs)

# merge process covariates into the survey data
merged_d <- d %>%
  left_join(process_covs)
#merged_d$month <- month(merged_d$end.date)
#merged_d$year <- year(merged_d$end.date)
#apply(merged_d, 2, function(x) mean(is.na(x)))

# merge all process covariates into survey data
merged_d2 <- d %>% 
  left_join(process_covs2)
#merged_d2 <- merge(d, process_covs2, all.x=TRUE, all.y=FALSE, by=c("FIPS", "timestep"))
merged_d <- merged_d2[!is.na(merged_d2$c_pasture), ]
#merged_d <- merged_d2

#assert_that(!any(is.na(merged_d$total.hydro.density)))
mean(is.na(merged_d$c_pasture))

# write this to csv so I can check
#write.csv(merged_d, "./data/covariates/merged_d.csv")

# find the area of each county
area_df <- read_excel(file.path('data', 'counties', 'area.xls')) %>%
  mutate(FIPS = STCOU)

merged_d <- merged_d %>%
  left_join(area_df)

# ensure that we have area for every county
assert_that(!any(is.na(merged_d$LND010190D))) # area in sq miles

# compute start and end indices for previous surveys
survey_d <- merged_d %>%
  # mutate(idx = 1:n()) %>%
  mutate(idx=1:nrow(merged_d)) %>%
  dplyr::select(idx, state, countyname, FIPS, start.date, end.date, timestep, order, 
         agrp_prp_id, trap_count, method, effort, y, scaled_effort, property_area_acres,
         order, Year,
         starts_with('c_'), LND010190D) %>%
  arrange(FIPS, agrp_prp_id, timestep, order) %>%
  mutate(method_factor = factor(method))
mean(is.na(survey_d$c_tree))

# Create a smaller data frame with spatiotemporal covs to avoid  ----------
# redundant computations
survey_d$FIPS_timestep <- paste(survey_d$FIPS, survey_d$timestep, sep = '_')
survey_d$property_factor <- as.factor(survey_d$agrp_prp_id)

set.seed(123)
st_d <- survey_d %>%
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
  mutate(group = sample(c('train', 'test'), 
                        size = n(), 
                        #size=nrow(st_d), 
                        replace = TRUE, 
                        prob = c(.7, .3))) %>% # TODO: change these later
  left_join(property_area_df) %>%
  filter(property_area_km2 > 1) %>%
  left_join(ecoregions) %>%
  left_join(timestep_df)



# Generate spatial neighbors ----------------------------------------------
st_d <- st_d %>%
  mutate(county_index = match(FIPS, as.character(shp@data$FIPS)))

nb <- poly2nb(shp, row.names = shp$FIPS)

island_idx <- which(card(nb) == 0)
n_island <- length(island_idx)

# generate neighborhood data for spatial prior
listw <- nb2listw(nb, style = 'B', zero.policy = TRUE)
B <- as(listw, 'symmetricMatrix')
# B is suitable for building N, N_edges, node1, and node2
# following http://mc-stan.org/users/documentation/case-studies/icar_stan.html

# B-splines for abundance over time
n_year <- length(unique(substr(timestep_df$start_dates, 
                               1, 4)))
short_basis <- splines::bs(st_d$timestep, 
                  df = n_year * 2, intercept = TRUE)
write_rds(short_basis, paste0(outDir, 'short_basis.rds'))

# Generating a design matrix
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
                    c_precipTV+
                    c_tminTV +
                    c_tmaxTV +
                    c_tree
                  , 
                  data = st_d, na.action="na.pass")

# find the number of NAs in each column
apply(X, 2, function(x){mean(is.na(x))})

# create a matrix that tells stan to ignore NAs
#X_use <- ifelse(is.na(X), 0, 1)
# then remove NAs from X
#X[is.na(X)] <- 0
# restore NA action
#options(na.action=current.na.action) # if this doesn't work (will throw error in get_survey_outputs) set: options(na.action="na.omit")

# Add observation covariates to the survey data
survey_d <- survey_d %>%
  left_join(obs_covs) %>%
  mutate(rural.road.density = mean_impute(rural.road.density), 
         prop.pub.land = mean_impute(prop.pub.land), 
         mean.ruggedness = mean_impute(mean.ruggedness), 
         mean.canopy.density = mean_impute(mean.canopy.density), 
         c_road_den = center_scale(rural.road.density),
         c_rugged = center_scale(mean.ruggedness), 
         c_canopy = center_scale(mean.canopy.density), 
         county_index = match(FIPS, as.character(shp@data$FIPS)))
  
assert_that(!any(is.na(survey_d$c_road_den)))
assert_that(!any(is.na(survey_d$county_index)))


#*** SAVE HERE FOR DIAGNOSIS
#save.image(file="./test/save_20190528.RData")

# need a function to create the following outputs from model:
# X_p
# start and end indices
# survey_idx to match space-time units
get_survey_outputs <- function(survey_d, group) {
  assert_that(group %in% c('train', 'test'))
  group_d <- survey_d %>%
    filter(FIPS_timestep %in% 
             st_d$FIPS_timestep[st_d$group == group]) %>%
    left_join(property_area_df)
  
  assert_that(all(group_d$property_area_km2 > 1))
  assert_that(nrow(group_d) < nrow(survey_d))
  
  # Generate start and end indices for previous surveys ---------------------
  group_d$start <- 0
  group_d$end <- 0
  
  pb <- txtProgressBar(max = nrow(group_d), style = 3)
  for (i in 1:nrow(group_d)) {
    if (group_d$order[i] > 1) {
      idx <- which(group_d$FIPS == group_d$FIPS[i] &
                     group_d$agrp_prp_id == group_d$agrp_prp_id[i] &
                     group_d$timestep == group_d$timestep[i] &
                     group_d$order < group_d$order[i])
      group_d$start[i] <- idx[1]
      group_d$end[i] <- idx[length(idx)]
      assert_that(identical(idx, group_d$start[i]:group_d$end[i]))
    }
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # indices to match surveys to rows in the spatiotemporal data
  group_d <- group_d %>%
    mutate(survey_idx = match(FIPS_timestep, st_d$FIPS_timestep), 
           effort_per = effort / trap_count)
  assert_that(!any(is.na(group_d$survey_idx)))
  assert_that(!any(st_d$property_area_acres <= 0))
  assert_that(!any(survey_d$effort == 0))
  
  # make a design matrix for detection probs
  X_p <- model.matrix(~ 0 +
                        method * c_road_den + 
                        method * c_rugged + 
                        method * c_canopy, 
                      data = group_d)
  assert_that(identical(nrow(X_p), nrow(group_d)))
  
  return(list(X_p = X_p,
              group_d = group_d))
}

train_d <- get_survey_outputs(survey_d, 'train')
dev_d <- get_survey_outputs(survey_d, 'test')

stan_d <- list(
  # spatiotemporal info
  n_st = nrow(st_d),
  m_n = ncol(X), 
  area_km2 = st_d$property_area_km2,
  n_property = length(unique(st_d$agrp_prp_id)), 
  property = as.numeric(st_d$property_factor), 
  X = X,
 # X_use = X_use,
  m_short = ncol(short_basis),
  X_short = short_basis,
  n_timestep = nrow(timestep_df),
  timestep = st_d$timestep,
  
  # data for spatial indexing and ICAR priors
  n_county = nrow(shp), 
  county_idx = st_d$county_index,
  n_edges = length(B@i), 
  node1 = B@i + 1, # add one to offset zero-based index
  node2 = B@j + 1,
  n_island = n_island, 
  island_idx = island_idx,
  
  # training data
  n_survey = nrow(train_d$group_d), 
  survey_idx = train_d$group_d$survey_idx, 
  y = train_d$group_d$y, 
  scaled_effort = train_d$group_d$scaled_effort,
  trap_count = train_d$group_d$trap_count,
  m_p = ncol(train_d$X_p),
  X_p = train_d$X_p, 
  order = train_d$group_d$order, 
  survey_area_km2 = train_d$group_d$property_area_km2,
  p_county_idx = train_d$group_d$county_index,
  p_property_idx = as.numeric(train_d$group_d$property_factor),
  start = train_d$group_d$start, 
  end = train_d$group_d$end, 
  n_method = length(levels(survey_d$method_factor)), 
  method = as.numeric(train_d$group_d$method_factor), 
  effort = train_d$group_d$effort,
  effort_per = train_d$group_d$effort_per,
  
  # dev data
  n_survey_dev = nrow(dev_d$group_d),
  survey_idx_dev = dev_d$group_d$survey_idx,
  y_dev = dev_d$group_d$y,
  scaled_effort_dev = dev_d$group_d$scaled_effort,
  trap_count_dev = dev_d$group_d$trap_count,
  X_p_dev = dev_d$X_p,
  order_dev = dev_d$group_d$order,
  survey_area_km2_dev = dev_d$group_d$property_area_km2,
  p_county_idx_dev = dev_d$group_d$county_index,
  p_property_idx_dev = as.numeric(dev_d$group_d$property_factor),
  start_dev = dev_d$group_d$start,
  end_dev = dev_d$group_d$end,
  method_dev = as.numeric(dev_d$group_d$method_factor),
  effort_dev = dev_d$group_d$effort,
  effort_per_dev = dev_d$group_d$effort_per
  )

# make sure the data contain no missing values
number_missing_values <- lapply(stan_d, FUN = function(x) mean(is.na(x))) %>%
  unlist %>%
  sum
no_missing_values <- number_missing_values == 0
assert_that(no_missing_values)

# make sure that all five methods made it in here
assert_that(all(levels(train_d$group_d$method_factor) %in% 
    c("firearms", "fixed wing", "helicopter", "snare", "trap")))

write_rds(merged_d, paste0(outDir, 'merged_d.rds'))
write_rds(st_d, paste0(outDir, 'st_d.rds'))
write_rds(survey_d, paste0(outDir, 'survey_d.rds'))
write_rds(train_d, paste0(outDir, 'train_d.rds'))
write_rds(dev_d, paste0(outDir, 'dev_d.rds'))
write_rds(stan_d, paste0(outDir, 'stan_d.rds'))
write_rds(process_covs, paste0(outDir, 'process_covs.rds'))
write_rds(timestep_df, paste0(outDir, 'timestep_df.rds'))

apply(X, 2, function(x) quantile(x, c(0.025,0.5,0.975)))
#-----------------------------------------------
# archive
#-----------------------------------------------
# determine if there is correlation among variables
head(process_covs)
cors <- cor(process_covs[,c(21:ncol(process_covs))], method="pearson")
write.csv(cors, "./data/covariates/processModelCovs_correlations.csv")
q
# create a matrix that tells stan to ignore NAs
vec <- c(1,2,3,NA)
mat <- matrix(sample(vec, replace=TRUE),
              nrow=10, ncol=10)
mat2 <- ifelse(is.na(mat), 0, 1)

tmp <- readRDS("stan_d.rds")
str(tmp)                 

#*** fix time varying covariates by doing it without dplyr
mergedTV <- merge(x=merged_d, y=timeVar, 
                  by=c("FIPS", "endYear", "endMonth"),
                  all=FALSE
                  #all.x=TRUE, all.y=FALSE
) # this df will have correct values for time varrying covs
mergedTV2 <- inner_join(merged_d, timeVar, by=c("FIPS", "endYear", "endMonth"))
dim(mergedTV2)

#*** Now add the time varying covariates because I have the timestep in merged_d
#*** this is where the joining is not happening properly
#*** none of these are working right
# merged_d2 <- merged_d %>%
#   #left_join(timeVar3, by=c("FIPS", "year", "month"))
# merged_d2 <- semi_join(merged_d, timeVar3, by=c("FIPS", "year", "month"))
# merged_d2 <- merge(merged_d, timeVar3, by=c("FIPS", "year", "month"))

#---- looking at test data
#---- trying to figure out why the model runs so much slower with time varrying covs
attach("/project/iwctml/mtabak/APHIS/abundance/wild-pigs/test/save_20190528.RData")
dim(process_covs)
length(unique(process_covs$FIPS))
# look at one county at a time
length(unique(merged_d$FIPS))
table(merged_d$FIPS) 
fips <- 54041 # choosing FIPS=49027, 48217, 20011, 54041
cc <- merged_d[merged_d$FIPS==fips, startsWith(colnames(merged_d), c("c_"))]
apply(cc, 2, function(x) length(unique(x)))
# there is only one value for each process cov, but multiple values for the TV covs, as expected

# now try looking at the values before merges
table(process_covs$FIPS)
c1 <- process_covs[process_covs$FIPS==fips, startsWith(colnames(process_covs), c("c_"))]
as.data.frame(cc[1,])
as.data.frame(c1)
c2 <- process_covs2[process_covs2$FIPS==fips, startsWith(colnames(process_covs2), c("c_"))]
as.data.frame(c2[1,])
apply(c2, 2, function(x) length(unique(x)))
c3 <- as.data.frame(process_covs2[process_covs2$FIPS==fips, ])
unique(c3$c_hetero)
nas <- c3[is.na(c3$c_hetero),]
nas <- process_covs2[(is.na(process_covs2$c_hetero)),]
apply(nas, 2, function(x) sum(!is.na(x)))
dim(nas)

# dims
dim(merged_d)
dim(merged_d2)
dim(X)
apply(X, 2, function(x){mean(is.na(x))})

# look at NAs
nas <- as.data.frame(merged_d[is.na(merged_d$c_pasture),])
head(nas)

# comare merged data files
Max <- readRDS("/project/iwctml/mtabak/APHIS/abundance/wild-pigs/merged_d.rds")
Mik <- readRDS("/project/iwctml/mtabak/APHIS/abundance/wild-pigs/reduced/7pars2/merged_d.rds")
all.equal(Max,Mik)
all.equal(Mik$State.x, Mik$State.y)
all.equal(Mik$STATE_NAME.x, Mik$STATE_NAME.y)
all.equal(Mik$STATE_NAME.x, Mik$STATE_NAME)
mikD <- Mik[Mik$FIPS==fips, startsWith(colnames(Mik), c("c_"))]
maxD <- Max[Max$FIPS==fips, startsWith(colnames(Max), c("c_"))]
