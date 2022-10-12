# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(

  # packages that your targets need to run
  packages = c("tidyverse",
               "lubridate",
               "spdep",
               "readxl",
               "spatialreg",
               "rgdal"),

  # default storage format
  format = "rds"

  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source() # will run all scripts in a specified directory
source("R/01_setup.R")
source("R/functions_data.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(

  tar_target(
    name = file.timestep,
    command = "data/timestep_df.rds",
    format = "file"
    ),
  tar_target(
    name = timestep_df,
    command = get_timestep_df(file.timestep)),

  tar_target(
    name = shp,
    command = get_shp()
  ),

  tar_target(
    name = file.traps,
    command = "data/insitu/feral.swine.effort.take.trap.ALL.daily.events.2021-03-25.csv",
    format = "file"),
  tar_target(
    name = data.traps,
    command = get_traps(file.traps)),

  tar_target(
    name = file.aerial,
    command = "data/insitu/feral.swine.effort.take.aerial.ALL.daily.events.2021-03-23.csv",
    format = "file"),
  tar_target(
    name = data.aerial,
    command = get_aerial(file.aerial)),

  tar_target(
    name = file.snares,
    command = "data/insitu/feral.swine.effort.take.snare.ALL.daily.events.2021-03-24.csv",
    format = "file"),
  tar_target(
    name = data.snares,
    command = get_snares(file.snares)),

  tar_target(
    name = file.firearms,
    command = "data/insitu/feral.swine.effort.take.firearms.ALL.daily.events.2021-03-24.csv",
    format = "file"),
  tar_target(
    name = data.firearms,
    command = get_firearms(file.firearms)),

  tar_target(
    name = insitu_df,
    command = insitu_all(data.traps, data.firearms, data.snares, data.aerial)),

  tar_target(
    name = property_area_df,
    command = resolve_duplicate(insitu_df)
  ),

  tar_target(
    name = insitu_property,
    command = join_area_filter(insitu_df, property_area_df)
  ),

  tar_target(
    name = file.fips,
    command = "data/fips/national_county.txt",
    format = "file"),
  tar_target(
    name = usa_fips,
    command = get_fips(file.fips)),

  tar_target(
    name = with_fips,
    command = merge_fips(insitu_property, usa_fips)
  ),

  tar_target(
    name = clean_d,
    command = clean_merge(with_fips)
  ),

  tar_target(
    name = order_df,
    command = order_interval(clean_d, property_area_df)
  ),

  tar_target(
    name = order_mid,
    command = order_stochastic(order_df)
  ),

  tar_target(
    name = order_event,
    command = order_of_events(order_mid, clean_d)
  ),

  tar_target(
    name = county_df,
    command = get_county_shp(usa_fips, shp)
  ),

  tar_target(
    name = file.spatial,
    command = "data/covariates/FINAL.Process.Model.nontime.varying.14Jan2018.csv",
    format = "file"),
  tar_target(
    name = file.crop,
    command = "data/covariates/FINAL.Process.Model.crop.acerage.08Jan2018.csv",
    format = "file"),
  tar_target(
    name = process_covs,
    command = get_center_nontime_covar(file.spatial, file.crop)),

  tar_target(
    name = file.observation.covar,
    command = 'data/covariates/FINAL.Process.Model.Observation.Covariates.12Jan2018.csv',
    format = "file"),
  tar_target(
    name = obs_covs,
    command = get_obs_covars(file.observation.covar)
  ),

  tar_target(
    name = merged_d,
    command = merge_survey_county(order_event, county_df, process_covs)
  ),

  tar_target(
    name = survey_d,
    command = start_end(merged_d)
  ),

  tar_target(
    name = st_d,
    command = train_test(survey_d, property_area_df, timestep_df, shp)
  ),

  tar_target(
    name = short_basis,
    command = spatial_neighbors(st_d, timestep_df, shp)
  ),

  tar_target(
    name = X,
    command = design_matrix(st_d)
  ),

  tar_target(
    name = survey_obs,
    command = merge_obs_covs(survey_d, obs_covs, shp)
  ),

  tar_target(
    name = group_d,
    command = get_survey_outputs(run_group, survey_obs, st_d, property_area_df)
  ),

  tar_target(
    name = nimble_lists,
    command = make_nimble_lists(st_d, X, short_basis, timestep_df, group_d,
                                survey_d, shp)
  )
)
