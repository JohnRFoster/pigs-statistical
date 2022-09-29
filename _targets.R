# Created by use_targets().
# Follow the comments below to fill in this target script.
# Then follow the manual to check and run the pipeline:
#   https://books.ropensci.org/targets/walkthrough.html#inspect-the-pipeline # nolint

# Load packages required to define the pipeline:
library(targets)
# library(tarchetypes) # Load other packages as needed. # nolint

# Set target options:
tar_option_set(
  packages = c("tidyverse", "lubridate"), # packages that your targets need to run
  format = "rds" # default storage format
  # Set other options as needed.
)

# tar_make_clustermq() configuration (okay to leave alone):
options(clustermq.scheduler = "multiprocess")

# tar_make_future() configuration (okay to leave alone):
# Install packages {{future}}, {{future.callr}}, and {{future.batchtools}} to allow use_targets() to configure tar_make_future() options.

# Run the R scripts in the R/ folder with your custom functions:
# tar_source() # will run all scripts in a specified directory
source("R/functions_data.R") # Source other scripts as needed. # nolint

# Replace the target list below with your own:
list(

  ## get data
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
    name = insitu,
    command = join_insitu(data.traps, data.firearms, data.snares, data.aerial))

)
