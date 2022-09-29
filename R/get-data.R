
# Get data for MIS analysis -----------------------------------------------

library(googledrive)
library(assertthat)
# provide API access to Google Drive
drive_auth(use_oob = TRUE)

#* need to set directory to repro location
# setwd("~/Desktop/APHIS/abundanceModeling/wild-pigs")

#### let's try the workflow without manually setting the wd


# Get county-level shapefile for USA --------------------------------------

# set path to where shapefiles are stored
shp_path <- file.path('data', 'counties', 'dtl.cnty.lower48.meters.shp')



# check if we have data
shp_missing <- !file.exists(shp_path)


if (shp_missing) {
  zip_loc_on_gdrive <- drive_find(pattern = 'Detailed.U.S.Counties')
  assert_that(nrow(zip_loc_on_gdrive) == 1)
  zip_path <- gsub(pattern = 'shp', x = shp_path, replacement = 'zip')
  drive_download(zip_loc_on_gdrive, path = zip_path)
  unzip(zip_path, exdir = file.path('data', 'counties'))
}



# Get ecoregion shapefile ------------------------------------------------
ecoregion_path <- file.path('data',
                            'ecoregions',
                            'Cnty.lower48.EcoRegions.Level2.shp')
ecoregions_missing <- !file.exists(ecoregion_path)
if (ecoregions_missing) {
  zip_loc_on_gdrive <- drive_find(pattern = 'Cnty.lower48.EcoRegions')
  assert_that(nrow(zip_loc_on_gdrive) == 1)
  zip_path <- gsub(pattern = 'shp', x = ecoregion_path, replacement = 'zip')
  drive_download(zip_loc_on_gdrive, path = zip_path)
  unzip(zip_path, exdir = file.path('data', 'ecoregions'))
}




# Get insitu data ---------------------------------------------------------

most_recent_insitu <- c(
  'feral.swine.effort.take.trap.ALL.daily.events.2021-03-25.csv',
  'feral.swine.effort.take.aerial.ALL.daily.events.2021-03-25.csv',
  'feral.swine.effort.take.snare.ALL.daily.events.2021-03-25.csv',
  'feral.swine.effort.take.firearms.ALL.daily.events.2021-03-25.csv'
)


insitu_prefix <- file.path('data', 'insitu')


maybe_get_csv <- function(filename, prefix) {
  insitu_path <- file.path(prefix, filename)
  if (!file.exists(insitu_path)) {
    drive_info <- drive_find(pattern = filename)
    drive_download(drive_info, path = insitu_path)
  }
}

# download insitu data if not already present
lapply(most_recent_insitu, maybe_get_csv, prefix = insitu_prefix)




# Fetch covariate data ----------------------------------------------------

most_recent_covariates <- c(
  'FINAL.Process.Model.nontime.varying.14Jan2018.csv',
  'FINAL.Process.Model.Observation.Covariates.12Jan2018.csv',
  'FINAL.Process.Model.crop.acerage.08Jan2018.csv'
)

cov_prefix <- file.path('data', 'covariates')

# download covariate data if not already present
lapply(most_recent_covariates, maybe_get_csv, prefix = cov_prefix)




# Official county areas ---------------------------------------------------

area_path <- file.path('data', 'counties', 'area.xls')
if (!file.exists(area_path)) {
  download.file("http://www2.census.gov/prod2/statcomp/usac/excel/LND01.xls",
                destfile = area_path)
}
