library(targets)
library(tidyverse)
library(lubridate)
library(nimble)

data_timestep <- read_csv("data/insitu/MIS.Effort.Take.All.Methods.Daily.Events.2021-03-25.csv") |>
  mutate(cnty_name = if_else(grepl("ST ", cnty_name), gsub("ST ", "ST. ", cnty_name), cnty_name),
         cnty_name = if_else(grepl("KERN", cnty_name), "KERN", cnty_name))

usa_fips <- read_rds("_targets/objects/usa_fips") |>
  rename(cnty_name = countyname,
         state_abr = state,
         st_gsa_state_cd = statefp) |>
  mutate(cnty_name = toupper(cnty_name),
         cnty_name = gsub(" COUNTY", "", cnty_name),
         cnty_name = gsub(" PARISH", "", cnty_name),
         st_gsa_state_cd = as.double(st_gsa_state_cd)) |>
  select(-classfp)

data_timestep <- left_join(data_timestep, usa_fips, by = c("st_gsa_state_cd", "cnty_name"))

dynamic_filter <- function(df, ts){
  df |>
    group_by(agrp_prp_id, .data[[ts]]) |>
    mutate(two_plus_takes = n() >= 2) |>
    filter(two_plus_takes) |>
    # group_by(.data[[ts]], .add = TRUE) |>
    # mutate(two_plus_timesteps = n() >= 2) |>
    # filter(two_plus_timesteps) |>
    group_by(agrp_prp_id) |>
    mutate(n_timesteps = length(unique(.data[[ts]]))) |>
    filter(n_timesteps >= 2) |>
    arrange(cnty_name, agrp_prp_id, .data[[ts]]) |>
    ungroup() |>
    mutate(p = 1:n())
}

property_filter <- data_timestep |>
  group_by(agrp_prp_id) |>
  summarise(n = sum(take)) |>
  filter(n > 0) |>
  pull(agrp_prp_id)

dynamic_4 <- data_timestep |>
  filter(agrp_prp_id %in% property_filter) |>
  dynamic_filter("timestep_4") |>
  filter(!is.na(timestep_4)) |>
  select(-timestep_1, -timestep_2, -timestep_3, -timestep_8, -timestep_12)

max_delta <- dynamic_4 |>
  filter(take > 0) |>
  select(agrp_prp_id, timestep_4) |>
  distinct() |>
  group_by(agrp_prp_id) |>
  reframe(delta = c(0, diff(timestep_4))) |>
  group_by(agrp_prp_id) |>
  filter(delta == max(delta),
         delta <= 20) |>
  pull(agrp_prp_id) |>
  unique()

data_delta <- dynamic_4 |>
  filter(agrp_prp_id %in% max_delta)

get_freqs <- function(df){
  df |>
    group_by(method) |>
    summarise(n_records = n(),
              n_take = sum(take),
              prop_records = n_records / nrow(df)) |>
    mutate(prop_take = n_take / sum(n_take)) |>
    rename(Method = method)
}

get_freqs(data_delta)


### filter data on an equal percent of methods used --------------------


percents <- data_delta |>
  select(agrp_prp_id, method) |>
  group_by(agrp_prp_id) |>
  summarise(percent_snare = length(which(method == "SNARE")) / n(),
            percent_traps = length(which(method == "TRAPS, CAGE")) / n(),
            percent_firearms = length(which(method == "FIREARMS")) / n(),
            percent_fixed = length(which(method == "FIXED WING")) / n(),
            percent_heli = length(which(method == "HELICOPTER")) / n())

shooting_filter <- percents |>
  filter(percent_firearms > 0 |
           percent_fixed > 0 |
           percent_heli > 0) |>
  filter(percent_traps < 0.4) |>
  filter(percent_snare < 0.4) |>
  pull(agrp_prp_id) |>
  unique()

fixed_filter <- percents |>
  filter(percent_fixed > 0) |>
  pull(agrp_prp_id) |>
  unique()

m_filter <- unique(c(shooting_filter, fixed_filter))

data_less_snare <- data_delta |>
  filter(agrp_prp_id %in% m_filter)

data_save <- data_less_snare |> dynamic_filter("timestep_4")
get_freqs(data_save)

data_save
write_csv(data_save, "data/insitu/sample_equal.csv")


### filter data without snares --------------------

data_no_snare <- data_delta |>
  filter(method != "SNARE") |>
  dynamic_filter("timestep_4")

get_freqs(data_no_snare)
