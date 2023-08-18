library(tidyverse)
library(data.table)
library(fst)

data <- read_csv("data/insitu/MIS_2020_timesteps.csv")
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

sample_units <- data |>
  dynamic_filter("timestep_4") |>
  # select(-alws_agrprop_id, -property.size, -n_timesteps, -p, -`...1`) |>
  arrange(agrp_prp_id, timestep_4) |>
  rename(timestep = timestep_4) |>
  filter(!is.na(timestep)) |>
  mutate(method = if_else(grepl("TRAPS", method), "TRAPS", method),
         st_gsa_state_cd = sprintf("%02d", st_gsa_state_cd),
         cnty_gsa_cnty_cd = sprintf("%03d", cnty_gsa_cnty_cd),
         fips = paste0(st_gsa_state_cd, cnty_gsa_cnty_cd))

min_date <- min(sample_units$start.date)
max_date <- max(sample_units$start.date)
interval <- 4
start_date <- seq(min_date, max_date, by = paste(interval, "week"))
end_date <- c(start_dates[-1] - 1, max_date)
timestep_df <- tibble(start_date, end_date) %>%
  mutate(PPNum = 1:n())
timestep_df$month <- month(timestep_df$end_date)
timestep_df$year <- year(timestep_df$end_date)

take <- sample_units |>
  select(fips, agrp_prp_id, take, timestep, state, cnty_name) |>
  group_by(fips, agrp_prp_id, timestep, state, cnty_name) |>
  summarise(y = sum(take)) |>
  ungroup()

n_timesteps_per_property <- take |>
  # filter(y > 0) |>
  group_by(agrp_prp_id) |>
  tally()
table(n_timesteps_per_property$n)

at_least_20_timesteps <- n_timesteps_per_property |>
  filter(n >= 20)

properties_for_eda <- at_least_20_timesteps$agrp_prp_id

take_by_property <- take |>
  filter(agrp_prp_id %in% properties_for_eda) |>
  mutate(fips = as.character(fips)) |>
  rename(PPNum = timestep)
take_by_property_m4 <- take_by_property |>
  mutate(PPNum = PPNum - 4)

conc_timesteps <- take

n_timesteps_per_county <- take |>
  # filter(y > 0) |>
  group_by(fips) |>
  tally()
table(n_timesteps_per_county$n)

PP_for_eda <- take |>
  filter(agrp_prp_id %in% properties_for_eda) |>
  pull(timestep) |>
  range()

dir_data <- "C:/Users/John.Foster/OneDrive - USDA/Desktop/fosteR/meteorology-intake/data"

data_files <- list.files(dir_data)

prcp_file <- grep("prcp.PPNum", data_files, value = TRUE)
prcp_file <- prcp_file[length(prcp_file)]
data_prcp <- read.fst(file.path(dir_data, prcp_file))

temp_file <- grep("temperature.fips", data_files, value = TRUE)
temp_file <- temp_file[length(temp_file)]
data_temp <- read.fst(file.path(dir_data, temp_file))

subset_weather <- function(df, pp_min, pp_max){
  setDT(df)
  df_county <- df[fips %in% take_by_property$fips]
  df_timestep <- df[PPNum >= pp_min & PPNum <= pp_max]
  as_tibble(df_timestep) |>
    mutate(PPNum = if_else(PPNum < 0, PPNum + 1, PPNum))
}

precip <- subset_weather(data_prcp, PP_for_eda[1]-7, PP_for_eda[2]) |>
  rename(p.total = value) |>
  arrange(fips, PPNum) |>
  group_by(fips) |>
  mutate(roll_sum_3 = RcppRoll::roll_sum(p.total, 3, fill = NA, align = "right"),
         roll_sum_4 = RcppRoll::roll_sum(p.total, 4, fill = NA, align = "right"),
         roll_sum_5 = RcppRoll::roll_sum(p.total, 5, fill = NA, align = "right"),
         roll_sum_6 = RcppRoll::roll_sum(p.total, 6, fill = NA, align = "right"),
         roll_sum_7 = RcppRoll::roll_sum(p.total, 7, fill = NA, align = "right"),
         roll_sum_8 = RcppRoll::roll_sum(p.total, 8, fill = NA, align = "right"),
         roll_sum_9 = RcppRoll::roll_sum(p.total, 9, fill = NA, align = "right"),
         roll_sum_10 = RcppRoll::roll_sum(p.total, 10, fill = NA, align = "right"),
         roll_sum_11 = RcppRoll::roll_sum(p.total, 11, fill = NA, align = "right"),
         roll_sum_12 = RcppRoll::roll_sum(p.total, 12, fill = NA, align = "right")) |>
  ungroup()

temperature <- subset_weather(data_temp, PP_for_eda[1]-7, PP_for_eda[2])
weather <- left_join(precip, temperature)
t.diff <- weather |> arrange(fips, PPNum) |> pull(t.mean) |> diff()

weather <- weather |>
  arrange(fips, PPNum) |>
  mutate(t.diff = c(0, t.diff),
         t.diff = if_else(PPNum == min(PPNum), 0, t.diff))

y.diff <- take_by_property |> arrange(fips, PPNum) |> pull(y) |> diff()
take_by_property <- take_by_property |>
  arrange(fips, PPNum) |>
  mutate(y.diff = c(0, y.diff),
         y.diff = if_else(PPNum == min(PPNum), 0, y.diff))

evi_file <- grep("evi.cnty", data_files, value = TRUE)
evi_file <- evi_file[length(evi_file)]
data_evi <- read_csv(file.path(dir_data, evi_file)) |>
  mutate(start_date = start_date + 11,
         end_date = end_date + 11)
e.diff <- data_evi |> arrange(fips, start_date) |> pull(mean) |> diff()
data_evi <- data_evi |> mutate(e.diff = c(0, e.diff))
data_evi <- left_join(data_evi, timestep_df, multiple = "all") |>
  filter(!is.na(PPNum))

take_evi <- left_join(take_by_property_m4, data_evi, multiple = "all") |>
  filter(!is.na(PPNum),
         !is.na(mean)) |>
  mutate(e.diff <- if_else(PPNum == 1, 0, e.diff))

data_evi |>
  filter(fips == "27077") |>
  ggplot() +
  aes(x = start_date, y = mean) +
  geom_point() +
  geom_line() +
  theme_bw()


take_evi |>
  group_by(fips, state) |>
  summarise(cor = cor(y, mean)) |>
  filter(!is.na(cor)) |>
  arrange(cor) |>
  ungroup() |>
  mutate(x = 1:n()) |>
  ggplot() +
  aes(x = x, y = cor, color = state) +
  geom_point() +
  geom_linerange(aes(ymin = 0, ymax = cor)) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(y = "Correlation",
       title = "EVI",
       color = "State") +
  theme_bw()  +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

center <- function(df){
  df |>
    mutate(p.total = (p.total - mean(p.total))/sd(p.total),
           roll_sum_3 = (roll_sum_3 - mean(roll_sum_3))/sd(roll_sum_3),
           roll_sum_4 = (roll_sum_4 - mean(roll_sum_4))/sd(roll_sum_4),
           roll_sum_5 = (roll_sum_5 - mean(roll_sum_5))/sd(roll_sum_5),
           roll_sum_6 = (roll_sum_6 - mean(roll_sum_6))/sd(roll_sum_6),
           roll_sum_7 = (roll_sum_7 - mean(roll_sum_7))/sd(roll_sum_7),
           roll_sum_8 = (roll_sum_8 - mean(roll_sum_8))/sd(roll_sum_8),
           roll_sum_9 = (roll_sum_9 - mean(roll_sum_9))/sd(roll_sum_9),
           roll_sum_10 = (roll_sum_10 - mean(roll_sum_10))/sd(roll_sum_10),
           roll_sum_11 = (roll_sum_11 - mean(roll_sum_11))/sd(roll_sum_11),
           roll_sum_12 = (roll_sum_12 - mean(roll_sum_12))/sd(roll_sum_12),
           t.mean = (t.mean - mean(t.mean))/sd(t.mean),
           t.min = (t.min - mean(t.min))/sd(t.min),
           t.max = (t.max - mean(t.max))/sd(t.max),
           t.diff = (t.diff - mean(t.diff))/sd(t.diff),
           t.range = (t.range - mean(t.range))/sd(t.range))
}

data_conc <- left_join(take_by_property, weather, multiple = "all") |>
  filter(!is.na(p.total)) |>
  group_by(fips) |>
  center() |>
  ungroup()

data_m4 <- left_join(take_by_property_m4, weather, multiple = "all") |>
  filter(!is.na(p.total)) |>
  group_by(fips) |>
  center() |>
  ungroup()



gg_point <- function(df, x){
  df |>
    ggplot() +
    aes(x = .data[[x]], y = y) +
    geom_point() +
    # geom_smooth() +
    theme_bw()
}

gg_point(data_conc, "p.total")
# gg_point(data_m4, "p.total")

gg_point(data_conc, "roll_sum_4")
gg_point(data_m4, "roll_sum_4")

gg_point(data_conc, "t.mean")
gg_point(data_m4, "t.mean")

gg_point(data_conc, "t.min")
gg_point(data_m4, "t.min")

gg_point(data_conc, "t.max")
gg_point(data_m4, "t.max")

gg_point(data_conc, "t.range")
gg_point(data_m4, "t.range")


gg_timeseries <- function(df, x){
  g1 <- df |>
    mutate(y = (y-mean(y))/sd(y)) |>
    ggplot() +
    aes(x = PPNum, y = y) +
    geom_line() +
    geom_line(aes(y = .data[[x]]), linetype = "dashed") +
    facet_wrap(~ agrp_prp_id, scales = "free") +
    theme_bw()
  g1
}
#
# gg_timeseries(data_conc, "p.total")
# gg_timeseries(data_conc, "roll_sum_4")
# gg_timeseries(data_conc, "t.mean")
#
# gg_timeseries(data_m4, "p.total")
# gg_timeseries(data_m4, "roll_sum_4")
#
# cor(data_conc$y, data_conc$p.total)
# cor(data_m4$y, data_m4$p.total)
#
# cor(data_conc$y, data_conc$t.mean)
# cor(data_m4$y, data_m4$t.mean)
#
# cor(data_conc$y, data_conc$t.min)
# cor(data_m4$y, data_m4$t.min)
#
# cor(data_conc$y, data_conc$t.max)
# cor(data_m4$y, data_m4$t.max)
#
# cor(data_conc$y, data_conc$t.range)
# cor(data_m4$y, data_m4$t.range)
#
# cor(data_conc$y, data_conc$t.diff)
#
#
# cor(data_conc$y, data_conc$roll_sum_4)
#
# data_conc |>
#   group_by(agrp_prp_id) |>
#   summarise(c = cor(y, roll_sum_4))

cor_summary <- function(df){
  df |>
    summarise(cor_t.mean = cor(y, t.mean),
              cor_t.max = cor(y, t.max),
              cor_t.min = cor(y, t.min),
              cor_t.range = cor(y, t.range),
              cor_roll_sum_3 = cor(y, roll_sum_3),
              cor_roll_sum_4 = cor(y, roll_sum_4),
              cor_roll_sum_5 = cor(y, roll_sum_5),
              cor_roll_sum_6 = cor(y, roll_sum_6),
              cor_roll_sum_7 = cor(y, roll_sum_7),
              cor_roll_sum_8 = cor(y, roll_sum_8),
              cor_roll_sum_9 = cor(y, roll_sum_9),
              cor_roll_sum_10 = cor(y, roll_sum_10),
              cor_roll_sum_11 = cor(y, roll_sum_11),
              cor_roll_sum_12 = cor(y, roll_sum_12),
              cor_p.total = cor(y, p.total))
}
cor_summary_diff <- function(df){
  df |>
    summarise(cor_t.mean_diff = cor(y.diff, t.mean),
              cor_t.max_diff = cor(y.diff, t.max),
              cor_t.min_diff = cor(y.diff, t.min),
              cor_t.range_diff = cor(y.diff, t.range),
              cor_roll_sum_3_diff = cor(y.diff, roll_sum_3),
              cor_roll_sum_4_diff = cor(y.diff, roll_sum_4),
              cor_roll_sum_5_diff = cor(y.diff, roll_sum_5),
              cor_roll_sum_6_diff = cor(y.diff, roll_sum_6),
              cor_roll_sum_7_diff = cor(y.diff, roll_sum_7),
              cor_roll_sum_8_diff = cor(y.diff, roll_sum_8),
              cor_roll_sum_9_diff = cor(y.diff, roll_sum_9),
              cor_roll_sum_10_diff = cor(y.diff, roll_sum_10),
              cor_roll_sum_11_diff = cor(y.diff, roll_sum_11),
              cor_roll_sum_12_diff = cor(y.diff, roll_sum_12),
              cor_p.total_diff = cor(y.diff, p.total))
}

cors <- data_conc |>
  group_by(agrp_prp_id, state, cnty_name) |>
  cor_summary() |>
  # mutate(across(cor_t.mean:cor_p.total, ~ if_else(abs(.x) < 0.1, 0, .x)))
  pivot_longer(cols = -c(agrp_prp_id, state, cnty_name),
               values_to = "Corr",
               names_to = "Weather") |>
  ungroup() |>
  filter(!is.na(Corr))
cors_diff <- data_conc |>
  group_by(agrp_prp_id, state, cnty_name) |>
  cor_summary_diff() |>
  # mutate(across(cor_t.mean:cor_p.total, ~ if_else(abs(.x) < 0.1, 0, .x)))
  pivot_longer(cols = -c(agrp_prp_id, state, cnty_name),
               values_to = "Corr",
               names_to = "Weather") |>
  ungroup() |>
  filter(!is.na(Corr))

cors <- bind_rows(cors, cors_diff)

cors_4 <- data_m4 |>
  group_by(agrp_prp_id, state, cnty_name) |>
  cor_summary() |>
  # mutate(across(cor_t.mean:cor_p.total, ~ if_else(abs(.x) < 0.1, 0, .x)))
  pivot_longer(cols = -c(agrp_prp_id, state, cnty_name),
               values_to = "Corr",
               names_to = "Weather") |>
  ungroup() |>
  filter(!is.na(Corr))

plot_cor_lollipop <- function(df, x){
  df |>
    filter(Weather == x) |>
    arrange(Corr) |>
    # group_by(cnty_name) |>
    mutate(x = 1:n()) |>
    ggplot() +
    aes(x = x, y = Corr, color = state) +
    geom_point() +
    geom_linerange(aes(ymin = 0, ymax = Corr)) +
    geom_hline(yintercept = 0) +
    # facet_wrap(~cnty_name) +
    coord_cartesian(ylim = c(-1, 1)) +
    labs(y = "Correlation",
         x = "",
         color = "State") +
    theme_bw()  +
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          legend.position = "none")
}

library(ggpubr)

g1 <- plot_cor_lollipop(cors, "cor_t.range_diff") + labs(title = "Concurrent mean temperature")
g2 <- plot_cor_lollipop(cors_4, "cor_t.mean") + labs(title = "Lag 4 mean temperature")
ggarrange(g1, g2, common.legend = TRUE, legend = "none")

g1 <- plot_cor_lollipop(cors, "cor_t.max_diff") + labs(title = "Concurrent max temperature")
g2 <- plot_cor_lollipop(cors_4, "cor_t.max") + labs(title = "Lag 4 max temperature")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_t.min_diff") + labs(title = "Concurrent min temperature")
g2 <- plot_cor_lollipop(cors_4, "cor_t.min") + labs(title = "Lag 4 min temperature")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_p.total") + labs(title = "Concurrent total precipitation")
g2 <- plot_cor_lollipop(cors_4, "cor_p.total") + labs(title = "Lag 4 total precipitation")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_roll_sum_3") + labs(title = "Total precipitation previous 3 months")
g2 <- plot_cor_lollipop(cors_4, "cor_roll_sum_3") + labs(title = "Lag 4 Total precipitation previous 3 months")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_roll_sum_4") + labs(title = "Total precipitation previous 4 months")
g2 <- plot_cor_lollipop(cors_4, "cor_roll_sum_4") + labs(title = "Lag 4 Total precipitation previous 4 months")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_roll_sum_5") + labs(title = "Total precipitation previous 5 months")
g2 <- plot_cor_lollipop(cors_4, "cor_roll_sum_5") + labs(title = "Lag 4 Total precipitation previous 5 months")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_roll_sum_6") + labs(title = "Total precipitation previous 6 months")
g2 <- plot_cor_lollipop(cors_4, "cor_roll_sum_6") + labs(title = "Lag 4 Total precipitation previous 6 months")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_roll_sum_7") + labs(title = "Total precipitation previous 7 months")
g2 <- plot_cor_lollipop(cors_4, "cor_roll_sum_7") + labs(title = "Lag 4 Total precipitation previous 7 months")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_roll_sum_8") + labs(title = "Total precipitation previous 8 months")
g2 <- plot_cor_lollipop(cors_4, "cor_roll_sum_8") + labs(title = "Lag 4 Total precipitation previous 8 months")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_roll_sum_9") + labs(title = "Total precipitation previous 9 months")
g2 <- plot_cor_lollipop(cors_4, "cor_roll_sum_9") + labs(title = "Lag 4 Total precipitation previous 9 months")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_roll_sum_10") + labs(title = "Total precipitation previous 10 months")
g2 <- plot_cor_lollipop(cors_4, "cor_roll_sum_10") + labs(title = "Lag 4 Total precipitation previous 10 months")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_roll_sum_11") + labs(title = "Total precipitation previous 11 months")
g2 <- plot_cor_lollipop(cors_4, "cor_roll_sum_11") + labs(title = "Lag 4 Total precipitation previous 11 months")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

g1 <- plot_cor_lollipop(cors, "cor_roll_sum_12") + labs(title = "Total precipitation previous 12 months")
g2 <- plot_cor_lollipop(cors_4, "cor_roll_sum_12") + labs(title = "Lag 4 Total precipitation previous 12 months")
ggarrange(g1, g2, common.legend = TRUE, legend = "bottom")

data_m4 |>
  group_by(agrp_prp_id) |>
  summarise(cor_t.mean = cor(y, t.mean),
            cor_t.max = cor(y, t.max),
            cor_t.min = cor(y, t.min),
            cor_roll_sum_4 = cor(y, roll_sum_4),
            cor_p.total = cor(y, p.total))





daily.mean.temp <- fst::read_fst("../meteorology-intake/data/temp.daily.means.fips.2023-05-12.fst")
setDT(daily.mean.temp)
temp <- daily.mean.temp[fips %in% take_by_property$fips]
temp <- unique(temp)
temp <- temp[,PPNum := as.numeric(PPNum)]
temp <- temp[PPNum >= -2]
temp <- temp[,date := lubridate::ymd(date)]
temp <- temp[,year := lubridate::year(date)]
temp <- temp[,value := value - 273.15]
temp <- temp[,gd := if_else(value > 0, value-0, 0)]
temp <- temp[order(fips, date)]
temp <- temp[,growing_deg := cumsum(gd), by = c("fips", "year")]
temp <- temp[date %in% weather$end_dates]
temp <- as_tibble(temp) |>
  select(fips, PPNum, growing_deg)





gd <- left_join(take_by_property, temp) |>
  filter(!is.na(growing_deg))

gd_cor <- gd |>
  group_by(agrp_prp_id, state, cnty_name) |>
  summarise(cor = cor(y, growing_deg)) |>
  filter(!is.na(cor)) |>
  ungroup()

gd_cor |>
  arrange(cor) |>
  mutate(x = 1:n()) |>
  ggplot() +
  aes(x = x, y = cor, color = state) +
  geom_point() +
  geom_linerange(aes(ymin = 0, ymax = cor)) +
  geom_hline(yintercept = 0) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(y = "Correlation",
       title = "Cumulative growing degrees w.r.t 10C",
       color = "State") +
  theme_bw()  +
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_blank())





