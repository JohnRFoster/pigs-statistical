## functions called throughout the target pipeline


## ---------------------------- Data ingest -----------------------------------

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

join_insitu <- function(traps, firearms, snares, aerial){

  # Merge all of the take data into one data frame
  insitu_all <- function(traps.data, firearms.data, snares.data, aerial.data){
    full_join(traps.data, firearms.data) %>%
      full_join(snares.data) %>%
      full_join(aerial.data) %>%
      distinct %>%
      mutate(insitu_id = 1:n(),
             method = tolower(method)) %>%
      filter(!is.na(y))
  # write_rds(insitu, paste0(outDir, 'insitu.rds'))
  }

  insitu_df <- insitu_all(traps, firearms, snares, aerial)


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

  property_area_df <- resolve_duplicate(insitu_df)

  # merge property areas back into insitu data, filter on area
  join_area_filter <- function(ins, prop){
    insitu_property <- ins %>%
      left_join(prop) %>%
      filter(!is.na(property_area_acres),
             property_area_km2 >= 1.8,
             effort > 0)
    return(insitu_property)
  }

  insitu_property <- join_area_filter(insitu_df, property_area_df)

  # fetch all fips data from census
  get_merge_fips <- function(ins.p){
    usa_fips <- read_csv("data/fips/national_county.txt",
                         col_names = c("state", "statefp", "countyfp",
                                       "countyname", "classfp"),
                         comment = '#') %>%
      mutate(st_gsa_state_cd = parse_number(statefp),
             cnty_gsa_cnty_cd = parse_number(countyfp))

    # merge fips data to insitu pig data
    merged_join <- left_join(ins.p, usa_fips) %>%
      arrange(countyname, start.date, end.date) %>%
      filter(st_gsa_state_cd != 15) # exclude Hawaii

    return(merged_join)
  }

  with_fips <- get_merge_fips(insitu_property)

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

  clean_d <- clean_merge(with_fips)

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

  order_df <- order_interval(clean_d, property_area_df)

  # impose stochastic ordering of events by adding jitter
  # we are assuming the following order of events when the events have the same midpoint
  # e.g., are on the same day:
  # 1. (trap or snare), with order random
  # 2. (heli or plane), with order random
  # 3. hunting
  order_stochastic <- function(order.df){
    order_df <- order.df
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
    return(order_df)
  }

  order_mid <- order_stochastic(order_df)

  # now compute orders of events based on jittered midpoints
  order_of_events <- function(order_df){
    order_event <- order_df %>%
      ungroup %>%
      group_by(st_gsa_state_cd, cnty_name, countyfp, agrp_prp_id, timestep) %>%
      mutate(order = order(jittered_midpoint),
             has_multi = any(order > 1),
             any_ties = any(duplicated(jittered_midpoint)),
             n_survey = n()) %>%
      arrange(st_gsa_state_cd, cnty_name, countyfp, agrp_prp_id, timestep, order) %>%
      ungroup()
    return(order_event)
  }

  order_event <- order_of_events(order_mid)

  # merge survey orders back into original data
  clean_out <- clean_d %>%
    left_join(order_df)

  return(clean_out)
}






