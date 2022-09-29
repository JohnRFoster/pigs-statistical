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
  full_join(traps, firearms) %>%
    full_join(snares) %>%
    full_join(aerial) %>%
    distinct %>%
    mutate(insitu_id = 1:n(),
           method = tolower(method)) %>%
    filter(!is.na(y))
  # write_rds(insitu, paste0(outDir, 'insitu.rds'))
}
