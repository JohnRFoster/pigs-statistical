# --------------------------------------------------------------------
#
# Generate psuedo-data for removal model simulations ----
#
# John Foster
#
# --------------------------------------------------------------------


library(tidyverse)

df <- read_rds("data/insitu/MIS_4weekPP.rds")
df <- df |>
  select(-method) |>
  rename(property = agrp_prp_id,
         method = Method)

# need a lookup table for property IDs and how may methods they employ ----
n_method_lookup <- df |>
  select(property, method) |>
  distinct() |>
  count(property)

# get the proportion of properties that use n methods
# for determining the number of properties to simulate
n_rel <- n_method_lookup |>
  count(n, name = "n_sum") |>
  mutate(rel_prop = n_sum / sum(n_sum),
         n_simulate = ceiling(rel_prop * 100))


n_pp <- 75                         # the number of primary periods to simulate
n_prop <- sum(n_rel$n_simulate)    # total number of properties to simulate

# function to get property IDs given number of methods used ----
get_properties <- function(n_method){
  n_method_lookup |>
    filter(n == n_method) |>
    pull(property)
}

# -----------------------------------------------------------------
# 1-method properties ----
# -----------------------------------------------------------------

n_one_method <- n_rel$n_simulate[1]
one_method_props <- get_properties(1)

# list of things we need to keep track of
property_attributes <- list(
  num = NULL,
  method_1 = NULL,
  area = NULL,
  effort = NULL
)

# initiate list with the number of properties we need ----
properties_one <- rep(list(property_attributes), n_one_method)

## function to place vectors into each property attributes' list ----
place_vec <- function(ls, n, where, what){
  require(purrr)
  ls |> map2(1:n, \(x, y) assign_in(x, where, what[y]))
}

properties_one <- place_vec(properties_one, n_one_method, "num", 1:n_one_method)

## relative proportion of properties that use a specific method ----
## for sampling from to assign 1-method properties a method
one_method_prob <- df |>
  filter(property %in% one_method_props) |>
  select(property, method) |>
  distinct() |>
  count(method) |>
  mutate(prob = n / sum(n))

## assign methods to 1-method properties ----
sample_one_method <- sample(one_method_prob$method,
                            n_one_method,
                            one_method_prob$prob,
                            replace = TRUE)

properties_one <- place_vec(properties_one, n_one_method, "method_1", sample_one_method)

## function to sample property area from a given set of properties and methods ----
get_area <- function(dat, m){
  dat |>
    filter(method == m) |>
    pull(property.size) |>
    sample(1)
}

one_method_areas <- df |>
  filter(property %in% one_method_props) |>
  select(property, method, property.size) |>
  distinct()

## assign areas to 1-method properties ----
areas <- sample_one_method |>
  map_vec(\(x) get_area(one_method_areas, x))

properties_one <- place_vec(properties_one, n_one_method, "area", round(areas, 2))

## now we need to determine which PPs are sampled for each property ----
### start with determining the first PP ----

### function to assign a stating PP given method ----
start_pp <- function() sample.int(65, n_one_method, replace = TRUE)
start <- start_pp()

### lookup table for return intervals by method for single method-properties ----
one_method_return <- df |>
  filter(property %in% one_method_props) |>
  select(property, timestep, method) |>
  distinct() |>
  group_by(property, method) |>
  mutate(delta = c(0, diff(timestep)),
         single_occurance = if_else(n() == 1, 1, 0),
         first_occurance = if_else(delta == 0 & single_occurance == 0, 1, 0))

#### function to assign sampling occasions ----
get_sample_occasions <- function(return_df, m, s, max_pp){
  sample_occasions <- s
  end_pp <- s
  while(end_pp <= max_pp){
    interval <- return_df |>
      filter(method %in% m,
             first_occurance == 0) |>
      pull(delta) |>
      sample(1)
    end_pp <- end_pp + interval
    sample_occasions <- c(sample_occasions, end_pp)

    # need to make sure we get at least two sample occasions
    if(length(sample_occasions) == 2 & end_pp > n_pp){
      sample_occasions <- sample_occasions[1]
      end_pp <- sample_occasions
    }
  }
  sample_occasions[sample_occasions <= n_pp] |> as.vector()
}

### lookup table for the number of repeat events in a PP by method for single-method properties ----
n_reps_single_method <- df |>
  filter(property %in% one_method_props) |>
  group_by(property, timestep, method) |>
  count() |>
  ungroup()

#### function to sample number of removal events in a PP given method
get_reps <- function(reps_df, m, size){
  reps_df |>
    filter(method == m) |>
    pull(n) |>
    sample(size, replace = TRUE)
}

### sample the number of observations and reps, place in properties list
pb <- txtProgressBar(min = 1, max = n_one_method, style = 3)
for(i in seq_len(n_one_method)){
  sample_occasions <- get_sample_occasions(one_method_return, sample_one_method[i], start[i], n_pp)
  n_reps <- get_reps(n_reps_single_method, sample_one_method[i], length(sample_occasions))

  effort <- tibble(sample_occasions = sample_occasions, n_reps = n_reps)

  testthat::expect(effort$sample_occasions[nrow(effort)] <= n_pp,
                   paste("Last PP exceeds boundary for 1-method property", i))

  properties_one <- assign_in(properties_one, list(i, "effort"), effort)
  setTxtProgressBar(pb, i)
}
close(pb)

# -----------------------------------------------------------------
# 2-method properties ----
# -----------------------------------------------------------------

n_two_method <- n_rel$n_simulate[2]
two_method_props <- get_properties(2)

# initiate list with the number of properties we need ----
# list of things we need to keep track of
property_attributes <- list(
  num = NULL,
  method_1 = NULL,
  method_2 = NULL,
  area = NULL,
  effort = NULL
)
properties_two <- rep(list(property_attributes), n_two_method)

properties_two <- place_vec(properties_two, n_two_method, "num", 1:n_two_method)

## count of method pairs by first method used then second ----
## for sampling from to assign 2-method properties their methods
m_vec <- c("Traps", "Snares", "Firearms", "Fixed wing", "Helicopter")

get_method_rel_freq <- function(prop_vec){
  temp <- df |>
    filter(property %in% prop_vec) |>
    select(property, method) |>
    distinct() |>
    group_by(property) |>
    mutate(m = paste0("m_", 1:n())) |>
    ungroup() |>
    pivot_wider(names_from = m,
                values_from = method) |>
    select(starts_with("m_"))

  m <- colnames(temp)

  temp |> count(!!!syms(m)) |>
    mutate(prob = n / sum(n))
}

### the probability of using a given method first for 2-method properties ----
two_method_prob <- get_method_rel_freq(two_method_props)

### sample the index of rows in probability table given probability column
pairs_sample <- sample(1:nrow(two_method_prob), n_two_method, prob = two_method_prob$prob, replace = TRUE)

### use index samples to select method pairs
sample_two_method <- two_method_prob |>
  slice(pairs_sample) |>
  select(starts_with("m_"))

m1 <- sample_two_method |> pull(m_1)
m2 <- sample_two_method |> pull(m_2)

## place the pair of methods in the 2-method data list ----
properties_two <- place_vec(properties_two, n_two_method, "method_1", m1)
properties_two <- place_vec(properties_two, n_two_method, "method_2", m2)

## function to sample property area from a given set of properties and methods ----
subset_area_n_method <- function(prop_vec){
  df |>
    filter(property %in% prop_vec) |>
    select(property, method, property.size) |>
    distinct() |>
    group_by(property) |>
    mutate(n = paste0("n_", 1:n())) |>
    ungroup() |>
    pivot_wider(names_from = n,
                values_from = method)
}

get_area_n <- function(area_df, m){

  area_df |>
    pivot_longer(cols = starts_with("n_"),
                 names_to = "n",
                 values_to = "method") |>
    group_by(property) |>
    filter(all(method %in% m)) |>
    ungroup() |>
    select(property, property.size) |>
    distinct() |>
    pull(property.size) |>
    sample(1)

}

## assign areas to 2-method properties ----
two_method_areas <- subset_area_n_method(two_method_props)
areas <- 1:n_two_method |>
  map(\(x) get_area_n(two_method_areas, c(m1[x], m2[x]))) |>
  list_c()

properties_two <- place_vec(properties_two, n_two_method, "area", round(areas, 2))

## joint and single return intervals ----
n_return_intervals <- function(prop_vec, n){
  temp <- df |>
    filter(property %in% prop_vec) |>
    select(property, timestep, method) |>
    distinct() |>
    pivot_wider(names_from = method,
                values_from = method) |>
    unite(method, m_vec, sep = "&", na.rm = TRUE) |>
    group_by(property, method) |>
    mutate(delta = c(0, diff(timestep))) |>
    # filter(delta > 0) |>
    ungroup() |>
    separate(method, paste0("method_", 1:n), sep = "&") |>
    suppressWarnings() # warns about not enough pieces, fill with NA is what we want

  n <- apply(temp[,paste0("method_", 1:n)], 1, function(x) sum(!is.na(x)))
  temp |> mutate(n = n)
}

all_return_intervals <- n_return_intervals(two_method_props, 2)

## function to create a return interval df given two methods ----
get_joint_return <- function(return_df, m){

  nm <- length(m)

  p <- return_df |>
    pivot_longer(cols = starts_with("method"),
                 names_to = "Mn",
                 values_to = "method",
                 values_drop_na = TRUE) |>
    select(property, method) |>
    distinct() |>
    group_by(property) |>
    filter(all(method %in% m)) |>
    ungroup() |>
    pull(property) |>
    unique()

  joint_delta <- return_df |>
    filter(property %in% p,
           delta > 0)

  # there are some combinations (i.e. snares and fixed wing) where all delta = 0
  # which means this combination only appears once in the time series for each property
  # to included these combinations let's just randomly assign a return interval
  if(all(is.na(joint_delta$method_2))){

    adjustment <- return_df |>
      filter(property %in% p,
             !is.na(method_2)) |>
      mutate(delta = sample.int(50, n()))

    return(bind_rows(joint_delta, adjustment))

  } else {

    return(joint_delta)

  }
}

## function to generate sample occasions of single or joint use of methods
get_sample_occasions_two <- function(joint_return_df, m1, m2, max_pp){
  start <- start_pp()

  obs <- joint_return_df |>
    slice(sample.int(nrow(joint_return_df), 1))
  effort <- obs |>
    select(starts_with("method")) |>
    mutate(sample_occasions = start)

  create_effort_df <- function(start_effort, start_pp){
    effort <- start_effort
    end_pp <- start_pp
    for(xx in start_pp:max_pp){
      obs <- joint_return_df |>
        slice(sample.int(nrow(joint_return_df), 1))

      interval <- obs |> pull(delta)
      end_pp <- end_pp + interval

      m <- obs |>
        select(starts_with("method")) |>
        mutate(sample_occasions = end_pp)
      effort <- bind_rows(effort, m)

      # need to make sure we get at least two sample occasions
      if(nrow(effort) == 2 & end_pp >= max_pp){
        effort <- effort |> slice(1)
        end_pp <- effort |> pull(sample_occasions)
      }

      if(end_pp > max_pp){
        effort <- effort |> slice(-nrow(effort))
      }

    }
    effort
  }

  effort_sample <- create_effort_df(effort, start)

  # need to make sure at least one occasion has two methods
  while(all(is.na(effort_sample$method_2))){
    effort_sample <- create_effort_df(effort, start)
  }
  effort_sample
}

get_n_reps_joint <- function(prop_vec){
  df |>
    filter(property %in% prop_vec) |>
    group_by(property, timestep, method) |>
    count() |>
    ungroup() |>
    pivot_wider(names_from = method,
                values_from = n)
}

n_reps_two_method <- get_n_reps_joint(two_method_props)

#### function to sample number of removal events in a PP given a pair of methods ----
get_reps_two <- function(m1, m2, ...){

  m <- c(m1, m2, ...)
  m <- m[complete.cases(m)]

  df <- n_reps_two_method |>
    select(all_of(m)) |>
    filter(if_all(all_of(m), ~ !is.na(.)))

  if(ncol(df) == 1){
    df <- df |> filter(if_all(all_of(m), ~ . >= 2))
    reps <- c(apply(df, 2, sample, 1), NA)
  } else {
    reps <- apply(df, 2, sample, 1)
  }
  tibble(n_reps_1 = reps[1], n_reps_2 = reps[2])
}


### sample the number of observations and reps, place in properties list
pb <- txtProgressBar(min = 1, max = n_two_method, style = 3)
for(i in seq_len(n_two_method)){
  sample_occasions <- get_joint_return(all_return_intervals, c(m1[i], m2[i])) |>
    get_sample_occasions_two(m1[i], m2[i], n_pp)

  m1_i <- sample_occasions |> pull(method_1)
  m2_i <- sample_occasions |> pull(method_2)

  n_reps <- map2(m1_i, m2_i, get_reps_two) |> list_rbind()
  effort <- bind_cols(sample_occasions, n_reps)

  testthat::expect(effort$sample_occasions[nrow(effort)] <= n_pp,
                   paste("Last PP exceeds boundary for 2-method property", i))
  testthat::expect(!all(is.na(effort$method_2)),
                   paste("Second method not used in 2-method property", i))

  properties_two <- assign_in(properties_two, list(i, "effort"), effort)
  setTxtProgressBar(pb, i)
}
close(pb)
str(properties_two)

# -----------------------------------------------------------------
# 3-method properties ----
# -----------------------------------------------------------------

n_three_method <- n_rel$n_simulate[3]
three_method_props <- get_properties(3)

# initiate list with the number of properties we need ----
# list of things we need to keep track of
property_attributes <- list(
  num = NULL,
  method_1 = NULL,
  method_2 = NULL,
  method_3 = NULL,
  area = NULL,
  effort = NULL
)

properties_three <- rep(list(property_attributes), n_three_method)
properties_three <- place_vec(properties_three, n_three_method, "num", 1:n_three_method)

c3 <- df |>
  filter(property %in% three_method_props) |>
  select(property, method) |>
  distinct() |>
  group_by(property) |>
  mutate(m = paste0("m_", 1:n())) |>
  ungroup() |>
  pivot_wider(names_from = m,
              values_from = method) |>
  select(starts_with("m_")) |>
  count(m_1, m_2, m_3) |>
  mutate(prob = n / sum(n))


draws <- sample(1:nrow(c3), n_three_method, prob = c3$prob, replace = TRUE)

### use index samples to select method combinations
sample_three_method <- c3 |>
  slice(draws) |>
  select(starts_with("m_"))

m1 <- sample_three_method |> pull(m_1)
m2 <- sample_three_method |> pull(m_2)
m3 <- sample_three_method |> pull(m_3)

## place the pair of methods in the 3-method data list ----
properties_three <- place_vec(properties_three, n_three_method, "method_1", m1)
properties_three <- place_vec(properties_three, n_three_method, "method_2", m2)
properties_three <- place_vec(properties_three, n_three_method, "method_3", m3)

## property area for 3-method properties ----
three_method_areas <- subset_area_n_method(three_method_props) |>
  mutate(id = 1:n())

## assign areas to 3-method properties ----
areas <- 1:n_three_method |>
  map(\(x) get_area_n(three_method_areas, c(m1[x], m2[x], m3[x]))) |>
  list_c()

properties_three <- place_vec(properties_three, n_three_method, "area", round(areas, 2))

## joint (2 and 3) and single return intervals ----
joint_return <- n_return_intervals(three_method_props, 3)

n_reps_three_method <- get_n_reps_joint(three_method_props)





