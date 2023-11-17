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

one_prop_n <- n_rel$n_simulate[1]
one_method_prop <- get_properties(1)

# list of things we need to keep track of
property_attributes <- list(
  num = NULL,
  method = NULL,
  area = NULL,
  observations = NULL,
  n_reps = NULL
)

# initiate list with the number of properties we need ----
properties_one <- rep(list(property_attributes), one_prop_n)

## function to place vectors into each property attributes' list ----
place_vec <- function(ls, n, where, what){
  require(purrr)
  ls |> map2(1:n, \(x, y) assign_in(x, where, what[y]))
}

properties_one <- place_vec(properties_one, one_prop_n, "num", 1:one_prop_n)

## relative proportion of properties that use a specific method ----
## for sampling from to assigned 1-method properties a method
one_method_prob <- df |>
  filter(property %in% one_method_prop) |>
  select(property, method) |>
  distinct() |>
  count(method) |>
  mutate(prob = n / sum(n))

## assign methods to 1-method properties ----
sample_one_method <- sample(one_method_prob$method,
                            one_prop_n,
                            one_method_prob$prob,
                            replace = TRUE)

properties_one <- place_vec(properties_one, one_prop_n, "method", sample_one_method)

## function to sample property area from a given set of properties and methods ----
get_area <- function(dat, m){
  dat |>
    filter(method == m) |>
    pull(property.size) |>
    sample(1)
}

one_method_areas <- df |>
  filter(property %in% one_method_prop) |>
  select(property, method, property.size) |>
  distinct()

## assign areas to 1-method properties ----
areas <- sample_one_method |>
  map_vec(\(x) get_area(one_method_areas, x))

properties_one <- place_vec(properties_one, one_prop_n, "area", round(areas, 2))

## now we need to determine which PPs are sampled for each property ----
### start with determining the first PP ----

### function to assign a stating PP given method ----
start_pp <- function(m){
  if(grepl("Helicopter", m) | grepl("Fixed", m)){
    high <- 0.45 * n_pp
  } else {
    high <- 0.75 * n_pp
  }
  round(runif(1, 0.5, high))
}

start <- vapply(sample_one_method, start_pp, 1)

### lookup table for return intervals by method for single method-properties ----
one_method_return <- df |>
  filter(property %in% one_method_prop) |>
  select(property, timestep, method) |>
  distinct() |>
  group_by(property, method) |>
  mutate(delta = c(0, diff(timestep)),
         single_occurance = if_else(n() == 1, 1, 0),
         first_occurance = if_else(delta == 0 & single_occurance == 0, 1, 0))

m <- sample_one_method[1]
s <- start[1]

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
one_method_reps <- df |>
  filter(property %in% one_method_prop) |>
  group_by(property, timestep, method) |>
  count() |>
  ungroup()

#### function to sample number of removal events in a PP given method
get_reps <- function(reps_df, m, size){
  n_reps_single_method |>
    filter(method == m) |>
    pull(n) |>
    sample(size, replace = TRUE)
}

### sample the number of observations and reps, place in properties list
for(i in seq_len(one_prop_n)){
  sample_occasions <- get_sample_occasions(one_method_return, sample_one_method[i], start[i], n_pp)
  properties_one <- assign_in(properties_one, list(i, "observations"), sample_occasions)

  n_reps <- get_reps(one_method_reps, sample_one_method[i], length(sample_occasions))
  properties_one <- assign_in(properties_one, list(i, "n_reps"), n_reps)
}


