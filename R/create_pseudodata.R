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
for(i in seq_len(n_one_method)){
  sample_occasions <- get_sample_occasions(one_method_return, sample_one_method[i], start[i], n_pp)
  n_reps <- get_reps(n_reps_single_method, sample_one_method[i], length(sample_occasions))

  effort <- tibble(sample_occasions = sample_occasions, n_reps = n_reps)

  properties_one <- assign_in(properties_one, list(i, "effort"), effort)
}


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
combos <- df |>
  filter(property %in% two_method_props) |>
  select(property, method) |>
  distinct() |>
  group_by(property) |>
  mutate(n = 1) |>
  ungroup() |>
  pivot_wider(names_from = method,
              values_from = n,
              values_fill = 0) |>
  select(all_of(m_vec)) |>
  as.matrix() |>
  crossprod()

combos[lower.tri(combos, diag = TRUE)] <- NA

method_pairs <- combos |>
  as_tibble() |>
  mutate(method = m_vec) |>
  pivot_longer(cols = -method,
               names_to = "method2",
               values_to = "count",
               values_drop_na = TRUE) |>
  filter(count > 0)

### the probability of using a given method first for 2-method properties ----
two_method_prob <- method_pairs |>
  mutate(prob = count / sum(count))

### sample the index of rows in probability table given probability column
pairs_sample <- sample.int(nrow(method_pairs), n_two_method, prob = two_method_prob$prob, replace = TRUE)

### use index samples to select method pairs
sample_two_method <- two_method_prob |>
  slice(pairs_sample) |>
  select(starts_with("method"))

m1 <- sample_two_method |> pull(method)
m2 <- sample_two_method |> pull(method2)

## place the pair of methods in the 2-method data list ----
properties_two <- place_vec(properties_two, n_two_method, "method_1", m1)
properties_two <- place_vec(properties_two, n_two_method, "method_2", m2)

## function to sample property area from a given set of properties and methods ----
two_method_areas <- df |>
  filter(property %in% two_method_props) |>
  select(property, method, property.size) |>
  distinct() |>
  group_by(property) |>
  mutate(n = paste0("n_", 1:n())) |>
  ungroup() |>
  pivot_wider(names_from = n,
              values_from = method)

get_area_2 <- function(m1, m2){
  two_method_areas |>
    filter((n_1 == m1 & n_2 == m2) | (n_1 == m2 & n_2 == m1)) |>
    pull(property.size) |>
    sample(1)
}

## assign areas to 2-method properties ----
areas <- map2_vec(m1, m2, get_area_2)

properties_two <- place_vec(properties_two, n_two_method, "area", round(areas, 2))

## joint and single return intervals ----
joint_return <- df |>
  filter(property %in% two_method_props) |>
  select(property, timestep, method) |>
  distinct() |>
  pivot_wider(names_from = method,
              values_from = method) |>
  unite(method, m_vec, sep = "&", na.rm = TRUE) |>
  group_by(property, method) |>
  mutate(delta = c(0, diff(timestep)),
         single_occurance = if_else(n() == 1, 1, 0),
         first_occurance = if_else(delta == 0 & single_occurance == 0, 1, 0)) |>
  ungroup() |>
  separate(method, c("method_1", "method_2"), sep = "&") |>
  mutate(joint = if_else(!is.na(method_2), 1, 0))

## function to create a return interval df given two methods ----
get_joint_return <- function(m1, m2){
  pairs_return <- joint_return |>
    filter(joint == 1,
           ((method_1 == m1 & method_2 == m2) |
              (method_1 == m2 & method_2 == m1)))

  single_return <- joint_return |>
    filter(joint == 0,
           method_1 == m1 | method_1 == m2)

  bind_rows(single_return, pairs_return) |>
    filter(delta > 0)
}


## function to generate sample occasions of single or joint use of methods
get_sample_occasions_two <- function(return_df, m1, m2, max_pp){
  start <- min(vapply(c(m1, m2), start_pp, 1))

  obs <- return_df |>
    slice(sample.int(nrow(return_df), 1))
  effort <- obs |>
    select(starts_with("method")) |>
    mutate(sample_occasions = start)

  # end_pp <- start

  create_effort_df <- function(start_effort, start_pp){
    effort <- start_effort
    end_pp <- start_pp
    for(xx in start_pp:max_pp){
      obs <- return_df |>
        slice(sample.int(nrow(return_df), 1))

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

      if(end_pp > max_pp) break

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

n_reps_two_method <- df |>
  filter(property %in% two_method_props) |>
  group_by(property, timestep, method) |>
  count() |>
  ungroup() |>
  pivot_wider(names_from = method,
              values_from = n)

#### function to sample number of removal events in a PP given a pair of methods ----
get_reps_two <- function(m1, m2){

  m <- c(m1, m2)
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
for(i in seq_len(n_two_method)){
  sample_occasions <- get_joint_return(m1[i], m2[i]) |>
    get_sample_occasions_two(m1[i], m2[i], n_pp)

  m1_i <- sample_occasions |> pull(method_1)
  m2_i <- sample_occasions |> pull(method_2)

  n_reps <- map2(m1_i, m2_i, get_reps_two) |> list_rbind()
  effort <- bind_cols(sample_occasions, n_reps)

  testthat::expect(effort$sample_occasions[nrow(effort)] <= n_pp,
                   paste("Last PP exceeds boundary for property", i))
  testthat::expect(all(!is.na(effort$method_2)),
                   paste("Second method not used in property", i))

  properties_two <- assign_in(properties_two, list(i, "effort"), effort)
  # glimpse(effort)
}

str(properties_two)

