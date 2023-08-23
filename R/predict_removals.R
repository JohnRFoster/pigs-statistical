
n_fx <- 3
n_reps <- 3

scenarios <- data.frame(t1 = c("N", "R"))
if(n_fx > 1){
  for(i in 2:n_fx){
    scenarios <- cbind(scenarios, c("N", "R"))
  }
  colnames(scenarios) <- paste0("t", 1:n_fx)
  scenarios <- expand.grid(scenarios)
}
scenarios <- t(scenarios)

rep_num <- 1 # starting with 13 there are 50 properties
out_dir <- "out/simulation"
model_dir <- "DM_recruitData_varyingEffort"
likelihood <- "poisson"
mcmc_config <- "customMCMC_conjugate"
rep <- paste0("simulation_", rep_num)

dest <- file.path(out_dir, model_dir, likelihood, mcmc_config, paste0("simulation_", 1))
rds <- read_rds(file.path(dest, "posterior.rds"))
samples <- rds$samples
sim_data <- read_rds(file.path(dest, "sim_data.rds"))
attach(sim_data$constants)
attach(sim_data$data)

sample_occasions <- tibble(
  property = property_x,
  PPNum = pp_x
) |>
  group_by(property, PPNum) |>
  mutate(timestep = cur_group_id()) |>
  ungroup() |>
  mutate(n_id = 1:n())

property_df <- tibble(
  y = sim_data$data$y,
  property = p_property_idx,
  sampled_pp = p_pp_idx,
  start = start,
  end = end,
  method = method,
  county = p_county_idx
) |>
  group_by(property, sampled_pp) |>
  mutate(timestep = cur_group_id()) |>
  ungroup() |>
  mutate(n_id = 1:n())

take <- sim_data$take


n_samples <- samples[, grep("xn[", colnames(samples), fixed = TRUE)]
n_mcmc <- nrow(n_samples)

pattern <- "(?<!beta_)p\\[\\d*\\]"
p_detect <- str_detect(colnames(samples), pattern)
p <- samples[,which(p_detect)]

mu_phi <- samples[,"logit_mean_phi"]
sigma_phi <- samples[,"sigma_phi"]
mean_ls <- samples[,"mean_ls"]
zeta <- mean_ls*28*1/365

log_rho <- samples[, grep("log_rho", colnames(samples))]
log_gamma <- samples[, grep("log_gamma", colnames(samples))]
p_unique <- ilogit(samples[, grep("p_mu", colnames(samples))])
beta_p <- samples[, grep("beta_p", colnames(samples))]

X <- unique(X_p)

# storage
N_predict <- tibble()
Y_predict <- tibble()

get_e <- function(e){
  lapply(m, function(i){
    E |> filter(method == i)
  }) |>
    bind_rows() |>
    pull(e)
}

for(s in 1:nrow(scenarios)){
  scenario_name <- as.vector(scenarios[,s])
  for(i in 1:n_property){

    prop <- property_df |>
      filter(property == i,
             timestep == 1)

    observations <- sample_occasions |>
      filter(property == i)

    start_pp <- min(observations$PPNum)
    end_pp <- max(observations$PPNum)

    n_id <- min(observations$n_id)
    n_node <- paste0("xn[", n_id, "]")
    N <- n_samples[,grep(n_node, colnames(n_samples), fixed = TRUE)]

    take_property <- take |>
      filter(property == i)

    E <- take_property |>
      group_by(method, area_property) |>
      summarise(effort_per = mean(effort_per),
                trap_count = mean(trap_count)) |>
      ungroup()

    c <- take_property$county[1]

    for(t in start_pp:end_pp){
      for(fx in 1:n_fx){


        # do we remove pigs?
        if(scenario_name[fx] == "R"){

          # removal methods based on frequency of use in a given property
          m <- sample(take_property$method, n_reps, replace = TRUE)

          # extract effort attributes from property
          ep <- get_e("effort_per")
          tc_m1 <- get_e("trap_count") - 1

          # the search area given method and effort
          log_potential_area <- calc_log_potential_area(
            method = m,
            log_rho = log_rho,
            log_gamma = log_gamma,
            p_unique = p_unique,
            effort_per = ep,
            n_trap_m1 = tc_m1
          )

          # probability of capture for each pass
          p <- calc_p(
            X = X[c,],
            beta = beta_p,
            log_potential_area = log_potential_area,
            area_property = E$area_property[1]
          )

          # remove
          y <- matrix(NA, n_mcmc, n_reps)
          for(j in seq_len(n_reps)){
            y[,j] <- rpois(n_mcmc, N * p[,j])
            y_rep <- tibble(
              y = y[,j],
              property = i,
              PPNum = t,
              method = m[j],
              effort_per = ep,
              trap_count = tc_m1,
              rep = j
            )
            Y_predict <- bind_rows(Y_predict, y_rep)
          }

          N <- N - rowSums(y)
        }

        phi <- ilogit(rnorm(n_mcmc, mu_phi, sigma_phi))
        S <- rbinom(n_mcmc, N, phi)
        R <- rpois(n_mcmc, N * zeta)
        N <- S + R

        # storage
        N_predict <- bind_rows(N_predict,
                               tibble(
                                 property = i,
                                 start_pp = t, # the PP we started from (before removals, if applicable)
                                 PPNum = t+fx, # the PP we are forecasting into
                                 horizon = fx, # the number of PPs being forecasted across
                                 removal = scenario_name[fx],
                                 S = S,
                                 R = R,
                                 N = N
                               ))
      } # fx







    }
  }
}

samps <- rds$samples[, grep("xn[", colnames(rds$samples), fixed = TRUE)]
N <- left_join(sim_data$county_to_property, sim_data$N)$abundance

# samps <- rds$y
# N <- sim_data$take$take

N1 <- apply(samps, 2, quantile, 0.025)
N2 <- apply(samps, 2, quantile, 0.975)
M <- apply(samps, 2, quantile, 0.5)
tibble(
  N1 = N1,
  N2 = N2,
  Predicted = M,
  Known = N
) |>
  ggplot() +
  aes(x = Known, y = Predicted, ymin = N1, ymax = N2) +
  # geom_linerange() +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(intercept = 0, slope = 1, color = "red") +
  theme_bw()
