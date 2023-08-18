
# n_fx <- 3
#
# expand_grid(
#   R = "R",
#   N = "N"
# )
#
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

for(m in 1:n_mcmc){
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

    m_freq <- table(take_property$method)/length(take_property$method)

    E <- take_property |>
      group_by(method, area_property) |>
      summarise(effort_per = mean(effort_per),
                trap_count = mean(trap_count))

    c <- take_property$county[1]

    for(t in start_pp:end_pp){

      # need frequencies of each sample method for each data set (or property?? or county??)
      # can use these frequencies to simulate removal events for each MCMC sample
      # for population abundance through time all we need is the total removed

      phi <- ilogit(rnorm(n_mcmc, mu_phi, sigma_phi))
      S <- rbinom(n_mcmc, N, phi)
      R <- rpois(n_mcmc, N * zeta)
      N <- S + R

      # storage
      N_predict <- bind_rows(N_predict,
                             tibble(
                               property = i,
                               PPNum = t,
                               S = S,
                               R = R,
                               N = N
                             ))

      # removals
      m <- sample.int(5, n_reps, prob = m_freq)

      log_potential_area <- matrix(NA, n_mcmc, n_reps)
      for(j in seq_len(n_reps)){
        if(m[j] == 4 | m[j] == 5){
          log_potential_area[,j] <- calc_log_potential_area_trap_snare(
            log_rho = log_rho[,m[j]],
            log_gamma = log_gamma[,m[j]-3],
            p_unique = p_unique[,m[j]],
            effort_per = E$effort_per[m[j]],
            n_trap_m1 = E$trap_count[m[j]] - 1
          )
        } else {
          log_potential_area[,j] <- calc_log_potential_area_shooting(
            log_rho = log_rho[,m[j]],
            p_unique = p_unique[,m[j]],
            effort_per = E$effort_per[m[j]],
            n_trap_m1 = E$trap_count[m[j]] - 1
          )
        }
      }

      p <- matrix(NA, n_mcmc, n_reps)
      for(mc in seq_len(n_mcmc)){
        log_theta <- numeric(n_reps)
        for(j in seq_len(n_reps)){

          # base probability of capture given an individual is the first survey
          log_theta[j] <- log(ilogit(inprod(X_p[c,], beta_p[mc,]))) +
            pmin(0, log_potential_area[mc, j] - log(E$area_property[m[j]]))

          # the probability an individual is captured on the first survey
          p[mc, j] <- exp(log_theta[1])

          # the probability an individual is captured after the first survey
          for(k in 2:n_reps){
            p[mc, k] <- exp(log_theta[1] +
                              sum(log(1 - exp(log_theta[1:(k-1)]))))
          }
        }
      }

      y <- matrix(NA, n_mcmc, n_reps)
      for(j in seq_len(n_reps)){
        y[,j] <- rpois(n_mcmc, N * p[,j])
        Y_predict <- bind_rows(
          Y_predict,
          tibble(
            y = y[,j],
            property = i,
            PPNum = t,
            method = m[j],
            effort_per = E$effort_per[m[j]],
            trap_count = E$trap_count[m[j]],
            rep = j
          )
        )
      }

      N <- N - rowSums(y)



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
