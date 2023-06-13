## functions to simulate pig population abundance

n_county <- 3
n_property <- 5 # number of properties in each county
n_pp <- 10
n_passes <- 5

c_road_den <- rnorm(n_county)
c_rugged <- rnorm(n_county)
c_canopy <- rnorm(n_county)

# ecological process
mean_r <- 1.1 # mean growth rate
sigma_proc <- 2.2 # process error

log_pi <- log(pi)

N <- array(NA, dim = c(n_county, n_property, n_pp)) # latent abunance [property, primary period] w/ stochasticity
x <- array(NA, dim = c(n_county, n_property, n_pp)) # expected log abundnace



for(i in seq_len(n_county)){
  for(j in seq_len(n_property)){

    # initial abundnace for each property
    x[i, 1, 1] <- log(max(round(rnorm(n_property, 50, 25)), 3))
    N[i, 1, 1] <- rpois(1, exp(x[i, 1, 1]))

    for(k in seq_len(n_pp-1)){

      mu <- x[i, j, k] + log(mean_r)              # process model
      x[i, j, k+1] <- rnorm(1, mu, sd_proc)       # add process error
      N[i, j, k+1] <- rpois(1, exp(x[i, j, k+1])) # add stochasticity

      for(l in seq_len(n_passes)){

        if(method[i, j, k, l] %in% 1:3){ # shooting data model
          log_potential_area <- log_pi +
            (2 * (log_rho[method[i, j, k, l]] +
                    log_effort_per[i, j, k, l] -
                    log(exp(log_gamma[method[i, j, k, l]]) +
                          effort_per[i, j, k, l]))) +
            log(1 + (p_unique[method[i, j, k, l]] * n_trap_m1[i, j, k, l]))
        } else { # trap/snare data model
          log_potential_area <- log_rho[method[i, j, k, l]] +
            log_effort_per[i, j, k, l] -
            log(exp(log_gamma[method[i, j, k, l]]) + effort_per[i, j, k, l]) +
            log(1 + (p_unique[method[i, j, k, l]] * n_trap_m1[i, j, k, l]))
        }

        if(l == 1){
          C[i, j, k, l] <- rbinom(1, N[i, j, k], p)
        } else {
          z <- N[i, j, k] - sum(C[i, j, k, 1:l])
          C[i, j, k, l] <- rbinom(1, z, p)
        }

      }
    }
  }
}
