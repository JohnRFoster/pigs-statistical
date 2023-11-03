inits <- function(data, constants, dir = NULL){
  with(append(data, constants), {

    if(is.null(dir)){
      beta1 <- rnorm(n_method, -1, 0.25)
      beta_p <- matrix(rnorm(m_p*n_method, 0, 0.1), n_method, m_p)
      p_mu <- rnorm(3)
      log_gamma <- log(runif(2, 0.1, 2))
      log_rho <- log(
        c(runif(1, 0.1, 3), runif(1, 8, 15), runif(1, 20, 25), runif(1, 0.75, 1.5), runif(1, 0.75, 1.5))
      )
      psi_phi <- runif(1, 2, 4)
      phi_mu <- runif(1, 0.7, 0.8)
      mean_ls <- round(runif(1, 5, 8))
    } else {
      post <- read_rds(file.path(dir, "posterior.rds"))
      init_mu <- apply(post$samples, 2, mean)
      beta_p <- matrix(jitter(init_mu[grep("beta_p[", names(init_mu), fixed = TRUE)]), constants$n_method, constants$m_p)
      beta1 <- jitter(init_mu[grep("beta1[", names(init_mu), fixed = TRUE)])
      p_mu <- jitter(init_mu[grep("p_mu[", names(init_mu), fixed = TRUE)])
      phi_mu <- jitter(init_mu[grep("phi_mu", names(init_mu), fixed = TRUE)])
      psi_phi <- jitter(init_mu[grep("psi_phi", names(init_mu), fixed = TRUE)])
      log_mean_ls <- jitter(init_mu[grep("log_mean_ls", names(init_mu), fixed = TRUE)])
      mean_ls <- exp(log_mean_ls)
      log_gamma <- jitter(init_mu[grep("log_gamma[", names(init_mu), fixed = TRUE)])
      log_rho <- jitter(init_mu[grep("log_rho[", names(init_mu), fixed = TRUE)])
    }

    n_init <- numeric(n_property)
    for(i in 1:n_property){

      p_survey <- which(p_property_idx == i)
      r <- which(min(p_pp_idx[p_survey]) == p_pp_idx[p_survey])

      nr <- log_theta <- numeric(length(r))
      for(j in 1:length(r)){
        idx <- r[j]

        log_potential_area <- calc_log_potential_area(
          log_rho = log_rho,
          log_gamma = log_gamma,
          p_unique = ilogit(p_mu),
          log_effort_per = log_effort_per[idx],
          effort_per = effort_per[idx],
          n_trap_m1 = n_trap_m1[idx],
          log_pi = log_pi,
          method = method[idx]
        )

        log_theta[j] <- log(ilogit(beta1[method[idx]] + inprod(X_p[idx, ], beta_p[method[idx],]))) +
          min(0, log_potential_area - log_survey_area_km2[idx])

        if(j == 1){
          p <- exp(log_theta[1])
        } else {
          p <- exp(log_theta[1] + sum(log(1 - exp(log_theta[1:(j-1)]))))
        }
        swine <- max(y[idx], 3)
        # print(swine)
        nr[j] <- swine / p

      }
      n_init[i] <- rpois(1, mean(nr))
    }

    # n_init <- rpois(n_property, property_area*8)

    a <- phi_mu * psi_phi
    b <- (1 - phi_mu) * psi_phi
    mean_lpy <- 1
    zeta <- mean_lpy / 365 * pp_len * mean_ls
    phi <- numeric(max(pH, na.rm = TRUE))
    dm <- matrix(NA, n_property, max(all_pp, na.rm = TRUE))
    N <- S <- R <- Z <- dm
    for(i in 1:n_property){
      N[i, all_pp[i, 1]] <- n_init[i]
      for(t in 2:n_pp_prop[i]){
        Z[i, t-1] <- N[i, all_pp[i, t-1]]
        S[i, t-1] <- Z[i, t-1] * rbeta(1, a, b)
        R[i, t-1] <- zeta * Z[i, t-1] / 2
        N[i, all_pp[i, t]] <- round(S[i, t-1] + R[i, t-1])
      }
    }

    buffer <- 500
    list(
      log_lambda_1 = log(n_init + buffer),
      beta_p = beta_p,
      beta1 = beta1,
      p_mu = p_mu,
      phi_mu = phi_mu,
      psi_phi = psi_phi,
      N = N + buffer,
      log_mean_ls = log(mean_ls),
      log_gamma = log_gamma,
      log_rho = log_rho
    )
  })

}




