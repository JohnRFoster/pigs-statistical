# figure out if we need to keep sampling
continue_mcmc <- function(mcmc, nodes, effective_size, max_psrf){
  require(coda)
  # check convergence and effective sample size on specified node
  subset_check_mcmc <- function(node){
    s <- mcmc[,grep(node, colnames(mcmc[[1]]), value = TRUE, fixed = TRUE)]
    df <- data.frame(
      psrf = gelman.diag(s, multivariate = FALSE)$psrf[,2], # checking the upper CI
      effective_samples = effectiveSize(s)
    )
    if(nrow(df) == 1) rownames(df) <- node
    return(df)
  }

  message("Checking convergence and sample size")
  checks <- map_dfr(lapply(nodes, subset_check_mcmc), as.data.frame)
  converged <- all(checks$psrf < 1.1)
  enough_samples <- all(checks$effective_samples >= effective_size)
  funky <- any(is.nan(checks$psrf)) | max(checks$psrf) > max_psrf

  message('Convergence = ', converged)
  message('Enough effective samples = ', enough_samples)
  message('Mixing = ', !funky)

  done <- converged & enough_samples
  if(done) message("MCMC complete!")

  # TODO determine burnin and effective sample size POST burnin

  if(funky){
    message("Something is wrong with the mcmc!")
    done <- FALSE
  }

  if(!done | funky) print(checks)

  return(list(
    done = done,
    checks = checks
  ))
}

subset_check_burnin <- function(node, plot = FALSE){
  s <- mcmc[,grep(node, colnames(mcmc[[1]]), value = TRUE, fixed = TRUE)]
  if(plot){
    GBR <- gelman.plot(s)
  } else {
    ff <- tempfile()
    png(filename = ff)
    GBR <- gelman.plot(s)
    dev.off()
    unlink(ff)
  }
  shrink <- GBR$shrink[, , 2]
  if(all(shrink < 1.1)){
    burnin <- 1
  } else {
    if(is.null(dim(s$chain1))){
      burnin <- GBR$last.iter[tail(which(shrink > 1.1), 1) + 1]
    } else {
      burnin <- GBR$last.iter[tail(which(apply(shrink > 1.1, 1, any)), 1) + 1]
    }
  }
  return(burnin)
}

# need to combine previous sample chunks if resetMV = TRUE
get_mcmc_chunks <- function(path, start = 1, stop = NULL){
  out.files <- list.files(path)
  mcmc_files <- grep("_mcmcSamples.rds", out.files, value = TRUE)
  if(is.null(stop)) stop <- length(mcmc_files)
  chunk_params <- list()
  chunk_predict <- list()
  pb <- txtProgressBar(max = length(start:stop), style = 3)
  j <- 1
  for(i in start:stop){
    cN <- read_rds(file.path(path, mcmc_files[[i]]))
    chunk_params[[j]] <- cN$params
    chunk_predict[[j]] <- cN$predict
    j <- j + 1
    setTxtProgressBar(pb, i)
  }
  close(pb)
  mcmc_params <- runjags::combine.mcmc(chunk_params)
  mcmc_predict <- runjags::combine.mcmc(chunk_predict)
  return(list(params = mcmc_params, predict = mcmc_predict))
}

split_out <- function(mcmc_out, state.col){
  require(coda)
  mat2mcmc_list <- function(w) {
    temp <- list()
    chain.col <- which(colnames(w) == "CHAIN")
    for (i in unique(w[, "CHAIN"])) {
      temp[[i]] <- coda:::as.mcmc(w[w[, "CHAIN"] == i, -chain.col])
    }
    return(as.mcmc.list(temp))
  }

  out <- list(params = NULL, predict = NULL)
  mfit <- as.matrix(mcmc_out, chains = TRUE)
  pred.cols <- grep(paste0(state.col, "["), colnames(mfit), fixed = TRUE)
  chain.col <- which(colnames(mfit) == "CHAIN")
  out$predict <- mat2mcmc_list(mfit[, c(chain.col, pred.cols)])
  out$params <- mat2mcmc_list(mfit[, -pred.cols])
  return(out)
}

tau_2_sigma <- function(tau){
  1/sqrt(tau)
}

sigma_2_tau <- function(sigma){
  1/sigma^2
}

sigma_2_var <- function(sigma){
  sigma^2
}

var_2_sigma <- function(var){
  sqrt(var)
}

var_2_tau <- function(var){
  1/var
}

tau_2_var <- function(tau){
  1/tau
}

sum_properties <- nimbleFunction(

  run = function(
    property = double(1),  # property vector index
    z = double(1)          # latent abundance vector
  ){
    N <- sum(z[property])
    return(N)
    returnType(double(0))
  },
  buildDerivs = TRUE
)

calc_log_potential_area <- nimbleFunction(
  run = function(
      log_rho = double(1),
      log_gamma = double(1),
      p_unique = double(1),
      log_effort_per = double(0),
      effort_per = double(0),
      n_trap_m1 = double(0),
      log_pi = double(0),
      method = double(0)
  ){
    m <- method

    if(m == 1){ # firearms
      log_potential_area <- log_rho[m] +
        log_effort_per -
        log(1 + (p_unique[m] * n_trap_m1))
    } else if(m == 2 | m == 3){ # fixed wing and helicopter
      log_potential_area <- log_rho[m] + log_effort_per
    } else if(m == 4 | m == 5){
      log_potential_area <- log_pi +
        (2 * (log_rho[m] + log_effort_per -
                log(exp(log_gamma[m-3]) + effort_per))) +
        log(1 + (p_unique[m-2] * n_trap_m1))
    }
    return(log_potential_area)
    returnType(double(0))
  },
  buildDerivs = TRUE
)


make_inits_function <- function(inits_dir, constants = NULL, data = NULL){
  if(is.null(inits_dir)){
    inits <- function(){
      N_init <- (rowSums(data$y, dims = 2, na.rm = TRUE) + 1)
      ls <- list(
        N = N_init,
        x = log(N_init),
        log_gamma = rnorm(5),
        log_rho = rnorm(5),
        beta_p = matrix(rnorm(constants$n_method*constants$n_beta_p, 0, 0.1), constants$n_method, constants$n_beta_p),
        p_mu = rnorm(constants$n_method),
        log_r_mu = runif(1, 0, 1),
        tau_proc = runif(1, 1, 10)
      )
      return(ls)
    }
  } else {
    samples <- read_rds(file.path(inits_dir, "thinnedSamples.rds"))
    params_mu <- apply(samples$params, 2, mean)
    states_mu <- apply(samples$predict, 2, mean)

    beta_p <- params_mu[grep("beta_p", names(params_mu))]
    log_gamma <- params_mu[grep("log_gamma", names(params_mu))]
    log_rho <- params_mu[grep("log_rho", names(params_mu))]
    log_r_mu <- params_mu[grep("log_r_mu", names(params_mu))]
    p_unique <- params_mu[grep("p_unique", names(params_mu))]
    tau_proc <- params_mu[grep("tau_proc", names(params_mu))]

    nim <- read_rds(file.path(inits_dir, "nimbleList.rds"))
    unit_lookup <- nim$unit_lookup

    N <- round(exp(states_mu)+1) |>
      as_tibble() |>
      mutate(property = unit_lookup$property_idx,
             timestep = unit_lookup$timestep) |>
      pivot_wider(names_from = timestep,
                  values_from = value) |>
      select(-property) |>
      as.matrix()

    inits <- function(){
      list(
        N = N,
        x = log(N + 1),
        log_gamma = rnorm(length(log_gamma), log_gamma, 0.05),
        log_rho = rnorm(length(log_rho), log_gamma, 0.05),
        p_mu = rnorm(length(p_unique), boot::logit(p_unique), 0.05),
        beta_p = matrix(rnorm(length(beta_p), beta_p, 0.05), 5, 4),
        log_r_mu = rnorm(length(log_r_mu), log_r_mu, 0.05),
        tau_proc = abs(rnorm(1, tau_proc, 0.05)),
        beta_r = rnorm(1, 0, 0.05)
      )
    }

  }
  return(inits)
}


make_inits_function_dm <- function(inits_dir, constants, data, demographic_stochasticity){

  # samples <- read_rds("out/test/dm/dm_static/thinnedSamples.rds")
  # nim <- read_rds("out/test/dm/dm_static/nimbleList.rds")
  # unit_lookup <- nim$unit_lookup
  # states_mu <- apply(samples$predict, 2, mean)
  # N_init <- round(states_mu) |>
  #   as_tibble() |>
  #   mutate(property = unit_lookup$property_idx,
  #          timestep = unit_lookup$timestep) |>
  #   pivot_wider(names_from = timestep,
  #               values_from = value) |>
  #   select(-property) |>
  #   as.matrix()

  if(is.null(inits_dir)){

    logit_mean_phi <- runif(1, 1.5, 2)
    sigma_phi <- runif(1, 0.5, 0.75)

    mean_lpy <- 1
    mean_ls <- round(runif(1, 3.5, 6.4))
    zeta <- rep(mean_lpy / 365 * constants$pp_len * mean_ls, constants$n_pp)

    p_mu <- jitter(nimble::logit(c(0.6, 0.4, 0.5, 0.3, 0.25)))
    dm <- matrix(NA, constants$n_property, max(constants$all_pp, na.rm = TRUE))
    phi_track <- 0
    S <- R <- dm
    for(i in 1:constants$n_property){
      dm[i, constants$all_pp[i, 1]] <- constants$N_init[i] + rpois(1, 200)
      for(t in 2:constants$n_pp_prop[i]){
        phi <- nimble::ilogit(rnorm(1, logit_mean_phi, sigma_phi))
        phi_track <- c(phi_track, phi)
        n_avail <- max(0, dm[i, constants$all_pp[i, t-1]] - data$rem[i, t-1])
        if(is.na(n_avail)) print(t)
        if(demographic_stochasticity){
          S[i, t-1] <- rbinom(1, n_avail, phi)
          R[i, t-1] <- rpois(1, zeta[constants$all_pp[i, t-1]] * n_avail/2)
        } else {
          S[i, t-1] <- round(n_avail * phi)
          R[i, t-1] <- round(zeta[constants$all_pp[i, t-1]] * n_avail/2)
        }
        dm[i, constants$all_pp[i, t]] <- S[i, t-1] + R[i, t-1]
      }
    }

    N <- matrix(NA, constants$n_property, constants$n_pp)
    for(i in 1:constants$n_property){
      for(t in 1:constants$n_timesteps[i]){ # loop through sampled PP only
        N[i, constants$PPNum[i, t]] <- dm[i, constants$PPNum[i, t]]
      }
    }
    # N[N == 0] <- 1
    inits <- function(){ list(
      N = N,
      lambda_1 = (apply(N, 1, function(x) x[min(which(!is.na(x)))])),
      log_gamma = log(runif(2, 0.5, 2)),
      log_rho = log(c(runif(1, 0.1, 3),
                      runif(1, 8, 15),
                      runif(1, 8, 15),
                      runif(1, 0.75, 1.5),
                      runif(1, 0.75, 1.5))),
      beta_p = matrix(15, 5, 3),
      beta1 = rnorm(5),
      p_mu = p_mu,
      logit_mean_phi = logit_mean_phi,
      tau_phi = sigma_2_tau(sigma_phi),
      # log_mean_lpy = log(mean_lpy),
      log_mean_ls = log(mean_ls),
      S = S,
      R = R
    )}
    return(inits)
  }
}
