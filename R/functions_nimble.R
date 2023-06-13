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





run_nimble_parallel <- function(cl, model_code, model_data, model_constants,
                                model_inits, n_iter, params_check, state.col,
                                model_flags, custom_samplers = NULL,
                                effective_size = 5000, n_burnin = NULL,
                                max_iter = NULL, max_psrf = 100,
                                calculate = TRUE, use_conjugacy = TRUE,
                                resetMV = FALSE,
                                save_iter = FALSE, dest = NULL){
  library(nimble)
  library(coda)

  if(is.null(n_burnin)) n_burnin <- round(n_iter / 4)
  if(is.null(max_iter)) max_iter <- n_iter * 100

  export <- c(
    "model_code",
    "model_data",
    "model_constants",
    "n_iter",
    "n_burnin",
    "calculate",
    "resetMV",
    "use_conjugacy",
    "params_check",
    "custom_samplers",
    "model_flags"
  )

  clusterExport(cl, export, envir = environment())

  for(i in seq_along(cl)){
    set.seed(i)
    init <- model_inits()
    clusterExport(cl[i], "init", envir = environment())
  }

  message("Compiling model and initial parallel sampling...")
  c <- 1
  start <- Sys.time()
  out <- clusterEvalQ(cl, {
    library(nimble)
    library(coda)
    # library(dplyr)

    list2env(model_flags, .GlobalEnv)
    likelihood_nb <- dplyr::if_else(likelihood == "nb", TRUE, FALSE)
    likelihood_binom <- dplyr::if_else(likelihood_nb, FALSE, TRUE)

    exponential <- dplyr::if_else(process_type == "exponential", TRUE, FALSE)
    ricker <- dplyr::if_else(process_type == "ricker", TRUE, FALSE)
    gompertz <- dplyr::if_else(process_type == "gompertz", TRUE, FALSE)
    jamiesonBrooks <- dplyr::if_else(process_type == "jamiesonBrooks", TRUE, FALSE)
    dennisTaper <- dplyr::if_else(process_type == "dennisTaper", TRUE, FALSE)

    zeta_constant <- dplyr::if_else(process_type == "static", TRUE, FALSE)
    zeta_pp <- dplyr::if_else(process_type == "zeta_pp", TRUE, FALSE)

    zeta_pp <- dplyr::if_else(grepl("zeta_pp", process_type), TRUE, FALSE)
    zeta_constant <- !zeta_pp

    phi_pp <- dplyr::if_else(grepl("phi_pp", process_type), TRUE, FALSE)
    phi_constant <- !phi_pp

    Rmodel <- nimbleModel(
      code = model_code,
      constants = model_constants,
      data = model_data,
      inits = init,
      calculate = calculate
    )

    Cmodel <- compileNimble(Rmodel)

    if(!calculate){
      calc <- Cmodel$calculate()
      if(is.infinite(calc) | is.nan(calc) | is.na(calc)){
        stop(paste0("Model log probability is ", calc))
      }
    }

    # default MCMC configuration
    mcmcConf <- configureMCMC(Rmodel,
                              useConjugacy = use_conjugacy,
                              monitors = c("xn", params_check))

    if(!is.null(custom_samplers)){
      for(i in seq_len(nrow(custom_samplers))){
        node <- custom_samplers$node[i]
        type <- custom_samplers$type[i]
        mcmcConf$removeSampler(node)
        mcmcConf$addSampler(node, type)
      }
    }

    Rmcmc <- buildMCMC(mcmcConf)
    Cmcmc <- compileNimble(Rmcmc)
    Cmcmc$run(niter = n_iter, nburnin = n_burnin)
    samples <- as.matrix(Cmcmc$mvSamples)
    return(samples)
  })
  message(n_iter, " iterations complete:")
  print(Sys.time() - start)


  # convert samples to mcmc list to check for convergence
  mcmc <- as.mcmc.list(lapply(out, as.mcmc))

  # figure out if we need to keep sampling
  continue_mcmc <- function(nodes){

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
      done <- NA
    }

    if(!done | funky) print(checks)

    return(list(
      done = done,
      checks = checks
    ))
  }

  diagnostic <- continue_mcmc(params_check)
  done <- diagnostic$done

  # if convergence is reached and we have enough samples, or something is wrong,
  # exit and return samples
  if(done | is.na(done)){
    mcmc_split <- split_out(mcmc, state.col)
    return(list(
      params = mcmc_split$params,
      predict = mcmc_split$predict,
      diagnostic = diagnostic$checks
    ))
  }

  save_samples <- function(save_iter, resetMV, state.col){

    mcmc_split <- split_out(mcmc, state.col)

    if(save_iter){
      out_rds <- list(
        params = mcmc_split$params,
        predict = mcmc_split$predict,
        diagnostic = diagnostic$checks
      )
      f <- "mcmcSamples.rds"
      if(resetMV){
        if(c < 10) c <- paste0("0", c)
        f <- paste0(c, "_", f)
      }
      message("  Writing samples")
      if(!dir.exists(dest)){
        dir.create(dest, recursive = TRUE, showWarnings = FALSE)
      }
      write_rds(out_rds, file.path(dest, f))
    }
  }

  save_samples(save_iter, resetMV, state.col)

  # if the script has made it to this point we need to keep sampling
  continue <- TRUE

  while(continue){

    # run mcmc
    c <- c + 1
    message("Parallel sampling ", c)
    start2 <- Sys.time()
    out2 <- clusterEvalQ(cl, {
      Cmcmc$run(n_iter, reset = FALSE, resetMV = resetMV)
      return(as.mcmc(as.matrix(Cmcmc$mvSamples)))
    })
    message(n_iter, " iterations complete:")
    print(Sys.time() - start2)

    mcmc_c <- as.mcmc.list(lapply(out2, as.mcmc))
    mcmc <- mcmc_c

    # save output if requested
    save_samples(save_iter, resetMV, state.col)

    if(resetMV) mcmc1 <- get_mcmc_chunks(dest)$params

    total_iter <- nrow(mcmc1[[1]])
    message("Total run time for ", total_iter, " iterations:")
    print(Sys.time() - start)

    diagnostic <- continue_mcmc(params_check)
    # save_samples(save_iter, resetMV)

    done <- diagnostic$done
    # if(done | is.na(done)){
    #   mcmc_split <- split_out(mcmc, state.col)
    #   return(list(
    #     params = mcmc_split$params,
    #     predict = mcmc_split$predict,
    #     diagnostic = diagnostic$checks
    #   ))
    # }

    continue_samp <- total_iter < max_iter
    continue <- if_else(continue_samp, !done, FALSE)

  }
  message("Maximum iterations reached. Returning current samples.")
  return(list(
    samples = as.mcmc.list(lapply(out2, as.mcmc)),
    diagnostic = diagnostic$checks
  ))
}

sum_properties <- nimbleFunction(
  run = function(
    property = double(1),  # property vector index
    z = double(1)          # latent property abundance vector
  ){
    N <- sum(z[property])
    return(N)
    returnType(double(0))
  }
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


make_inits_function_dm <- function(inits_dir, process_type, constants, data,
                                   phi_constant, phi_pp, zeta_constant, zeta_pp){

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

    if(phi_constant){
      phi <- rep(runif(1, 0.85, 0.9), constants$n_pp)
    } else if(phi_pp){
      logit_phi_global <- logit(runif(1, 0.9, 0.95))
      tau_phi <- jitter(3)
      logit_phi <- rnorm(constants$n_pp, logit_phi_global, 1/sqrt(tau_phi))
      phi <- ilogit(logit_phi)
    }

    if(zeta_constant){
      zeta <- rep(runif(1, 0.1, 0.11), constants$n_pp)
    } else if(zeta_pp){
      log_zeta_global <- log(runif(1, 0.1, 0.11))
      tau_zeta <- jitter(3)
      log_zeta <- rnorm(constants$n_pp, log_zeta_global, 1/sqrt(tau_zeta))
      zeta <- exp(log_zeta)
    }

    p_mu <- rnorm(constants$n_method, -1, 1)
    N_init <- as.matrix(round((rowSums(data$y_wide, dims = 2, na.rm = TRUE)+400) / 0.2))
    N <- N_init
    dm <- matrix(NA, constants$n_property, max(constants$all_pp, na.rm = TRUE))
    S <- R <- dm
    for(i in 1:constants$n_property){
      if((constants$n_timesteps[i]+1) < ncol(N_init)){
        N_init[i, (constants$n_timesteps[i]+1):ncol(N_init)] <- 0
      }

      dm[i, constants$all_pp[i, 1]] <- N_init[i, 1]
      for(t in 2:constants$n_pp_prop[i]){
        S[i, t] <- rbinom(1, dm[i, constants$all_pp[i, t-1]], phi[constants$all_pp[i, t-1]])
        R[i, t] <- rpois(1, zeta[constants$all_pp[i, t-1]] * dm[i, constants$all_pp[i, t-1]])
        dm[i, constants$all_pp[i, t]] <- S[i, t] + R[i, t]
      }
    }

    inits <- function(){ list(
      N = N_init,
      z1 = rpois(constants$n_property, N_init[,1]),
      size = runif(constants$n_county, 0.01, 5),
      log_gamma = rnorm(5),
      log_rho = rnorm(5),
      beta_p = matrix(rnorm(constants$n_method*constants$n_beta_p, 0, 0.1), constants$n_method, constants$n_beta_p),
      p_mu = p_mu,
      logit_phi = logit(phi[1]),
      log_zeta = log(zeta[1]),
      # log_zeta_global = log_zeta_global,
      # tau_zeta = tau_zeta,
      # logit_phi_global = logit_phi_global,
      # tau_phi = tau_phi,
      S = S,
      R = R
    )}
    return(inits)
  }
}

