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

    list2env(model_flags, .GlobalEnv)
    exponential <- dplyr::if_else(process_type == "exponential", TRUE, FALSE)
    ricker <- dplyr::if_else(process_type == "ricker", TRUE, FALSE)
    gompertz <- dplyr::if_else(process_type == "gompertz", TRUE, FALSE)
    jamiesonBrooks <- dplyr::if_else(process_type == "jamiesonBrooks", TRUE, FALSE)
    dennisTaper <- dplyr::if_else(process_type == "dennisTaper", TRUE, FALSE)

    Rmodel <- nimbleModel(
      code = model_code,
      constants = model_constants,
      data = model_data,
      inits = init,
      calculate = calculate
    )

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
    Cmodel <- compileNimble(Rmodel)
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

# function from https://groups.google.com/g/nimble-users/c/i8U54pSGSLs/m/QlCS8hjEBQAJ
count_by_site <- nimbleFunction(
  run = function(
    group = double(1),
    z = double(1),
    ngroup = integer(0)
  ) {

    M <- length(group)
    N <- numeric(ngroup)

    for(i in 1:M) {
      if(z[i] == 1) {
        N[group[i]] <- N[group[i]] + 1
      }
    }

    return(N) # brackets are only necessary in model code, not nimbleFunction
    returnType(double(1))
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


