# need to combine previous sample chunks if resetMV = TRUE
get_mcmc_chunks <- function(path, start = 1, stop = NULL){
  out.files <- list.files(path)
  mcmc_files <- grep("_mcmcSamples.rds", out.files, value = TRUE)
  if(is.null(stop)) stop <- length(mcmc_files)
  chunk_params <- list()
  chunk_predict <- list()
  pb <- txtProgressBar(max = length(start:stop), style = 3)
  for(i in start:stop){
    cN <- read_rds(file.path(path, mcmc_files[[i]]))
    chunk_params[[i]] <- cN$params
    chunk_predict[[i]] <- cN$predict
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
      if(is.null(dim(s$chain1))) rownames(df) <- node
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


# make_inits_function <- function(inits_dir, constants = NULL, data = NULL){
#   if(is.null(inits_dir)){
#     inits <- function(){
#       N_init <- (rowSums(data$y, dims = 2, na.rm = TRUE) + 1) * 5
#       ls <- list(
#         N =  N_init,
#         x =  log(N_init),
#         log_gamma = rnorm(1),
#         log_rho = rnorm(1),
#         beta_p = matrix(rnorm(constants$n_method*constants$n_beta_p, 0, 0.1), constants$n_method, constants$n_beta_p),
#         p_unique = rbeta(constants$n_method, 1, 1),
#         r_mu = rexp(1, 1),
#         tau_proc = runif(1, 1, 10)
#       )
#       return(ls)
#     }
#   } else {
#     samples <- read_rds(file.path(inits_dir, "thinnedSamples.rds"))
#     params_mu <- apply(samples$params, 2, mean)
#     states_mu <- apply(samples$predict, 2, mean)
#
#     inits
#
#   }
#   return(inits)
# }


