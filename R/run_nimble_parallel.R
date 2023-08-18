run_nimble_parallel <- function(cl, model_code, model_data, model_constants,
                                model_inits, n_iter, params_check, state.col,
                                model_flags, monitors_add = NULL, custom_samplers = NULL,
                                effective_size = 5000,
                                max_iter = NULL, max_psrf = 100,
                                calculate = TRUE, use_conjugacy = TRUE,
                                resetMV = FALSE,
                                save_iter = FALSE, dest = NULL){
  require(nimble)
  require(coda)
  source("R/functions_nimble.R")

  if(is.null(max_iter)) max_iter <- n_iter * 100

  export <- c(
    "model_code",
    "model_data",
    "model_constants",
    "n_iter",
    "calculate",
    "resetMV",
    "use_conjugacy",
    "monitors_add",
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
  start <- Sys.time()
  out <- clusterEvalQ(cl, {
    library(nimble)
    library(nimbleHMC)
    library(coda)
    # library(dplyr)

    list2env(model_flags, .GlobalEnv)
    likelihood_nb <- ifelse(likelihood == "nb", TRUE, FALSE)
    likelihood_binom <- ifelse(likelihood == "binomial", TRUE, FALSE)
    likelihood_poisson <- ifelse(likelihood == "poisson", TRUE, FALSE)

    setwd("C:/Users/John.Foster/OneDrive - USDA/Desktop/fosteR/pigs-statistical/")
    source("R/functions_nimble.R")

    Rmodel <- nimbleModel(
      code = model_code,
      constants = model_constants,
      data = model_data,
      inits = init,
      calculate = calculate
    )


    if(!calculate){
      calc <- Rmodel$calculate()
      if(is.infinite(calc) | is.nan(calc) | is.na(calc)){
        stop(paste0("Model log probability is ", calc))
      }
    }

    # default MCMC configuration
    mcmcConf <- configureMCMC(Rmodel,
                              useConjugacy = use_conjugacy)


    if(!is.null(custom_samplers)){
      for(i in seq_len(nrow(custom_samplers))){
        node <- custom_samplers$node[i]
        type <- custom_samplers$type[i]
        mcmcConf$removeSampler(node)
        mcmcConf$addSampler(node, type)
      }
    }
    # mcmcConf <- configureHMC(Rmodel)

    if(!is.null(monitors_add)){
      mcmcConf$addMonitors(monitors_add)
    }

    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc)
    Cmcmc$run(niter = 1000, nburnin = 0)
    samples <- as.matrix(Cmcmc$mvSamples)
    return(samples)
  })
  message("Model compiltion and initial ", 1000, " iterations complete:")
  print(Sys.time() - start)

  # convert samples to mcmc list to check for convergence
  mcmc <- as.mcmc.list(lapply(out, as.mcmc))

  # diagnostic <- continue_mcmc(params_check)
  # done <- diagnostic$done

  # if convergence is reached and we have enough samples, or something is wrong,
  # exit and return samples
  # if(done | is.na(done)){
  #   mcmc_split <- split_out(mcmc, state.col)
  #   return(list(
  #     params = mcmc_split$params,
  #     predict = mcmc_split$predict,
  #     diagnostic = diagnostic$checks
  #   ))
  # }

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

  # save_samples(save_iter, resetMV, state.col)

  # if the script has made it to this point we need to keep sampling
  continue <- TRUE
  c <- 0
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

    mcmc <- as.mcmc.list(lapply(out2, as.mcmc))

    # if(resetMV & c > 1){
    #   mcmc1 <- get_mcmc_chunks(dest)$params
    # } else {
    #   mcmc1 <- mcmc
    # }
    #
    total_iter <- c * n_iter
    message("Total run time for ", total_iter, " iterations:")
    print(Sys.time() - start)

    diagnostic <- continue_mcmc(mcmc, params_check, effective_size, max_psrf)

    if(c %% 10 == 0){
      message("Checking convergence on all iterations...")
      mm <- get_mcmc_chunks(path = dest, start = 1)
      check_mm <- continue_mcmc(mm$params, params_check, effective_size, max_psrf)
    }

    # save output if requested
    save_samples(save_iter, resetMV, state.col)

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
