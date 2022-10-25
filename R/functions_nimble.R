library(nimble)
library(coda)



run_nimble_parallel <- function(cl, model_code, model_data, model_constants,
                                model_inits, n_iter, params_check,
                                effective_size = 5000, n_burnin = NULL,
                                max_iter = NULL, max_psrf = 1000, calculate = TRUE,
                                use_conjugacy = FALSE){

  if(is.null(n_burnin)) n_burnin <- round(n_iter / 4)
  if(is.null(max_iter)) max_iter <- n_iter * 100

  export <- c(
    "model_code",
    "model_data",
    "model_constants",
    "n_iter",
    "n_burnin",
    "calculate",
    "use_conjugacy"
  )

  clusterExport(cl, export, envir = environment())

  for(i in seq_along(cl)){
    set.seed(i)
    init <- model_inits()
    clusterExport(cl[i], "init", envir = environment())
  }

  message("Initial parallel sampling...")
  c <- 1
  start <- Sys.time()
  out <- clusterEvalQ(cl, {
    library(nimble)
    library(coda)

    Rmodel <- nimbleModel(
      code = model_code,
      constants = model_constants,
      data = model_data,
      inits = init,
      calculate = calculate
    )

    # default MCMC configuration
    mcmcConf <- configureMCMC(Rmodel,
                              useConjugacy = use_conjugacy)

    # these print statements will not display when running in parallel
    mcmcConf$printMonitors()
    mcmcConf$printSamplers(byType = TRUE)

    mcmcConf$addMonitors("z_shortR")

    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc)
    Cmcmc$run(niter = n_iter, nburnin = n_burnin)
    samples <- as.mcmc(as.matrix(Cmcmc$mvSamples))
    return(samples)
  })
  message(n_iter, " iterations complete:")
  print(Sys.time() - start)


  # convert samples to mcmc list to check for convergence
  mcmc <- as.mcmc.list(out)

  # figure out if we need to keep sampling
  # TODO write code to save mcmc intermitently
  continue_mcmc <- function(nodes){

    # check convergence and effective sample size on specified node
    subset_check_mcmc <- function(node){
      s <- mcmc[,grep(node, colnames(mcmc[[1]]), value = TRUE, fixed = TRUE)]
      data.frame(
        psrf = gelman.diag(s, multivariate = FALSE)$psrf[,2], # checking the upper CI
        effective_samples = effectiveSize(s)
      )
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
    return(list(
      samples = mcmc,
      diagnostic = diagnostic$checks
    ))
  }

  # if the script has made it to this point we need to keep sampling
  continue <- TRUE

  while(continue){
    c <- c + 1
    message("Parallel sampling ", c)
    start2 <- Sys.time()
    out2 <- clusterEvalQ(cl, {
      Cmcmc$run(n_iter, reset = FALSE)
      return(as.mcmc(as.matrix(Cmcmc$mvSamples)))
    })
    message(n_iter, " iterations complete:")
    print(Sys.time() - start2)

    mcmc <- as.mcmc.list(lapply(out2, as.mcmc))
    total_iter <- nrow(mcmc[[1]])
    message("Total run time for ", total_iter, " iterations:")
    print(Sys.time() - start)

    diagnostic <- continue_mcmc(params_check)
    done <- diagnostic$done
    if(done | is.na(done)){
      return(list(
        samples = mcmc,
        diagnostic = diagnostic$checks
      ))
    }

    continue_samp <- total_iter < max_iter
    continue <- if_else(continue_samp, !done, FALSE)

  }
  message("Maximum iterations reached. Returning current samples.")
  return(mcmc)
}





