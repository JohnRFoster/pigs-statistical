fit_simulations_parallel <- function(cl,
                                     rep_num,
                                     phi_mu,
                                     psi_phi,
                                     sigma_dem,
                                     n_pp,
                                     modelCode,
                                     n_iter, n_chains, params_check,
                                     monitors_add = NULL, custom_samplers = NULL,
                                     inits_dir = NULL, calculate = TRUE,
                                     use_conjugacy = TRUE, dest = NULL){
  require(nimble)
  require(coda)
  require(tidyverse)
  require(lubridate)
  source("R/functions_nimble.R")

  export <- c(
    "rep_num",
    'phi_mu',
    'psi_phi',
    "sigma_dem",
    "n_pp",
    "modelCode",
    "n_iter",
    "n_chains",
    "calculate",
    "use_conjugacy",
    "monitors_add",
    "params_check",
    "custom_samplers",
    "inits_dir",
    "dest"
  )

  clusterExport(cl, export, envir = environment())

  out <- clusterEvalQ(cl, {
    library(nimble)
    library(coda)
    library(tidyverse)
    library(targets)

    task_id <- Sys.getpid()
    rep <- paste0("simulation_", rep_num, "_", task_id)

    dest <- file.path(dest, rep)
    if(!dir.exists(dest)) dir.create(dest, recursive = TRUE, showWarnings = FALSE)

    a_phi <- phi_mu * psi_phi
    b_phi <- (1 - phi_mu) * psi_phi

    property_attributes <- expand_grid(
      property_density = round(runif(5, 0.3, 6), 2),
      proportion_county_surveyed = seq(0.1, 0.9, length.out = 5)
    ) |>
      mutate(
        county = as.numeric(as.factor(proportion_county_surveyed)),
        phi_mu = phi_mu,
        psi_phi = psi_phi,
        area_property = round(runif(n(), 5, 250), 2),
        initial_abundnace = round(property_density * area_property),
        n_sample_occasions = round(runif(n(), 2.5, 12.4))
      ) |>
      group_by(county) |>
      mutate(
        area_county_surveyed = sum(area_property),
        area_county = area_county_surveyed / proportion_county_surveyed,
        area_not_sampled = area_county - area_county_surveyed
      ) |>
      ungroup() |>
      arrange(county) |>
      mutate(property = 1:n())

    n_county <- max(property_attributes$county)
    c_road_den <- rnorm(n_county)
    c_rugged <- rnorm(n_county)
    c_canopy <- rnorm(n_county)
    beta_p <- matrix(rnorm(20), 5, 4)

    source("R/functions_nimble.R")
    source("R/functions_simulation.R")

    sim_data <- simulate_dm(
      property_data = property_attributes,
      likelihood = "poisson",
      phi_mu = phi_mu,
      psi_phi = psi_phi,
      sigma_dem = sigma_dem,
      n_pp = n_pp,
      c_road_den = c_road_den,
      c_rugged = c_rugged,
      c_canopy = c_canopy,
      beta_p = beta_p,
      demographic_stochasticity = TRUE,
      plot = TRUE
    )

    source("R/priors.R")
    phi_prior <- create_surv_prior()

    constants <- sim_data$constants
    constants$phi_mu_a <- phi_prior$alpha
    constants$phi_mu_b <- phi_prior$beta
    data <- sim_data$data

    spatial <- FALSE

    sim_data$params_check <- params_check
    sim_data$custom_samplers <- custom_samplers
    sim_data$phi_mu <- phi_mu
    sim_data$psi_phi <- psi_phi
    sim_data$sigma_dem <- sigma_dem

    write_rds(sim_data,
              file = file.path(dest, "sim_data.rds"))
    write_rds(property_attributes,
              file = file.path(dest, "property_attributes.rds"))

    source("R/dm_inits.R")

    Rmodel <- nimbleModel(
      code = modelCode,
      constants = constants,
      data = data,
      inits = inits(data, constants, inits_dir),
      calculate = TRUE
    )

    for(i in 1:constants$n_survey){
      N_model <- Rmodel$N[constants$p_property_idx[i], constants$p_pp_idx[i]]
      n <- N_model - data$y_sum[i]
      if(n < 0){
        Rmodel$N[constants$p_property_idx[i], constants$p_pp_idx[i]] <- N_model + n^2
      }
    }

    Rmodel$initializeInfo()

    # default MCMC configuration
    mcmcConf <- configureMCMC(Rmodel, useConjugacy = TRUE)
    mcmcConf$addMonitors(monitors_add)

    if(!is.null(custom_samplers)){
      for(i in seq_len(nrow(custom_samplers))){
        node <- custom_samplers$node[i]
        type <- custom_samplers$type[i]
        mcmcConf$removeSampler(node)
        mcmcConf$addSampler(node, type)
      }
    }

    for(i in 1:5){
      node <- paste0("beta_p[", i, ", ", 1:constants$m_p, "]")
      node <- c(paste0("beta1[", i, "]"), node)
      mcmcConf$removeSampler(node)
      mcmcConf$addSampler(node, "AF_slice")
    }

    mcmcConf$printSamplers(byType = TRUE)

    Rmcmc <- buildMCMC(mcmcConf)
    Cmodel <- compileNimble(Rmodel)
    Cmcmc <- compileNimble(Rmcmc)

    # n_iter <- 5
    # n_chains <- 1

    samples <- runMCMC(
      Cmcmc,
      niter = n_iter,
      nchains = n_chains,
      nburnin = n_iter / 2,
      samplesAsCodaMCMC = TRUE
    )

    all_nodes <- colnames(samples[[1]])

    j <- unlist(lapply(params_check, function(x) grep(x, all_nodes)))

    params <- samples[,j]

    message("Calculating PSRF...")
    psrf <- gelman.diag(params, multivariate = FALSE)
    print(psrf)

    message("Calculating effective samples...")
    effective_samples <- effectiveSize(params)
    print(effective_samples)

    message("Calculating burnin...")
    # ff <- tempfile()
    # png(filename = ff)
    GBR <- gelman.plot(params)
    # dev.off()
    # unlink(ff)

    burnin <- GBR$last.iter[tail(which(apply(GBR$shrink[, , 2] > 1.1, 1, any)), 1) + 1]
    print(burnin)

    if(is.na(burnin) | burnin >= 0.9*nrow(params[[1]])){
      burnin <- round(n_iter / 4)
    }

    samples_burn_mcmc <- window(samples, start = burnin)
    params_burn <- samples_burn_mcmc[, j]

    message("Creating traceplots...")
    png(filename = file.path(dest, "mcmcTimeseries%03d.png"))
    plot(params_burn)
    dev.off()
    message("  done")

    samples_burn_mat <- as.matrix(samples_burn_mcmc)
    draws <- sample.int(nrow(samples_burn_mat), 1000, replace = TRUE)

    samples_draw <- as.data.frame(samples_burn_mat[draws, ])


    source("R/functions_predict.R")
    post <- data_posteriors(samples_draw, constants, data)

    write_csv(samples_draw, file = file.path(dest, "posteriorSamples.csv"))
    write_csv(as.data.frame(post$y), file = file.path(dest, "posteriorPredictedTake.csv"))
    write_csv(as.data.frame(post$potential_area), file = file.path(dest, "posteriorPredictedPotentialArea.csv"))

    write_rds(
      list(psrf = psrf$psrf,
           effective_samples = effective_samples,
           burnin = burnin,
           converged = max(psrf$psrf[,"Upper C.I."]) <= 1.1,
           bad_mcmc = any(is.na(psrf$psrf)),
           task_id = paste0(rep_num, "_", task_id)),
      file = file.path(dest, "posteriorEval.rds")
    )

  })
}

