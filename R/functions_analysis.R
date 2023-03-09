
######
# check burnin on a subset of nodes
######

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



######
# load output
######
get_samples <- function(path, ...){
  # source("R/functions_nimble.R")
  out_mcmc <- get_mcmc_chunks(path, ...)
  model_inputs <- read_rds(file.path(path, "nimbleList.rds"))
  return(
    list(
      params = out_mcmc$params,
      predict = out_mcmc$predict,
      model_inputs = model_inputs
    )
  )
}

######
# check convergence and burnin
######
convergence_check <- function(mcmc, params_check, verbose = TRUE){
  psrf <- coda::gelman.diag(mcmc, multivariate = FALSE)$psrf
  converged <- max(psrf[, "Upper C.I."]) < 1.1
  message("Convergence = ", converged)
  if(converged){
    burnin <- subset_check_burnin(params_check)
  } else {
    if(verbose) {
      message("Non-converged nodes:")
      psrf |>
        as_tibble() |>
        mutate(node = rownames(psrf)) |>
        select(node, `Point est.`, `Upper C.I.`) |>
        filter(`Upper C.I.` >= 1.1) |>
        as.data.frame() |>
        print()
    }
    message("Setting burnin start to iteration 1000")
    burnin <- 1000
  }
  return(
    list(
      psrf = psrf,
      burnin = burnin
    )
  )
}

######
# get a random sample of posteriors to work with later
######
thin_samples <- function(path, start, override_trace = FALSE, ...){

  out_mcmc <- get_samples(path, start = start)

  params <- out_mcmc$params
  predict <- out_mcmc$predict
  params_check <- out_mcmc$model_inputs$params_check

  traceplots <- any(grepl("traceplots", list.files(path)))
  plot <- if_else(traceplots, override_trace, TRUE)

  if(plot){
    message("creating traceplots...")
    pdf(
      file = file.path(path, "traceplots.pdf"),
      onefile = TRUE
    )
    plot(params)
    dev.off()
    message("  done")
  }

  converge <- convergence_check(params, params_check)
  burnin <- converge$burnin

  param_burn <- window(params, start = burnin)
  params_mat <- as.matrix(param_burn)

  predict_burn <- window(predict, start = burnin)
  predict_mat <- as.matrix(predict_burn)

  draws <- sample.int(nrow(params_mat), 5000, replace = TRUE)

  params_samps <- params_mat[draws,]
  predict_samps <- predict_mat[draws,]

  write_rds(list(params = params_samps, predict = predict_samps),
            file.path(path, "thinnedSamples.rds"))
}


######
# combine and rename parameters
######
make_parameters_tb <- function(mcmc_mat, spatial, model_name, quants = TRUE, county_lookup = NULL, dest){
  out <- mcmc_mat |>
    as_tibble() |>
    pivot_longer(cols = everything(),
                 names_to = "node")

  mutate_method <- function(df){
    df |>
      mutate(
        method = if_else(grepl("[1", node, fixed = TRUE), "Firearms", "none"),
        method = if_else(grepl("[2", node, fixed = TRUE), "Fixed wing", method),
        method = if_else(grepl("[3", node, fixed = TRUE), "Helicopter", method),
        method = if_else(grepl("[4", node, fixed = TRUE), "Snare", method),
        method = if_else(grepl("[5", node, fixed = TRUE), "Trap", method)
      )
  }


  obs_betas <- out |>
    filter(grepl("beta_p", node)) |>
    mutate(
      name = "Observation covariate",
      parameter = if_else(grepl(", 1]", node, fixed = TRUE), "Intercept", node),
      parameter = if_else(grepl(", 2]", node, fixed = TRUE), "Slope (road density)", parameter),
      parameter = if_else(grepl(", 3]", node, fixed = TRUE), "Slope (ruggedness)", parameter),
      parameter = if_else(grepl(", 4]", node, fixed = TRUE), "Slope (canopy cover)", parameter),
      value = if_else(parameter == "Intercept", boot::inv.logit(value), value)
    ) |>
    mutate_method()

  gamma <- out |>
    filter(grepl("log_gamma", node)) |>
    mutate(
      name = "Saturation constant",
      parameter = "Saturation constant",
      value = exp(value)
    ) |>
    mutate_method()

  rho <- out |>
    filter(grepl("log_rho", node)) |>
    mutate(
      name = "Max area",
      parameter = "Max area",
      value = exp(value)
    ) |>
    mutate_method()

  p <- out |>
    filter(grepl("p_unique", node)) |>
    mutate(
      name = "Proportion unique area",
      parameter = "Proportion unique area"
    ) |>
    mutate_method()

  lambda <- out |>
    filter(node == "log_r_mu") |>
    mutate(
      name = "Population growth",
      parameter = "Population growth",
      value = exp(value))

  tau_proc <- out |>
    filter(node == "tau_proc") |>
    mutate(
      name = "Process error",
      parameter = "Process error",
      value = 1/sqrt(value))

  beta_r <- tibble()
  if("beta_r" %in% unique(out$node)){
    beta_r <- out |>
      filter(node == "beta_r") |>
      mutate(name = "Density dependence effect",
             parameter = "Density dependence effect")
  }

  all_params <- bind_rows(
    obs_betas,
    gamma,
    rho,
    p,
    lambda,
    tau_proc,
    beta_r
  )

  if(spatial){
    tau_car <- out |>
      filter(node == "tau_car") |>
      mutate(
        name = "CAR SD",
        parameter = "CAR SD",
        value = 1/sqrt(value))

    s <- out |>
      filter(grepl("s[", node, fixed = TRUE)) |>
      mutate(
        name = "County spatial autocorrelation",
        parameter = "County spatial autocorrelation",
        cnty_idx = as.numeric(str_extract(node, "\\d*(?=\\])")),
        value = log(value))

    all_params <- bind_rows(
      all_params,
      tau_car,
      s
    )
  }

  if(quants){
    if("cnty_idx" %in% colnames(all_params)){
      params_group <- all_params |>
        left_join(county_lookup) |>
        group_by(node, name, parameter, method, cnty_idx, cnty_name, state, island)
    } else {
      params_group <- all_params |>
        group_by(node, name, parameter, method)
    }

    params_stats <- params_group |>
      summarise(lower95 = quantile(value, 0.025),
                lower75 = quantile(value, 0.125),
                median = quantile(value, 0.5),
                upper75 = quantile(value, 0.875),
                upper95 = quantile(value, 0.975)) |>
      ungroup()
  } else {
    params_stats <- all_params
  }

  xx <- params_stats |>
    mutate(model = model_name) |>
    ungroup()

  write_rds(xx, dest)
  return(xx)

}

######
# combine and score predictive posteriors
######
rmse <- function(actual, predicted){
  sqrt(mean((actual - predicted)^2))
}

pred_post <- function(constants, data, pp_N, pp_yp, model_name, dest){
  data$N <- pp_N
  data$yp <- pp_yp
  data$model_name <- model_name
  ls <- purrr::prepend(constants, data)
  pp <- with(ls, {
    abundance <- pred_obs <- tibble()
    for(i in 1:n_property){
      for(t in 1:n_timesteps[i]){
        value <- N[i, t, ]
        if(any(is.na(value))) message(i, " ", t)
        qs_a <- tibble(
          lower95 = quantile(value, 0.025),
          lower75 = quantile(value, 0.125),
          median = quantile(value, 0.5),
          upper75 = quantile(value, 0.875),
          upper95 = quantile(value, 0.975),
          observed = sum(y[i, t, 1:n_reps[i, t]]),
          rmse = rmse(observed, value),
          crps = scoringRules::crps_sample(observed, value),
          property_idx = i,
          timestep_idx = t,
          model = model_name
        )
        abundance <- bind_rows(abundance, qs_a)

        for(j in 1:n_reps[i, t]){
          value <- yp[i, t, j, ]
          if(any(is.na(value))) message(i, " ", t, " ", j)
          qs_y <- tibble(
            lower95 = quantile(value, 0.025),
            lower75 = quantile(value, 0.125),
            median = quantile(value, 0.5),
            upper75 = quantile(value, 0.875),
            upper95 = quantile(value, 0.975),
            observed = data$y[i, t, j],
            rmse = rmse(observed, value),
            crps = scoringRules::crps_sample(observed, value),
            property_idx = i,
            timestep_idx = t,
            pass = j,
            model = model_name
          )
          pred_obs <- bind_rows(pred_obs, qs_y)
        }
      }
    }
    list(abundance = abundance, pred_obs = pred_obs)
  })
  write_rds(pp, dest)
  return(pp)
}



#######################################################################
# plot quantiles
#######################################################################

# plot line ranges
gg_stats <- function(gg, s = 2, w75 = 1, w95 = 0.5){
  gg +
    geom_point(aes(y = median), size = s) +
    geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = w95) +
    geom_linerange(aes(ymin = lower75, ymax = upper75), linewidth = w75)
}

######
# population constant population growth
######
gg_lambda_constant <- function(df, ...){
  gg <- df |>
    filter(parameter == "Population growth") |>
    ggplot() +
    aes(x = model, color = model)

  gg_stats(gg, ...) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme_bw()

}

######
# observation model intercepts
######
gg_obs_intercepts <- function(df, ...){
  df |>
    filter(parameter == "Intercept") |>
    ggplot() +
    aes(x = method, color = model) +
    geom_point(aes(y = median), size = 2, position = position_dodge(width=0.5)) +
    geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = 0.5, position = position_dodge(width=0.5)) +
    geom_linerange(aes(ymin = lower75, ymax = upper75), linewidth = 1, position = position_dodge(width=0.5)) +
    theme_bw() +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "Base capture rate (random intercepts by method)",
         y = "Probability of capture | individual in the surveyed area",
         x = "Method")

}

######
# observation model slopes
######
gg_obs_slopes <- function(df, ...){
  df |>
    filter(grepl("Slope", parameter)) |>
    ggplot() +
    aes(x = parameter, color = model) +
    geom_point(aes(y = median), size = 2, position = position_dodge(width=0.5)) +
    geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = 0.5, position = position_dodge(width=0.5)) +
    geom_linerange(aes(ymin = lower75, ymax = upper75), linewidth = 1, position = position_dodge(width=0.5)) +
    theme_bw() +
    facet_wrap(~ method) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(title = "Effect of land cover on observation probability",
         y = "Effect",
         x = "Land cover covariate")
}

######
# observation model slopes
######
gg_gamma <- function(df, ...){
  df |>
    filter(grepl("gamma", node)) |>
    ggplot() +
    aes(x = parameter, color = model) +
    geom_point(aes(y = median), size = 2, position = position_dodge(width=0.5)) +
    geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = 0.5, position = position_dodge(width=0.5)) +
    geom_linerange(aes(ymin = lower75, ymax = upper75), linewidth = 1, position = position_dodge(width=0.5)) +
    theme_bw() +
    facet_wrap(~ method, scales = "free") +
    labs(title = "Saturation constant",
         y = "",
         x = "")
}
######
# observation model slopes
######
gg_rho <- function(df, ...){
  df |>
    filter(grepl("rho", node)) |>
    ggplot() +
    aes(x = parameter, color = model) +
    geom_point(aes(y = median), size = 2, position = position_dodge(width=0.5)) +
    geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = 0.5, position = position_dodge(width=0.5)) +
    geom_linerange(aes(ymin = lower75, ymax = upper75), linewidth = 1, position = position_dodge(width=0.5)) +
    theme_bw() +
    facet_wrap(~ method, scales = "free") +
    labs(title = "Maximum area effected (km^2)",
         y = "Maximum area effected (km^2)",
         x = "")
}


gg_proc_error <- function(df, ...){
  df |>
    filter(parameter == "Process error") |>
    ggplot() +
    aes(x = model, color = model) +
    geom_point(aes(y = median), size = 2, position = position_dodge(width=0.5)) +
    geom_linerange(aes(ymin = lower95, ymax = upper95), linewidth = 0.5, position = position_dodge(width=0.5)) +
    geom_linerange(aes(ymin = lower75, ymax = upper75), linewidth = 1, position = position_dodge(width=0.5)) +
    theme_bw() +
    labs(title = "Process error",
         y = "",
         x = "")

}










