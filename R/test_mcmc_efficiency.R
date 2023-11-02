library(compareMCMCs)
init_vals <- inits()

nimbleMCMCdefs <- list(

  slice_scalars = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c("tau_phi", "beta1", "log_mean_ls", "logit_mean_phi"))
    mcmcConf$addSampler(target = c("tau_phi"), type = "slice")
    mcmcConf$addSampler(target = c("beta1"), type = "slice")
    mcmcConf$addSampler(target = c("log_mean_ls"), type = "slice")
    mcmcConf$addSampler(target = c("logit_mean_phi"), type = "slice")
    mcmcConf
  },

  slice_scalars1 = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c("tau_phi", "beta1", "log_mean_ls", "logit_mean_phi"))
    mcmcConf$addSampler(target = c("tau_phi"), type = "slice")
    mcmcConf$addSampler(target = c("beta1"), type = "slice")
    mcmcConf$addSampler(target = c("log_mean_ls"), type = "slice")
    mcmcConf$addSampler(target = c("logit_mean_phi"), type = "ess")
    mcmcConf
  },

  slice_scalars2 = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c("tau_phi", "beta1", "log_mean_ls"))
    mcmcConf$addSampler(target = c("tau_phi"), type = "slice")
    mcmcConf$addSampler(target = c("beta1"), type = "slice")
    mcmcConf$addSampler(target = c("log_mean_ls"), type = "slice")
    mcmcConf
  },

  slice_scalars3 = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c("tau_phi", "beta1", "log_mean_ls", "beta_p"))
    mcmcConf$addSampler(target = c("tau_phi"), type = "slice")
    mcmcConf$addSampler(target = c("beta1"), type = "slice")
    mcmcConf$addSampler(target = c("log_mean_ls"), type = "slice")
    mcmcConf$addSampler(target = c("beta_p"), type = "RW_block")
    mcmcConf
  },

  slice_scalars4 = function(model){
    mcmcConf <- configureMCMC(model)
    mcmcConf$removeSamplers(c("tau_phi", "beta1", "log_mean_ls", "logit_mean_phi", "beta_p"))
    mcmcConf$addSampler(target = c("tau_phi"), type = "slice")
    mcmcConf$addSampler(target = c("beta1"), type = "slice")
    mcmcConf$addSampler(target = c("log_mean_ls"), type = "slice")
    mcmcConf$addSampler(target = c("logit_mean_phi"), type = "slice")
    mcmcConf$addSampler(target = c("beta_p"), type = "RW_block")
    mcmcConf
  }

)


mcmcResults_nimble <- compareMCMCs(
  modelInfo = list(code = modelCode,
                   data = data,
                   constants = constants, # centered
                   inits = init_vals),
  MCMCs = c('slice_scalars',
            'slice_scalars1',
            'slice_scalars2',
            'slice_scalars3',
            'slice_scalars4'
  ),
  nimbleMCMCdefs = nimbleMCMCdefs,
  monitors = c("beta_p", "beta1", "log_gamma", "log_rho", "log_mean_ls",
               "p_mu", "tau_phi", "logit_mean_phi"),
  MCMCcontrol = list(niter = 5000, burnin = 100)
)

# MCMC efficiency for a parameter is defined as the effective sample size divided by computation time in seconds
# It is the number of effectively independent samples generated per second.
res <- combineMetrics(mcmcResults_nimble)
res$byMCMC
least_efficient <- res$byParameter |>
  as_tibble() |>
  group_by(MCMC) |>
  filter(efficiency == min(efficiency))
least_efficient

res$byParameter |>
  as_tibble() |>
  filter(MCMC == "RWblock_BR_slice_tau_beta1_ls") |>
  arrange(efficiency)
