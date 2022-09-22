install.packages("rstan", repos="http://cran.us.r-project.org")
install.packages("V8", repos="http://cran.us.r-project.org")
library(rstan)
rstan_options(auto_write = TRUE) #R/fit-model20210413.R
#options(mc.cores = parallel::detectCores())
options(mc.cores=4)
setwd('/project/iwctml/mtabak/APHIS/abundance/wild-pigs/')
#stan_d <- readRDS('/project/iwctml/mtabak/APHIS/abundance/wild-pigs/stan_d.rds')
stan_d <- readRDS('./stan_d.rds')


pars <- c('beta', 
          'beta_p', 
          'gamma', 
          'log_rho',
          'p_unique',
          'log_pr_area_sampled', 
          'log_pr_area_sampled_dev', 
          'log_p', 
          'log_p_dev',
          'log_lambda', 
          'sigma_property',
          'z_short',
          'sigma_short',
          'loglik_train',
          'loglik_dev',
          'y_rep',
          'y_rep_dev',
          'eps_property',
          'sigma_property_p', 
          'eps_property_pR', 
          'eps_stR', 
          'eta', 
          'sigma_st0', 
          'sigma_st'
          )

#m_init <- stan_model('/project/iwctml/mtabak/APHIS/abundance/wild-pigs/models/saturating.stan')
m_init <- stan_model('./models/saturating.stan')

m_fit <- vb(m_init, data = stan_d, pars = pars,
            init = 0, tol_rel_obj = 0.007,
            adapt_engaged = FALSE, eta = 0.2)

#saveRDS(m_fit, file = "/project/iwctml/mtabak/APHIS/abundance/wild-pigs/fitted_mods/saturating_fit.rds")
saveRDS(m_fit, file = "./fitted_mods/saturating_fit_VA_20210413.rds")

