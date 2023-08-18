data {
  // spatiotemporal data
  int<lower = 1> n_st;                       // number of spatiotemporal units
  int<lower = 1> m_n;                        // number of abundance parameters
  vector[n_st] area_km2;                     // property area for each st unit
  int n_property;
  int<lower = 1, upper = n_property> property[n_st];
  int<lower = 1> n_timestep;
  int<lower = 1, upper = n_timestep> timestep[n_st];

  // sparse matrix for spatiotemporal process params (fixed effects)
  matrix[n_st, m_n] X;
  matrix[n_st, m_n] X_use; // matrix telling stan to ignore NAs
  int<lower = 1> m_short;
  matrix[n_st, m_short] X_short;

  // data for spatial indexing, neighbors, ICAR priors
  int<lower = 1> n_county;
  int<lower = 0> n_island;
  int<upper = n_county> island_idx[n_island];
  int<lower = 1, upper = n_county> county_idx[n_st];
  int<lower = 0> n_edges;
  int<lower = 1, upper = n_county> node1[n_edges];
  int<lower = 1, upper = n_county> node2[n_edges];

  // training data: take of wild pigs
  int<lower = 1> n_survey;                   // number of surveys
  int<lower = 1, upper = n_st> survey_idx[n_survey];
  int<lower = 0> y[n_survey];                // number of pigs caught
  vector[n_survey] scaled_effort;
  int m_p;                              // num. cols in detection design matrix
  matrix[n_survey, m_p] X_p;
  int<lower = 1> order[n_survey];            // survey order (1st, 2nd, ...)
  vector<lower = 0>[n_survey] survey_area_km2;          // property area on each survey
  int<lower = 1, upper = n_county> p_county_idx[n_survey];
  int<lower = 1, upper = n_property> p_property_idx[n_survey];
  int<lower = 0, upper = n_survey> start[n_survey]; // start index: past surveys
  int<lower = 0, upper = n_survey> end[n_survey];   // end index: past surveys
  int<lower = 1> n_method;
  int<lower = 1, upper = n_method> method[n_survey];
  vector<lower = 0>[n_survey] effort;
  vector<lower = 0>[n_survey] effort_per;
  vector<lower = 1>[n_survey] trap_count;

  // dev data: take of wild pigs
  int<lower = 1> n_survey_dev;
  int<lower = 1, upper = n_st> survey_idx_dev[n_survey_dev];
  int<lower = 0> y_dev[n_survey_dev];
  vector[n_survey_dev] scaled_effort_dev;
  matrix[n_survey_dev, m_p] X_p_dev;
  int<lower = 1> order_dev[n_survey_dev];
  vector[n_survey_dev] survey_area_km2_dev;
  int<lower = 1, upper = n_county> p_county_idx_dev[n_survey_dev];
  int<lower = 1, upper = n_property> p_property_idx_dev[n_survey_dev];
  int<lower = 0, upper = n_survey_dev> start_dev[n_survey_dev];
  int<lower = 0, upper = n_survey_dev> end_dev[n_survey_dev];
  int<lower = 1, upper = n_method> method_dev[n_survey_dev];
  vector<lower = 0>[n_survey_dev] effort_dev;
  vector<lower = 0>[n_survey_dev] effort_per_dev;
  vector<lower = 1>[n_survey_dev] trap_count_dev;
}

transformed data {
  int first_survey[n_survey];
  int first_survey_dev[n_survey_dev];
  vector<lower = 0>[n_survey] n_trap_m1 = trap_count - 1;
  vector<lower = 0>[n_survey_dev] n_trap_m1_dev = trap_count_dev - 1;
  vector[n_survey] log_effort_per = log(effort_per);
  vector[n_survey_dev] log_effort_per_dev = log(effort_per_dev);
  real log_pi = log(pi());
  vector[n_survey] log_survey_area_km2 = log(survey_area_km2);
  vector[n_survey_dev] log_survey_area_km2_dev = log(survey_area_km2_dev);
  vector[n_st] log_area_km2 = log(area_km2);

  for (i in 1:n_survey) first_survey[i] = order[i] == 1;
  for (i in 1:n_survey_dev) first_survey_dev[i] = order_dev[i] == 1;
}

parameters {
  vector[m_p] beta_p;
  vector[m_n] beta;

  // abundance ranefs
  vector[n_property] eps_propertyR;
  real<lower = 0> sigma_property;

  vector[n_property] eps_property_pR;
  real<lower = 0> sigma_property_p;

  // michaelis menten parameters for area sampled as a function of effort
  vector<lower = 0>[2] gamma; // snares and traps
  vector[n_method] log_rho;
  vector<lower = 0, upper = 1>[n_method] p_unique;

  // B-spline spatial coefficients
  matrix[m_short, n_county] z_shortR;
  real<lower = 0> sigma_short;

  // spatiotemporal random effect
  real<lower = 0> sigma_st;
  real<lower = 0> sigma_st0;
  real<lower = 0, upper = 1> eta; // autoregressive param
  matrix[n_timestep, n_county] eps_stR;
}

transformed parameters {
  vector<upper = 0>[n_survey] log_pr_area_sampled;
  vector[n_survey] log_p;
  vector[2] log_gamma = log(gamma);
  matrix[m_short, n_county] z_short;
  vector[n_st] log_lambda;
  vector[n_property] eps_property;

  eps_property = eps_propertyR * sigma_property;

  {
    real log_potential_area;
    for (i in 1:n_survey) {
      if (method[i] > 3) {
        // trap & snare- need to convert radius to area
        log_potential_area = log_pi
          + 2 * (log_rho[method[i]]
                 + log_effort_per[i]
                 - log_sum_exp(log_gamma[method[i] - 3], log_effort_per[i]))
          + log1p(p_unique[method[i]] * n_trap_m1[i]);
      } else {
        // otherwise, for shooting and aerial, the prior is directly on area
        log_potential_area = log_rho[method[i]] + log_effort_per[i]
            + log1p(p_unique[method[i]] * n_trap_m1[i]);
      }
      if (log_potential_area > log_survey_area_km2[i]) {
        log_pr_area_sampled[i] = 0;
      } else {
        log_pr_area_sampled[i] = log_potential_area - log_survey_area_km2[i];
      }
    }
  }

  {
    vector[n_survey] log_theta;
    log_theta = log_inv_logit(X_p * beta_p
        + eps_property_pR[p_property_idx] * sigma_property_p) // inv_logit(...) = pr_det
                + log_pr_area_sampled;
    for (i in 1:n_survey) {
      if (first_survey[i]) {
        log_p[i] = log_theta[i];
      } else {
        log_p[i] = log_theta[i] + sum(log1m_exp(log_theta[start[i]:end[i]]));
      }
    }
  }


  // log(E(# wild pigs))
  z_short[1, ] = z_shortR[1, ];
  for (i in 2:m_short) z_short[i, ] = z_shortR[i - 1, ] + z_shortR[i, ] * sigma_short;

  {
    matrix[n_timestep, n_county] eps_st;
    eps_st[1] = eps_stR[1] * sigma_st0;
    for (t in 2:n_timestep) {
      eps_st[t] = eta * eps_st[t - 1] + eps_stR[t] * sigma_st;
    }

    for (i in 1:n_st) {
      // add if statement for X!=NA
      for(j in 1:m_n){
        if(X_use[i,j] == 1){
                log_lambda[i] = X[i, ] * beta  // I might need to add j indices on these
                      + X_short[i, ] * z_short[, county_idx[i]]
                      + eps_property[property[i]]
                      + eps_st[timestep[i]][county_idx[i]]
                      + log_area_km2[i];
        }
      }
    }
  }
}

model {
  // priors
  beta ~ normal(0, 1);
  beta_p ~ normal(0, 1);
  p_unique ~ beta(2, 2);
  sigma_property ~ normal(0, 1);
  eps_propertyR ~ normal(0, 1);
  sigma_property_p ~ normal(0, 1);
  eps_property_pR ~ normal(0, 1);


  sigma_st ~ normal(0, 1);
  sigma_st0 ~ normal(0, 1);
  eta ~ beta(9, 1);
  for (i in 1:n_timestep) eps_stR[i, ] ~ normal(0, 1);

  sigma_short ~ normal(0, 1);
  for (i in 1:m_short) {
    target += -0.5 * dot_self(z_shortR[i][node1] - z_shortR[i][node2]);
    for (j in 1:n_island) {
      z_shortR[i, island_idx[j]] ~ normal(0, 1);
    }
    sum(z_shortR[i]) ~ normal(0, .001 * n_county);
  }

  gamma[1] ~ gamma(7.704547, 4.41925);         // snare
  gamma[2] ~ gamma(3.613148, 3.507449);        // trap

  log_rho[1] ~ normal(-1.654013, 0.8092821);   // firearms
  log_rho[2] ~ normal(2.74274, 0.5976205);     // fixed wing
  log_rho[3] ~ normal(2.74274, 0.5976205);     // helicopter
  log_rho[4] ~ normal(0.2777449, 0.2255966);   // snare
  log_rho[5] ~ normal(0.6438721, 0.2225168);   // trap

  // likelihood
  y ~ poisson_log(log_lambda[survey_idx] + log_p);
}

generated quantities {
  vector[n_survey_dev] log_pr_area_sampled_dev;
  int y_rep[n_survey]; // is this predicted abundance? or harvest
  vector[n_survey] loglik_train;
  int y_rep_dev[n_survey_dev];
  vector[n_survey_dev] loglik_dev;
  vector[n_survey_dev] log_p_dev;
  vector[n_survey_dev] eps_surveyR_dev;

  {
    real log_potential_area;
    for (i in 1:n_survey_dev) {
      if (method_dev[i] > 3) {
        // trap & snare- need to convert radius to area
        log_potential_area = log_pi
          + 2 * (log_rho[method_dev[i]]
                 + log_effort_per_dev[i]
                 - log_sum_exp(log_gamma[method_dev[i] - 3], log_effort_per_dev[i]))
          + log1p(p_unique[method_dev[i]] * n_trap_m1_dev[i]);
      } else {
        // otherwise, for shooting and aerial, the prior is directly on area
        log_potential_area = log_rho[method_dev[i]] + log_effort_per_dev[i]
            + log1p(p_unique[method_dev[i]] * n_trap_m1_dev[i]);
      }
      if (log_potential_area > log_survey_area_km2_dev[i]) {
        log_pr_area_sampled_dev[i] = 0;
      } else {
        log_pr_area_sampled_dev[i] = log_potential_area - log_survey_area_km2_dev[i];
      }
    }
  }

  {
    vector[n_survey_dev] log_theta_dev;
    real loglam; // posterior predicted abundance

    log_theta_dev = log_inv_logit(X_p_dev * beta_p
            + eps_property_pR[p_property_idx_dev] * sigma_property_p) // pr_det
                    + log_pr_area_sampled_dev;

    for (i in 1:n_survey_dev) {
      if (first_survey_dev[i]) {
        log_p_dev[i] = log_theta_dev[i];
      } else {
        log_p_dev[i] = log_theta_dev[i]
                       + sum(log1m_exp(log_theta_dev[start_dev[i]:end_dev[i]]));
      }
    }

    for (i in 1:n_survey) {
      loglam = log_lambda[survey_idx[i]] + log_p[i];
      y_rep[i] = poisson_log_rng(loglam);
      loglik_train[i] = poisson_log_lpmf(y[i] | loglam);
    }
    for (i in 1:n_survey_dev) {
      loglam = log_lambda[survey_idx_dev[i]] + log_p_dev[i];
      y_rep_dev[i] = poisson_log_rng(loglam);
      loglik_dev[i] = poisson_log_lpmf(y_dev[i] | loglam);
    }
  }
}
