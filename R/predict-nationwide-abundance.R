# srun --pty --account="iwctml" -t 0-4:00 --mem=128000 /bin/bash
# module load swset/2018.05  gcc/7.3.0 r/3.5.2-py27 r-rgdal/1.2-16-py27 r-sf/0.5-5-py27 r-rcpp/1.0.0-py27

#module load swset/2018.05  gcc/7.3.0 r/4.0.5-py27 #r-rstan/2.18.2-py27 r-sf/0.5-5-py27 r-rgdal/1.2-16-py27 r-rcpp/1.0.2-py27

local=FALSE
if(local){
  setwd("~/Desktop/APHIS/abundanceModeling/wild-pigs")
  library(raster)
  library(lme4)
  library(reshape2) 
  library(rstan)
  library(splines)
  library(tidyverse)
  library(readxl)
  library(sf)
  library(assertthat)
  library(data.table)
}else{
  setwd("/project/iwctml/mtabak/APHIS/abundance/wild-pigs")
  library(lme4)
  library(reshape2) 
  library(rstan)
  library(splines)
  library(tidyverse)
  library(readxl)
  library(assertthat)
}


# Simulating posterior draws for abundance across the U.S. -----------------

# 1. Load up the NLCD data
nlcd <- read_csv('./data/out/nlcd-counties.csv')
colnames(nlcd) <- c("n", "cover_class", "year","FIPS")

nlcd_classes <- count(nlcd, cover_class) %>%
  arrange(cover_class)

class_dict <- c('unknown' = 0, 
                'open_water' = 11, 
                'perennial_ice/snow' = 12, 
                'developed_open' = 21, 
                'developed_low' = 22, 
                'developed_med' = 23, 
                'developed_high' = 24, 
                'barren' = 31,
                'deciduous' = 41, 
                'evergreen' = 42, 
                'mixed_forest' = 43, 
                'shrub/scrub' = 52, 
                'grassland' = 71, 
                'pasture' = 81, 
                'cultivated' = 82, 
                'woody_wetland' = 90, 
                'emergent_herbaceous_wetlands' = 95)

nlcd_classes <- nlcd_classes %>%
  mutate(class_name = names(class_dict)[match(nlcd_classes$cover_class, 
                                              class_dict)])

nlcd <- nlcd %>%
  left_join(nlcd_classes, by="cover_class")

# determine which habitats may be pig habitat
valid_habitat <- c('developed_open', 
                   'deciduous', 
                   'evergreen', 
                   'mixed_forest', 
                   'shrub/scrub', 
                   'grassland', 
                   'pasture', 
                   'cultivated')

# nlcd <- nlcd %>%
#   mutate(is_valid = class_name %in% valid_habitat)
nlcd$is_valid <- ifelse(nlcd$class_name %in% valid_habitat, TRUE, FALSE)

# this is not working properly - I do it properly below
# nlcd_summary <- nlcd %>%
#   group_by(year, FIPS) %>%
#   summarize(pr_valid = mean(is_valid)) %>%
#   group_by(FIPS) %>%
#   summarize(pr_valid = mean(pr_valid))

#nlcd_sum <- aggregate(nlcd$is_valid ~ nlcd$FIPS, FUN=mean)

# why isnt this working
#quantile(nlcd_summary$pr_valid)
#quantile(nlcd_sum[,2])
#tmp <- cbind(nlcd$class_name, nlcd$is_valid)
#tmp[100:150,]

# try doing this with base
unique_FIPS <- unique(nlcd$FIPS)
years <- unique(nlcd$year)
nlcd_sum <- matrix(NA, nrow=length(unique_FIPS), ncol=2)
for(i in seq_along(unique_FIPS)){
    #d1 <- nlcd[nlcd$FIPS==i & nlcd$year==j, ]
    d1 <- nlcd[nlcd$FIPS==unique_FIPS[i], ]
    n_tot <- sum(d1$n.x)
    d2 <- d1[d1$is_valid == TRUE, ]
    pr_valid <- sum(d2$n.x)/n_tot
    nlcd_sum[i,2] <- pr_valid
    nlcd_sum[i,1] <- unique_FIPS[i]
}
nlcd_summary <- data.frame(nlcd_sum[,1], as.numeric(nlcd_sum[,2]))
colnames(nlcd_summary) <- c("FIPS", "pr_valid")

fips_df <- data.frame("fips_reduced" = 1:length(unique(nlcd$FIPS)),
"fips_full" =  unique_FIPS)

# 2. Compute the total area of potential habitat for each county
 # shp <- read_sf('data/counties/dtl.cnty.lower48.meters.shp') %>%
 #   left_join(read_rds('process_covs.rds')) %>%
 #   mutate(area_km2 = AREA_GEO, 
 #          id = 1:n()) %>%
 #   left_join(nlcd_summary, by="FIPS") %>%
 #   mutate(valid_km2 = area_km2 * pr_valid)

#^this doesn't work. do it without all the tidyr
#*** do this separately, can't load sf and tidyverse
if(local){
  shp <- read_sf('data/counties/dtl.cnty.lower48.meters.shp') %>%
    left_join(read_rds('process_covs.rds')) %>%
    left_join(nlcd_summary, by="FIPS") 
  shp$area_km2 <- shp$AREA_GEO
  shp$valid_km2 <- shp$area_km2 *shp$pr_valid
  shp$id <- 1:nrow(shp)
  
  shp <- as.data.frame(shp)
  saveRDS(shp, 'shp.rds')
} else{
  shp <- read_rds('shp.rds')
}


st_d <- read_rds('st_d.rds')


# 3. Write a function that allocates total area to properties
property_model <- rstan::extract(read_rds('property-size-fit.rds'))

simulate_properties <- function(df_row, property_model, n_iter) {
  # for one county, simulate property area vectors from the posterior of the 
  # property area model
  assert_that(nrow(df_row) == 1)
  mu_area <- property_model$mu_area[1:n_iter, df_row$id]
  property_areas <- vector(mode = 'list', length = n_iter)
  n_property <- rep(NA, n_iter)
  for (i in 1:n_iter) {
    # adjust for the fact that the valid km^2 is a maximum area: some fraction
    # of this will actually be occupied. 
    fraction_occupied <- rbeta(1, 1, 2)
    area_to_allocate <- df_row$valid_km2 * fraction_occupied
    area_left <- area_to_allocate
    property_areas[[i]] <- NA
    counter <- 0
    while(area_left > 0) {
      counter <- counter + 1
      potential_value <- exp(rnorm(1, 
                                   mu_area[i], 
                                   property_model$sigma_area[i]))
      if (potential_value <= area_left) {
        property_areas[[i]][counter] <- potential_value
      } else {
        property_areas[[i]][counter] <- area_to_allocate
      }
      area_left <- area_left - property_areas[[i]][counter]
    }
    assert_that(area_to_allocate - sum(property_areas[[i]]) < 1e-10)
  }
  return(property_areas)
}

# Load pig abundance model posterior
post <- rstan::extract(readRDS('fitted_mods/saturating_fit_VA_20210413.rds'),
#post <- rstan::extract(readRDS('saturating_fit.rds'),                       
                       pars = c('sigma_property', 'sigma_st0', 'sigma_st', 
                                'eps_stR', 'eta', 'beta', 'z_short', 
                                'eps_property'))

# Load data required for fixed effects etc.
stan_d <- read_rds('stan_d.rds')
X_full <- model.matrix(~ 0 +
                         c_hydroden + 
                         c_lewis +
                         c_hetero + 
                         c_carnrich +
                         c_crop + 
                         c_pasture + 
                         evergreen +
                         deciduous +
                         Cereals +
                         c_tree 
                         , 
                       data = shp)
assert_that(nrow(X_full) == nrow(shp))
#X_full <- X_full[, -c(7:9)]

short_basis <- read_rds('short_basis.rds')
timestep_df <- read_rds('timestep_df.rds')
assert_that(all(timestep_df$timestep == 1:nrow(timestep_df)))
X_short <- predict(short_basis, newx = timestep_df$timestep)

# save progress for easier loading
#save.image("workspace_20190821.RData")

# Function to take one instatiation of a property area vector and simulate
# abundance over time
sim_abundance <- function(shp_row, post, X_full, shp, n_iter = 10) {
  property_areas <- simulate_properties(shp[shp_row, ], property_model, n_iter)
  
  # get property adjustments for log lambda
  n_properties <- lapply(property_areas, length) %>%
    unlist
  
  eps_prop_df <- tibble(n_properties = n_properties, 
                        property_areas = property_areas, 
                        sigma_property = post$sigma_property[1:n_iter]) %>%
    rowwise() %>%
    mutate(eps_property_plus_offset = list(
      rnorm(n_properties, log(property_areas), sigma_property))) %>%
    ungroup
  
  fixef_vals <- tibble(iter = 1:n_iter, 
                       fixef = c(X_full[shp_row, ] %*% t(post$beta[1:n_iter, ])))

  eps_st <- matrix(nrow = n_iter, ncol = stan_d$n_timestep)
  eps_stR <- post$eps_stR[1:n_iter, , shp_row]
  eps_st[, 1] = eps_stR[, 1] * post$sigma_st0[1:n_iter]
  for (t in 2:stan_d$n_timestep) {
    eps_st[, t] = post$eta[1:n_iter] * eps_st[, t - 1] + 
                    eps_stR[, t] * post$sigma_st[1:n_iter]
  }
  
  timestep_adj <- X_short %*% t(post$z_short[1:n_iter, , shp_row]) %>%
    melt(varnames = c('timestep', 'iter'), value.name = 'bspline_adj') %>%
    as_tibble %>%
    arrange(iter, timestep) %>%
    bind_cols(melt(eps_st, varnames = c('iter', 'timestep'), 
                   value.name = 'eps_st')) %>%
    mutate(adj = bspline_adj + eps_st) %>%
    select(timestep, iter, adj)
  
  total_abundance <- eps_prop_df %>%
    bind_cols(fixef_vals) %>%
    mutate(iter = 1:n_iter) %>%
    full_join(timestep_adj) %>%
    mutate(log_lam_part = fixef + adj) %>%
    rowwise() %>%
    mutate(total_pigs = sum(rpois(n_properties, 
                                   exp(log_lam_part + 
                                         eps_property_plus_offset)))) %>%
    ungroup %>%
    filter(!is.na(total_pigs)) 
  total_abundance
}

#- find properties that are sampled so that we can only estimate density in these
# counties
sampled_subset <- st_d %>%
  group_by(property_factor) %>%
  mutate(n_obs = n()) %>%
  filter(n_obs >= 20)
FIPS_use <- unique(sampled_subset$FIPS)

# Function to take one instatiation of a property area vector and simulate
# density over time
sim_abundance <- function(shp_row, post, X_full, shp, n_iter = 10) {
  
  # only do this function if the county is in the sampled_subset
  #if(shp[shp_row, "FIPS"]$FIPS %in% FIPS_use){
    property_areas <- simulate_properties(shp[shp_row, ], property_model, n_iter)
    
    # get property adjustments for log lambda
    n_properties <- lapply(property_areas, length) %>%
      unlist
    
    eps_prop_df <- tibble(n_properties = n_properties, 
                          property_areas = property_areas, 
                          sigma_property = post$sigma_property[1:n_iter]) %>%
      rowwise() %>%
      mutate(eps_property_plus_offset = list(
        rnorm(n_properties, log(property_areas), sigma_property))) %>%
      ungroup
    
    fixef_vals <- tibble(iter = 1:n_iter, 
                         fixef = c(X_full[shp_row, ] %*% t(post$beta[1:n_iter, ])))
    
    eps_st <- matrix(nrow = n_iter, ncol = stan_d$n_timestep)
    eps_stR <- post$eps_stR[1:n_iter, , shp_row]
    eps_st[, 1] = eps_stR[, 1] * post$sigma_st0[1:n_iter]
    for (t in 2:stan_d$n_timestep) {
      eps_st[, t] = post$eta[1:n_iter] * eps_st[, t - 1] + 
        eps_stR[, t] * post$sigma_st[1:n_iter]
    }
    
    timestep_adj <- X_short %*% t(post$z_short[1:n_iter, , shp_row]) %>%
      melt(varnames = c('timestep', 'iter'), value.name = 'bspline_adj') %>%
      as_tibble %>%
      arrange(iter, timestep) %>%
      bind_cols(melt(eps_st, varnames = c('iter', 'timestep'), 
                     value.name = 'eps_st')) %>%
      mutate(adj = bspline_adj + eps_st) %>%
      rename(timestep = timestep...1, iter=iter...2) %>%
      select(timestep, iter, adj)
    
    prop_occ <- 1 #rbeta(1,1,10)
    
    total_abundance <- eps_prop_df %>%
      bind_cols(fixef_vals) %>%
      mutate(iter = 1:n_iter) %>%
      full_join(timestep_adj) %>%
      mutate(log_lam_part = fixef + adj) %>%
      rowwise() %>%
      mutate(total_pigs = sum(rpois(n_properties, 
                                    lambda=(exp(log_lam_part + 
                                                 eps_property_plus_offset)*prop_occ)))) %>%
      # convert this to density
      mutate(pig_density = total_pigs/sum(unlist(property_areas))) %>%
      ungroup %>%
      filter(!is.na(total_pigs)) 
    
    total_abundance
  #}
}
#as.data.frame(total_abundance[1729,])


abundance_ts <- vector(mode = 'list', length = nrow(shp))
pb <- txtProgressBar(max = nrow(shp), style = 3)
for (i in 1:nrow(shp)) {
  if (is.null(abundance_ts[[i]])) {
    abundance_ts[[i]] <- sim_abundance(i, post, X_full, shp = shp, 
                                       n_iter = 498)
  }
  setTxtProgressBar(pb, i)
  print(paste('i =', i, " of ", nrow(shp)))
}
warnings()
# remove null elements of list
#abundance_ts2 <- abundance_ts[-which(sapply(abundance_ts, is.null))]
#names(abundance_ts2) <- FIPS_use
#abundance_ts <- abundance_ts2

#write_rds(abundance_ts, path = './data/out/abundance_ts_20200925.rds')
#abundance_ts <- read_rds(path = './data/out/abundance_ts_20190828d.rds')

all_df <- abundance_ts %>%
  bind_rows(.id = "FIPS") %>%
  #bind_rows() %>%
  select(FIPS, timestep, iter, total_pigs, pig_density)
quantile(all_df$pig_density, c(0.025,0.25, 0.5, 0.75, 0.975))
write_csv(all_df, "./data/out/all_df_20210413.csv")

# write rds to fst
library(fst)
write.fst(abundance_ts, './data/out/abundance_ts.fst', 100)
write.fst(all_df, './data/out/all_df_density_20210413.fst', 100)
all_df <- read.fst('./data/out/all_df_fromCSV_20210413.fst')

# testing something here
all_df <- abundance_ts %>%
  bind_rows(.id = "FIPS") %>%
  #bind_rows() %>%
  select(FIPS, timestep, iter, total_pigs, pig_density, property_areas)
quantile(all_df$pig_density, c(0.025,0.25, 0.5, 0.75, 0.975))
saveRDS(all_df, "./data/out/all_df_20210413.rds")



# bind_rows(.id) only works sometimes, so do this with base
all_df$fips <- unique_FIPS[as.numeric(all_df$FIPS)]

test <- all_df[sample(nrow(all_df),100),]
test$fips <- unique_FIPS[as.numeric(test$FIPS)]

write_rds(all_df, path = './data/out/all_df_density_20210413.rds')
#all_df <- read_rds(path = './data/out/all_df_density_20190828d.rds')

# save workspace because all_df is messing up the FIPS
#save.image("workspace_201908229.RData")

# write rds to fst
library(fst)
write.fst(all_df, './data/out/all_df_density_20210413.fst', 100)

# for density
total_abundance_summary <- all_df %>%
  group_by(timestep, iter) %>%
  summarize(dens = median(pig_density)) %>%
  group_by(timestep) %>%
  filter(!is.na(dens)) %>%
  summarize(med = median(dens), 
            lo = quantile(dens, .025), 
            hi = quantile(dens, .975)
            #lo = quantile(total, .49), 
            #hi = quantile(total, .51)
  )
cnty_abundance_summary <- all_df %>%
  group_by(FIPS, iter) %>%
  summarize(abun = median(total_pigs)) %>%
  group_by(timestep) %>%
  filter(!is.na(dens)) %>%
  summarize(med = median(abun), 
            lo = quantile(abun, .025), 
            hi = quantile(abun, .975)
            #lo = quantile(total, .49), 
            #hi = quantile(total, .51)
  )

pdf('./fig/total-density_VA_20210413.pdf', width = 6, height = 3)
total_abundance_summary %>%
  left_join(timestep_df) %>%
  ggplot(aes(start_dates, med)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .3) +
  ylab('National-scale wild pig density') + 
  theme_minimal() + 
  xlab('')
dev.off()


# for abundance
{
  total_abundance_summary <- all_df %>%
    group_by(timestep, iter) %>%
    summarize(total = sum(total_pigs)) %>%
    group_by(timestep) %>%
    filter(!is.na(total)) %>%
    summarize(med = median(total), 
              lo = quantile(total, .025), 
              hi = quantile(total, .975)
              #lo = quantile(total, .49), 
              #hi = quantile(total, .51)
    )
}
pdf('./fig/total-abundance_VA_20210413.pdf', width = 6, height = 3)
total_abundance_summary %>%
  left_join(timestep_df) %>%
  ggplot(aes(start_dates, med)) +
  geom_line() + 
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .3) +
  ylab('Total wild pig abundance') + 
  theme_minimal() + 
  xlab('')
dev.off()
#ggsave(filename = 'fig/total-abundance.png', width = 6, height = 3)

#-- test something
all_df$dens2 <- all_df$total_pigs/sum(unlist(all_df[,3][[1]]))

tmp.shit <- apply(all_df, 1, function(x) (all_df$total_pigs/sum(unlist(all_df[x,3][[1]]))))
tmp <- apply(all_df, 1, function(x) (x$total_pigs/sum(unlist(x[,3][[1]]))))

#------------------------------------------------------------
# get state data
#------------------------------------------------------------
# library(noncensus)
# state <- zip_codes[as.numeric(all_df$fips) %in% zip_codes$fips , "state"]
# all_df$state <- state
all_df <- read_rds(path = './data/out/all_df_density_20210413.rds')
state_fips <- read.csv("./data/state_fips.csv")
state_fips[,1] <- as.character(state_fips[,1])
state_fips[,2] <- as.character(state_fips[,2])
state_fips[56,] <- c(("District of Columbia"), ("DC"), 11)
state_fips[,3 ] <- as.numeric(state_fips[,3])
all_df$st_code <- as.numeric(substring(all_df$fips, 1,2))


# get a state name into all_df
# try merge: This works!!!
all_df <- merge(all_df, state_fips, by.x=7, by.y=3, all.x=TRUE, all.y=FALSE)
all_df$state <- all_df$Name

# write this out before it kills R
#write_rds(all_df, path = './data/out/all_df_density_20190828e.rds')
all_df <- read_rds(path = './data/out/all_df_density_20210413.rds')

# check how this worked
states <- unique(all_df$state)
uFIPS <- unique(all_df$st_code)
difs_states <- setdiff(state_fips$Name, states)
difs <- setdiff(state_fips$FIPS, uFIPS)
state_fips[which(state_fips$FIPS %in% difs),1] == difs_states

# calculate density by timestep and state
total_density_summary <- all_df %>%
  group_by(timestep, iter, state) %>%
  summarize(dens = median(pig_density)) %>%
  group_by(timestep) %>%
  filter(!is.na(dens)) %>%
  summarize(med = median(dens), 
            lo = quantile(dens, .025), 
            hi = quantile(dens, .975)
            #lo = quantile(total, .49), 
            #hi = quantile(total, .51)
  )
# using aggregate
dens_sum <- aggregate(all_df$pig_density, by=list(all_df$state, all_df$timestep), FUN=
                        median)
ag.lo <- aggregate(all_df$pig_density, by=list(all_df$state, all_df$timestep), FUN=
                   function(x) quantile(x, 0.025))
ag.hi <- aggregate(all_df$pig_density, by=list(all_df$state, all_df$timestep), FUN=
                     function(x) quantile(x, 0.975))
dens_summary <- data.frame(dens_sum, ag.lo[,3], ag.hi[,3])
colnames(dens_summary) <- c("state", "timestep", "median", "lo", "hi")
# add timestep
dens_sum2 <- merge(x=dens_summary, y=timestep_df, by="timestep", all.x=TRUE, all.y=FALSE)
write_csv(dens_sum2, path = './data/out/dens_summary_20210413.csv')

# density by county: using FIPS
dens_sum <- aggregate(all_df$pig_density, by=list(all_df$fips, all_df$timestep), FUN=
                        function(x) quantile(x, c(0.025, 0.5, 0.975))

                      )
dens_summary <- data.frame(dens_sum$Group.1, dens_sum$Group.2, dens_sum[,1], dens_sum[,2], dens_sum[,3])
dens_summary2 <- dens_summary[,-(1:2)]
colnames(dens_summary2) <- c("fips", "timestep", "lo", "med", "hi")
# add timestep
dens_sum2 <- merge(x=dens_summary2, y=timestep_df, by="timestep", all.x=TRUE, all.y=FALSE)
dens_sum2$fips <- as.character(dens_sum2$fips)

# add state
all_state <- all_df[!duplicated(all_df$fips),c(7,10)]
dens_sum3 <- merge(x=dens_sum2, y=all_state, by.x="fips", by.y="fips", all.x=TRUE, all.y=FALSE)

write_csv(dens_sum3, path = './data/out/dens_summary_cnty_20210413.csv')
getwd()
# get property data using loaded workspace
load("workspace_20190829.RData")
str(abundance_ts)

#-----------------------------------------------------
# archive
#-----------------------------------------------------

# this requires 128 GB
lookup1 <- setNames(state_fips$Name, state_fips$FIPS)
#all_df$state <- sapply(all_df$st_code, function(x) lookup1[x])
tmp<- lapply(all_df$st_code, function(x) lookup1[x])
tmp <- apply(all_df$st_code, 1, function(x) lookup1[x])
all_df$state <- unlist(tmp)

# try dplyr
all_df %>%
  gather(key = "st_code") %>%
  left_join(state_fips, by = "st_code") %>%
  spread(key = FIPS, value = NAME)

# not working. Use ifelse
all_df$state <- ifelse(all_df$st_code == 1, "Alabama", 
                       ifelse(all_df$st_code == 2, "Alaska", 
                              ifelse(all_df$st_code == 4, "Arizona",
                                     ifelse(all_df$st_code == 5, "Arkansas", 
                                            ifelse(all_df$st_code == 6, "California", 
                                                   ifelse(all_df$st_code == 8, "Colorado", 
                                                          ifelse(all_df$st_code == 9, "Connecticut", 
                                                                 ifelse(all_df$st_code == 10, "Delaware", 
                                                                        ifelse(all_df$st_code == 12, "Florida", 
                                                                               ifelse(all_df$st_code == 13, "Georgia", 
                                                                                      ifelse(all_df$st_code == 15, "Hawaii", 
                                                                                             ifelse(all_df$st_code == 16, "Idaho", 
                                                                                                    ifelse(all_df$st_code == 17, "Illinois", 
                                                                                                           ifelse(all_df$st_code == 18, "Indiana", 
                                                                                                                  ifelse(all_df$st_code == 19, "Iowa", 
                                                                                                                         ifelse(all_df$st_code == 20, "Kansas", 
                                                                                                                                ifelse(all_df$st_code == 21, "Kentucky", 
                                                                                                                                       ifelse(all_df$st_code == 22, "Louisiana", 
                                                                                                                                              ifelse(all_df$st_code == 23, "Maine", 
                                                                                                                                                     ifelse(all_df$st_code == 24, "Maryland",
                                                                                                                                                            ifelse(all_df$st_code == 25, "Massachusetts", 
                                                                                                                                                                   ifelse(all_df$st_code == 26, "Michigan", 
                                                                                                                                                                          ifelse(all_df$st_code == 27, "Minnesota", 
                                                                                                                                                                                 ifelse(all_df$st_code == 28, "Mississippi", 
                                                                                                                                                                                        ifelse(all_df$st_code == 29, "Missouri", 
                                                                                                                                                                                               ifelse(all_df$st_code == 30, "Montana", 
                                                                                                                                                                                                      ifelse(all_df$st_code == 31, "Nebraska", 
                                                                                                                                                                                                             ifelse(all_df$st_code == 32, "Nevada", 
                                                                                                                                                                                                                    ifelse(all_df$st_code == 33, "New Hampshire", 
                                                                                                                                                                                                                           ifelse(all_df$st_code == 34, "New Jersey", 
                                                                                                                                                                                                                                  ifelse(all_df$st_code == 35, "New Mexico", 
                                                                                                                                                                                                                                         ifelse(all_df$st_code == 36, "New York", 
                                                                                                                                                                                                                                                ifelse(all_df$st_code == 37, "North Carolina", 
                                                                                                                                                                                                                                                       ifelse(all_df$st_code == 38, "North Dakota", 
                                                                                                                                                                                                                                                              ifelse(all_df$st_code == 39, "Ohio", 
                                                                                                                                                                                                                                                                     ifelse(all_df$st_code == 40, "Oklahoma", 
                                                                                                                                                                                                                                                                            ifelse(all_df$st_code == 41, "Oregon", 
                                                                                                                                                                                                                                                                                   ifelse(all_df$st_code == 42, "Pennsylvania", 
                                                                                                                                                                                                                                                                                          ifelse(all_df$st_code == 44, "Rhode Island", 
                                                                                                                                                                                                                                                                                                 ifelse(all_df$st_code == 45, "South Carolina", 
                                                                                                                                                                                                                                                                                                        ifelse(all_df$st_code == 46, "South Dakota", 
                                                                                                                                                                                                                                                                                                               ifelse(all_df$st_code == 47, "Tennessee", 
                                                                                                                                                                                                                                                                                                                      ifelse(all_df$st_code == 48, "Texas", 
                                                                                                                                                                                                                                                                                                                             ifelse(all_df$st_code == 49, "Utah", 
                                                                                                                                                                                                                                                                                                                                    ifelse(all_df$st_code == 50, "Vermont", 
                                                                                                                                                                                                                                                                                                                                           ifelse(all_df$st_code == 51, "Virginia", 
                                                                                                                                                                                                                                                                                                                                                  ifelse(all_df$st_code == 53, "Washington", 
                                                                                                                                                                                                                                                                                                                                                         ifelse(all_df$st_code == 54, "West Virginia", 
                                                                                                                                                                                                                                                                                                                                                                ifelse(all_df$st_code == 55, "Wisconsin", 
                                                                                                                                                                                                                                                                                                                                                                       ifelse(all_df$st_code == 56, "Wyoming", 
                                                                                                                                                                                                                                                                                                                                                                              ifelse(all_df$st_code == 60, "American Samoa", 
                                                                                                                                                                                                                                                                                                                                                                                     ifelse(all_df$st_code == 66, "Guam", 
                                                                                                                                                                                                                                                                                                                                                                                            ifelse(all_df$st_code == 69, "Northern Mariana Islands", 
                                                                                                                                                                                                                                                                                                                                                                                                   ifelse(all_df$st_code == 72, "Puerto Rico", 
                                                                                                                                                                                                                                                                                                                                                                                                          ifelse(all_df$st_code == 78, "Virgin Islands", 
                                                                                                                                                                                                                                                                                                                                                                                                                 ifelse(all_df$st_code == 11, "District of Columbia", NA)
                                                                                                                                                                                                                                                                                                                                                                                                          )
                                                                                                                                                                                                                                                                                                                                                                                                   )
                                                                                                                                                                                                                                                                                                                                                                                            )
                                                                                                                                                                                                                                                                                                                                                                                     )
                                                                                                                                                                                                                                                                                                                                                                              )
                                                                                                                                                                                                                                                                                                                                                                       )
                                                                                                                                                                                                                                                                                                                                                                )
                                                                                                                                                                                                                                                                                                                                                         )
                                                                                                                                                                                                                                                                                                                                                  )
                                                                                                                                                                                                                                                                                                                                           )
                                                                                                                                                                                                                                                                                                                                    )
                                                                                                                                                                                                                                                                                                                             )
                                                                                                                                                                                                                                                                                                                      )
                                                                                                                                                                                                                                                                                                               )
                                                                                                                                                                                                                                                                                                        )
                                                                                                                                                                                                                                                                                                 )
                                                                                                                                                                                                                                                                                          )
                                                                                                                                                                                                                                                                                   )
                                                                                                                                                                                                                                                                            )
                                                                                                                                                                                                                                                                     )
                                                                                                                                                                                                                                                              )
                                                                                                                                                                                                                                                       )
                                                                                                                                                                                                                                                )
                                                                                                                                                                                                                                         )
                                                                                                                                                                                                                                  )
                                                                                                                                                                                                                           )
                                                                                                                                                                                                                    )
                                                                                                                                                                                                             )
                                                                                                                                                                                                      )
                                                                                                                                                                                               )
                                                                                                                                                                                        )
                                                                                                                                                                                 )
                                                                                                                                                                          )
                                                                                                                                                                   )
                                                                                                                                                            )
                                                                                                                                                     )
                                                                                                                                              )
                                                                                                                                       )
                                                                                                                                )
                                                                                                                         )
                                                                                                                  )
                                                                                                           )
                                                                                                    )
                                                                                             )
                                                                                      )
                                                                               )
                                                                        )
                                                                 )
                                                          )
                                                   )
                                            )
                                     )
                              )
                       )
)


