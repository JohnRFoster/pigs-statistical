# srun --pty --account="iwctml" -t 0-2:00 --mem=128000 /bin/bash
# module load swset/2018.05  gcc/7.3.0 r/3.5.2-py27 r-rgdal/1.2-16-py27 r-sf/0.5-5-py27 r-rcpp/1.0.0-py27

library(readr)
library(lubridate)
library(maps)

libs <- c("rgdal", "maptools", "gridExtra","classInt","Rcpp","raster","maps","mapdata",
          "dismo","rgeos", "RColorBrewer")
lapply(libs, require, character.only = TRUE)


setwd("/project/iwctml/mtabak/APHIS/abundance/wild-pigs")
setwd("~/Desktop/APHIS/abundanceModeling/wild-pigs")

# read in data
st_d <- read_rds('st_d.rds')
m_fit <- read_rds("./fitted_mods/saturating_fit_20190416.rds")
fit <- "fit20190416" # this is the location name in ./fitted_mods/
stan_d <- read_rds("stan_d.rds")
timestep_df <- read_rds('timestep_df.rds')
cnty_dens_meth <- "weighted"

# extract abunaces
log_lambda <- data.frame(rstan::extract(m_fit, "log_lambda"))
# for every row in st_d, there is a column in log_lambda that represents the 
# posterior distribution of the abundance for that spatio temporal unit

# convert to densities for these
dens <- matrix(NA, nrow(log_lambda), ncol(log_lambda))
lam <- exp(log_lambda) # abundance (not on log scale)
for(i in 1:ncol(lam)){
  dens[,i] <- lam[,i]/st_d$property_area_km2[i]
}
#dens[,100] == lam[,100] /  st_d$property_area_km2[100] # checking on this

###*** Temporary
###* Temporaility use only lambda (abundance) instead of density
dens <- lam

#- calculate county-level density
# first get the row indeces in st_d
st_index <- split(seq_along(st_d$FIPS_timestep), st_d$FIPS_timestep)
n_st <- length(st_index)
sts <- unique(st_d$FIPS_timestep) # these are the FIPS_timesteps that correspond to the st_means

st_means <- matrix(NA, nrow(log_lambda), ncol=length(unique(st_d$FIPS_timestep)))
cnty_dens <- matrix(NA, nrow(log_lambda), ncol=length(unique(st_d$FIPS_timestep)))
if(cnty_dens_meth=="weighted"){
  # using weighted mean
  for(i in 1:n_st){
    if(length(st_index[[i]]) > 1){
      density_cols <- dens[,st_index[[i]][1:length(st_index[[i]])]]
      areas <- st_d$property_area_km2[st_index[[i]][1:length(st_index[[i]])]]
      weights <- areas/sum(areas)
      tmpMeans <- matrix(NA, nrow(density_cols), ncol(density_cols))
      for(j in 1:length(weights)){
        tmpMeans[, j] <- density_cols[,j] * weights[j]
      }
      st_means[,i] <- rowMeans(tmpMeans) 
    } else {
      st_means[,i] <- dens[, st_index[[i]][1]]
    }
  }
} else {
  # using just the mean
  for(i in 1:n_st){
    if(length(st_index[[i]]) > 1){
      st_means[,i] <- rowMeans(dens[,st_index[[i]][1:length(st_index[[i]])]]) #apply(dens, 2, function(x){rowMeans(x[st_index[[i]][1:length(st_index[[i]])]])}) 
    } else {
      st_means[,i] <- dens[, st_index[[i]][1]]
    }
  }
}
# st_means is the mean density at each spatio-temporal unit
# the st units are defined by sts
cnty_dens <- t(st_means) # here st units are rows and columns are MCMC iterations
#write_csv(data.frame(cnty_dens), "./fitted_mods/fit20190416/cnty_dens_posteriors.csv")
#write_csv(data.frame(cnty_dens), "./fitted_mods/fit20190416/cnty_abundance_posteriors.csv")


# downsize st_d so there is just one row for every st unit
# (removing duplicate rows for mulitple surveys because these are now accounted for in st_means)
cnty_d <- st_d[!duplicated(st_d$FIPS_timestep), ]
# mean and median posterior density for each st unit
cnty_d$meanDens <- apply(cnty_dens, 1, mean) 
cnty_d$medDens <- apply(cnty_dens, 1, function(x) quantile(x, 0.5))
cnty_d$loDens <- apply(cnty_dens, 1, function(x) quantile(x, 0.025))
cnty_d$hiDens <- apply(cnty_dens, 1, function(x) quantile(x, 0.975))

# write out cnty_d
#write.csv(cnty_d, "./fitted_mods/fit20190416/cnty_data_abundance.csv")
#cnty_d <- read_csv("./fitted_mods/fit20190416/cnty_data.csv")

# get densities over time
ag1 <- aggregate(medDens ~ end_dates, data=cnty_d, FUN=median)
ag2 <- aggregate(loDens ~ end_dates, data=cnty_d, FUN=median)
ag3 <- aggregate(hiDens ~ end_dates, data=cnty_d, FUN=median)
all(ag1$end_dates == ag2$end_dates)
all(ag1$end_dates == ag3$end_dates)

dens_sum <- data.frame(ag1, loDens=ag2$loDens, hiDens=ag3$hiDens)

# plot density over time
#pdf("./fitted_mods/fit20190416/density_unweighted_20190821.pdf", width = 6, height = 6)
plot(dens_sum$end_dates, dens_sum$medDens, ylab="national wild pig density (pigs/km^2)",
     ylim=c(0, max(dens_sum$hiDens)), type="l", lty=1, lwd=2
)
lines(dens_sum$end_dates, dens_sum$loDens, lty=2)
lines(dens_sum$end_dates, dens_sum$hiDens, lty=2)
dev.off()

# write density summary
#write.csv(dens_sum, "./fitted_mods/fit20190416/density_summary_weighted_20190821.csv")

#- try using whole posteriors in summaries
cnty_d$year <- year(cnty_d$end_dates)
cnty_d$month <- month(cnty_d$end_dates)
cnty_d$mo <- ifelse(nchar(cnty_d$month) == 2, cnty_d$month, paste0("0", cnty_d$month))
cnty_d$ym <- paste0(cnty_d$year, cnty_d$mo)
cnty_post <- data.frame(cnty_d, cnty_dens)
# none of that works. Use a loop
ends <- unique(cnty_post$end_dates)
dens_sum2 <- matrix(NA, length(ends), 4)
for(i in seq_along(ends)){
  tmp <- cnty_post[cnty_post$end_dates==ends[i], (ncol(cnty_d) +1):ncol(cnty_post)]
  dens_sum2[i,2:4] <- quantile(unlist(tmp), c(0.025, 0.5, 0.975))
}
dens_sum2[,1] <- as.character(ends)
write.csv(dens_sum2, "./fitted_mods/fit20190416/density_summary_20190822.csv")

#- aggregate by month
ends <- unique(cnty_post$ym)
dens_sum3 <- matrix(NA, length(ends), 4)
for(i in seq_along(ends)){
  tmp <- cnty_post[cnty_post$ym==ends[i], (ncol(cnty_d) +1):ncol(cnty_post)]
  dens_sum3[i,2:4] <- quantile(unlist(tmp), c(0.025, 0.5, 0.975))
}
dens_sum3[,1] <- as.character(ends)

write.csv(dens_sum3, "./fitted_mods/fit20190416/density_byMonth_20190822.csv")

#- make a map with posterior median density for each county
library(maps)
library(dplyr)
data(county.fips)

# TMP: read in cnty data
#cnty_d <- read.csv("./data/cnty_d_20190410.csv")

# split up data by temporal units
#year(cnty_d$start_dates)
cnty_d$year <- substr(cnty_d$end_dates, start = 1, stop = 4)
cnty_d$month <- substr(cnty_d$end_dates, start = 6, stop = 7)
cnty_d2 <- cnty_d[!duplicated(data.frame(cnty_d$FIPS, cnty_d$year)), ]
yrs <- as.numeric(unique(cnty_d2$year))
for(i in yrs){
  assign(paste0("d", i, "1"), eval(parse(text=paste0("cnty_d2[cnty_d2$year == ", i, " ,]"))))
  assign(paste0("d", i), eval(parse(text=paste0("d", i, "1[!is.na(d", i, "1$FIPS),]"))))
  len <- eval(parse(text=paste0("nrow(d", i, ")")))
  cat(paste0(i, " has ", len, " records \n"))
}


#------------------------------------------------------------------------------
#--- DEVIANCE: look at model fit
#------------------------------------------------------------------------------
# loo_train <- loo(m_fit, "loglik_train")
# saveRDS(loo_train, "loo_train_20190411.rds")
# loo_dev <- loo(m_fit, "loglik_dev")
# saveRDS(loo_dev, "loo_dev_20190411.rds")

# extract deviance
dev_d <- readRDS("dev_d.rds")$group_d
group_d <- readRDS("dev_d.rds")$group_d
y_dev_hat <- data.frame(rstan::extract(m_fit, "y_rep_dev"))
y_dev <- stan_d$y_dev
y_dev_mat <- t(replicate(nrow(y_dev_hat), y_dev)) # make this for elementwise computation
dev <- (y_dev_mat - y_dev_hat)/(y_dev_mat+1) # need to add 1 to denominator because sometimes its 0
meds <- apply(dev,2, median)
dev_d$dev <- meds # median deviance
dev_d$y_hat <- apply(y_dev_hat, 2, median)
dev_d$dev2 <- (dev_d$y - dev_d$y_hat)/(dev_d$y +1)

#--- link deviance up with COUNTY DATA
survey_d <- readRDS("survey_d.rds")
# determine number of surveys in each county
n_samp_cnty <- as.data.frame(table(survey_d$FIPS))
colnames(n_samp_cnty)
FIPS_idx <- split(seq_along(survey_d$FIPS), survey_d$FIPS) # row indeces in survey_d for each FIPS code
n_samp_cty <- lengths(FIPS_idx) # number of samples per county
# data from dev dataset
dev_idx <- split(seq_along(group_d$FIPS), group_d$FIPS) # these are the columns I want for each county in the dev posteriors
dev_cty <- unique(group_d$FIPS)
cty_devs <- matrix(NA, nrow=length(dev_cty), ncol=7) # these are the quantiles of deviance for each county
cty_devs[,1] <- dev_cty
for(i in dev_cty){
  post_cty <- as.vector(as.matrix(data.frame(dev[, eval(parse(text=paste0("dev_idx$`", i, "`")))])))
  #cty_devs[which(i %in% dev_cty), 2] <- quantile(post_cty, c(0.025))
  cty_devs[match(i,dev_cty), 2:6] <- quantile(post_cty, c(0.025, 0.25, 0.5, 0.75, 0.975))
}
# now try breaking this out by n_samp_cty
cty_devs[,7] <- n_samp_cty_dev <- n_samp_cty[match(dev_cty, n_samp_cnty[,1])] # number of samples per county for county in the deviance dataset
devs_bySamp <- matrix(NA, nrow=length(unique(n_samp_cty_dev)), ncol=7) 
devs_bySamp.025 <- aggregate(as.numeric(cty_devs[,2]) ~ as.numeric(cty_devs[,7]), FUN=mean)
devs_bySamp.25 <- aggregate(as.numeric(cty_devs[,3]) ~ as.numeric(cty_devs[,7]), FUN=mean)
devs_bySamp.5 <- aggregate(as.numeric(cty_devs[,4]) ~ as.numeric(cty_devs[,7]), FUN=mean)
devs_bySamp.75 <- aggregate(as.numeric(cty_devs[,5]) ~ as.numeric(cty_devs[,7]), FUN=mean)
devs_bySamp.975 <- aggregate(as.numeric(cty_devs[,6]) ~ as.numeric(cty_devs[,7]), FUN=mean)
devs_bySamp <- data.frame(devs_bySamp.025, devs_bySamp.25[,2], devs_bySamp.5[,2],
                          devs_bySamp.75[,2], devs_bySamp.975[,2])
colnames(devs_bySamp) <- c("n_surveys", "deviance.025", "deviance.25", "deviance.5", 
                           "deviance.75", "deviance.975")
# make plot
thick95 <- 1
#postscript("dev_county_20190411.eps")
pdf(paste0("./fitted_mods/", fit, "/dev_county.pdf"))
plot(log(devs_bySamp$n_surveys), abs(devs_bySamp$deviance.5), pch=16,
     xlab="log number of surveys in county",
     ylab="deviance", axes=FALSE)
#segments(x0=log(devs_bySamp$n_surveys), y0=devs_bySamp$deviance.025, y1=devs_bySamp$deviance.975,
#         lwd=thick95, col="grey")
axis(1)
axis(2, las=2)
dev.off()

#---- estimate deviance by PROPERTY DATA
n_samp_prop <- as.data.frame(table(survey_d$agrp_prp_id))
prp_idx <- split(seq_along(survey_d$agrp_prp_id), survey_d$agrp_prp_id) # row indeces in survey_d for each property
n_samp_prp <- lengths(prp_idx) # number of samples per property
# get data from dev dataset and posteriors
dev_idx <- split(seq_along(group_d$agrp_prp_id), group_d$agrp_prp_id) # these are the columns I want for each county in the dev posteriors
dev_prp <- unique(group_d$agrp_prp_id)
prp_devs <- matrix(NA, nrow=length(dev_prp), ncol=7) # these are the quantiles of deviance for each county
prp_devs[,1] <- dev_prp
for(i in dev_prp){
  post_prp <- as.vector(as.matrix(data.frame(dev[, eval(parse(text=paste0("dev_idx$`", i, "`")))])))
  prp_devs[match(i,dev_prp), 2:6] <- quantile(post_prp, c(0.025, 0.25, 0.5, 0.75, 0.975))
}
# break out by property
prp_devs[,7] <- n_samp_prp_dev <- n_samp_prp[match(dev_prp, n_samp_prop[,1])] # number of samples per county for county in the deviance dataset
devs_bySamp <- matrix(NA, nrow=length(unique(n_samp_prp_dev)), ncol=7) 
devs_bySamp.025 <- aggregate(as.numeric(prp_devs[,2]) ~ as.numeric(prp_devs[,7]), FUN=mean)
devs_bySamp.25 <- aggregate(as.numeric(prp_devs[,3]) ~ as.numeric(prp_devs[,7]), FUN=mean)
devs_bySamp.5 <- aggregate(as.numeric(prp_devs[,4]) ~ as.numeric(prp_devs[,7]), FUN=mean)
devs_bySamp.75 <- aggregate(as.numeric(prp_devs[,5]) ~ as.numeric(prp_devs[,7]), FUN=mean)
devs_bySamp.975 <- aggregate(as.numeric(prp_devs[,6]) ~ as.numeric(prp_devs[,7]), FUN=mean)
devs_bySamp <- data.frame(devs_bySamp.025, devs_bySamp.25[,2], devs_bySamp.5[,2],
                          devs_bySamp.75[,2], devs_bySamp.975[,2])
colnames(devs_bySamp) <- c("n_surveys", "deviance.025", "deviance.25", "deviance.5", 
                           "deviance.75", "deviance.975")
# make plot
thick95 <- 1
#postscript("dev_property_20190411.eps")
pdf(paste0("./fitted_mods/", fit, "/dev_property.pdf"))
plot(log(devs_bySamp$n_surveys), abs(devs_bySamp$deviance.5), pch=16,
     xlab="log number of surveys at property",
     ylab="deviance", axes=FALSE)
axis(1)
axis(2, las=2)
#segments(x0=log(devs_bySamp$n_surveys), y0=devs_bySamp$deviance.025, y1=devs_bySamp$deviance.975,
#         lwd=thick95, col="grey")
dev.off()

#---- examine deviance by time since last survey
survey_d$dev <- NA
survey_d$dev2 <- NA
for(i in 1:nrow(dev_d)){
  survey_d$dev[match(dev_d$idx[i], survey_d$idx)] <- dev_d$dev[i]
  survey_d$dev2[match(dev_d$idx[i], survey_d$idx)] <- dev_d$dev2[i]
}
# In survey_d I need to calculate time since previous survey if it was the same property
surv_sort <- survey_d[order(survey_d$agrp_prp_id, survey_d$start.date), ]
surv_sort$timeSinceSurvey <- NA
for(i in 2:nrow(surv_sort)){ # i can exclude row 1 because it was only survey at this prop
  if(surv_sort$agrp_prp_id[i] == surv_sort$agrp_prp_id[i-1]){
    surv_sort$timeSinceSurvey[i] <- surv_sort$start.date[i] - surv_sort$start.date[i-1]
    # need to use start date for the time difference because there were overlapping surveys,
    # so if we look at the end date of the previous survey, there will be negative times
  }
}
dim(surv_sort[!is.na(surv_sort$timeSinceSurvey),])
min(surv_sort$timeSinceSurvey, na.rm=TRUE) >=0 # (NAs are places surveyed only one time)
# ^ should be TRUE

# I need to get the full posterior for each timeSince Survey insted of using just the median from all of them
# first put the timeSinceSurvey column into the dev dataset using idx
#dev_d$timeSinceSurvey <- surv_sort[match(dev_d$idx, surv_sort$idx), "timeSinceSurvey"]
idx2 <- surv_sort[match(dev_d$idx, surv_sort$idx), c("idx", "timeSinceSurvey")]
ttt <- idx2[match(dev_d$idx, idx2$idx), "timeSinceSurvey"]
write.csv(ttt, "dev_survSort_idx.csv") # need to do this to fix the tidyverse structure imposed on these dfs
ttt <- read.csv("dev_survSort_idx.csv")
write.csv(dev_d, "dev_d2.csv")
dev_d <- read.csv("dev_d2.csv")
dev_d$timeSinceSurvey <- ttt$timeSinceSurvey
time_d <- dev_d[!is.na(dev_d$timeSinceSurvey), ]
time_idx <- split(seq_along(time_d$timeSinceSurvey), time_d$timeSinceSurvey)
days_since <- unique(time_d$timeSinceSurvey)
days_devs <- matrix(NA, nrow=length(days_since), ncol=2) # only taking the median now but could take other credible intervals
days_devs[,1] <- days_since
for(i in days_since){
  post_day <- as.vector(as.matrix(data.frame(dev[, eval(parse(text=paste0("time_idx$`", i, "`")))])))
  days_devs[match(i, days_since), 2] <- median(post_day)
}

# plot the relationship between time since survey and deviance
# first aggregate
#ag_t_dev <- aggregate(surv_sort$dev2 ~ surv_sort$timeSinceSurvey, FUN=mean)
# then plot
#postscript("dev_timeSinceSurvey_20190411.eps")
pdf(paste0("./fitted_mods/", fit, "/dev_timeSinceSurvey.pdf"))
#plot(surv_sort$timeSinceSurvey, surv_sort$dev2, xlab="time since last survey at property (days)",
#     ylab="deviance", axes=FALSE, pch=16)
plot(log(days_devs[,1]), abs(days_devs[,2]), xlab="time since previous survey at property (log of days)",
     ylab="deviance", axes=FALSE, pch=16)
axis(1)
axis(2, las=2)
dev.off()

# look at really bad model fits
bad_fits <- dev_d[abs(dev_d$dev2) > 20, ]
cbind(bad_fits$y, bad_fits$y_hat)

#---------------------------------------------------------------
# compare posterior predicted harvest with 
# observed harvest
#---------------------------------------------------------------
# extract some data
y_rep <- data.frame(rstan::extract(readRDS('./fitted_mods/saturating_fit_20190416.rds'), pars = "y_rep"))
y_rep_dev <- data.frame(rstan::extract(readRDS('./fitted_mods/saturating_fit_20190416.rds'), pars = "y_rep_dev"))

y_obs <- stan_d$y
y_rep_sum <- apply(y_rep, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))
y_rep_sum2 <- data.frame(t(y_rep_sum))
colnames(y_rep_sum2) <- c("lo", "med", "hi")
y_rep_sum2$y_obs <- y_obs
y_rep_sum2$scaled_effort <- stan_d$scaled_effort
y_rep_sum2$effort <- stan_d$effort
y_rep_sum2$effort_per <- stan_d$effort_per
y_rep_sum2$trap_count <- stan_d$trap_count

y_rep_all <- cbind(y_rep, y_rep_dev)
y_rep_all_sum <- data.frame(t(apply(y_rep_all, 2, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
colnames(y_rep_all_sum) <- c("lo", "med", "hi")
y_obs_all <- c(stan_d$y, stan_d$y_dev)

# predicted abunance
lam <- data.frame(rstan::extract(readRDS('fitted_mods/saturating_fit_20190416.rds'), pars = "log_lambda"))
log_p <- data.frame(rstan::extract(readRDS('fitted_mods/saturating_fit_20190416.rds'), pars = "log_p"))
lam_sum <- apply(lam, 2, median)
p_sum <- apply(log_p, 2, median)
lam_sum2 <- data.frame(t(lam_sum), t(p_sum))

# put these together by survey index
# there can be multiple ys for each survey_idx
survey_idx <- stan_d$survey_idx
n <- length(survey_idx)
n_surv <- table(survey_idx)
survey_num <- 1:max(survey_idx)
surv_list <- list() # for every survey conducted, these are the indices in y, y_rep, log_p, that correspond 
for(i in 1:length(survey_num)){
  surv_list[[i]] <- which(survey_idx %in% survey_num[i])
}

# calculate median predicted abundance
u_surv <- unique(survey_idx)
pred_abun <- rep(NA, length(survey_num)) # an abundance for every unique survey
for(i in survey_idx){
  pred_abun[i] <-  exp(lam_sum[i] + median(p_sum[surv_list[[i]]]))
}
pred_abun2 <- pred_abun[!is.na(pred_abun)] 

# calculate median predicted harvest
num_survs <- lapply(surv_list, length)
max(unlist(num_survs))
pred_harv <- matrix(NA, nrow=length(survey_num), ncol=172)
for(i in survey_idx){
  pred_harv[i,1:length(surv_list[[i]])] <- y_rep_sum2$med[surv_list[[i]]]
} 

# get effort for each of these
scaled_effort <- effort <- matrix(NA, nrow=length(survey_num), ncol=172)
for(i in survey_idx){
  scaled_effort[i,1:length(surv_list[[i]])] <- y_rep_sum2$scaled_effort[surv_list[[i]]]
  effort[i,1:length(surv_list[[i]])] <- y_rep_sum2$effort[surv_list[[i]]]
} 

# export
out <- data.frame(pred_abun=pred_abun,
                  pred_harv=pred_harv, 
                  scaled_effort=scaled_effort,
                  effort=effort
)
write_csv(out, "./fitted_mods/fit20190416/post_predicted_harvestCaptureEffort.csv")

write_csv(y_rep_sum2, "./fitted_mods/fit20190416/post_predicted_capture.csv")
write_csv(lam_sum2, "./fitted_mods/fit20190416/post_predicted_lambda.csv")
write_csv(y_rep_all_sum, "./fitted_mods/fit20190416/y_rep_all_sum.csv")

# method
method <- stan_d$method
method_mat <- matrix(NA, nrow=length(survey_num), ncol=172)
for(i in survey_idx){
  method_mat[i,1:length(surv_list[[i]])] <- method[surv_list[[i]]]
} 
write_csv(as.data.frame(method_mat), "./fitted_mods/fit20190416/method_matrix.csv")

#*** loading this on local machine
y_sum <- read_csv("~/Desktop/APHIS/abundanceModeling/wild-pigs/fitted_mods/fit20190416/post_predicted_capture.csv")
hist(y_sum$scaled_effort)

# plot all
plot(y_sum$y_obs, y_sum$med, xlab="observed harvest", ylab="median posterior predicted harvest", axes=FALSE)
axis(1)
axis(2, las=2)

# split by poorly sampled and well sampled
y_poor <- y_sum[y_sum$scaled_effort <= 0,]
y_well <- y_sum[y_sum$scaled_effort > 0, ]

plot(y_well$y_obs, y_well$med, xlab="observed harvest", ylab="median posterior predicted harvest", axes=FALSE,
     col="orange", cex=.2)
axis(1)
axis(2, las=2)
points(y_poor$y_obs, y_poor$med, col="blue", cex=.2)
legend("topleft", legend=c("well-sampled", "poor-sampled"), col=c("orange", "blue"), pch=16, bty="n")

# now deviance
y_sum$dev <- abs(y_sum$y_obs - y_sum$med)
plot(y_sum$scaled_effort,y_sum$dev, ylab="deviance from perfect fit", xlab="scaled effort", axes=FALSE)
axis(1)
axis(2, las=2)

# now compare posterior predicted harvest with capture
harv <- read.csv("~/Desktop/APHIS/abundanceModeling/wild-pigs/fitted_mods/fit20190416/post_predicted_harvestCaptureEffort.csv")
plot(harv$pred_abun, grep("pred_harv", harv))
plot(harv$pred_abun, harv$pred_harv.1)
pred_harv <- harv[, grep("pred_harv", colnames(harv))]
#matplot(harv$pred_abun, pred_harv)

# make plot
plot(harv$pred_abun, pred_harv[,1], xlab="posterior predicted abundance", 
     ylab="posterior predicted harvest", axes=FALSE, cex=.4)
axis(1)
axis(2, las=2)
for(i in 2:ncol(pred_harv)){
  points(harv$pred_abun, pred_harv[,i], cex=.4)
}
abline(a=0,b=1, lty=2, col="grey")

# split by poorly sampled and well sampled
poor <- harv[harv$scaled_effort.1 <= 0,]
well <- harv[harv$scaled_effort.1 > 0, ]

# make plot - well sampled
pred_harv <- well[, grep("pred_harv", colnames(well))]
plot(well$pred_abun, pred_harv[,1], xlab="posterior predicted abundance", 
     ylab="posterior predicted harvest", axes=FALSE, cex=.4, col="dark orange")
axis(1)
axis(2, las=2)
for(i in 2:ncol(pred_harv)){
  points(well$pred_abun, pred_harv[,i], cex=.4, col="dark orange")
}

# add poorly sampled
pred_harv <- poor[, grep("pred_harv", colnames(poor))]
for(i in 1:ncol(pred_harv)){
  points(poor$pred_abun, pred_harv[,i], cex=.4, col="blue")
}

abline(a=0,b=1, lty=2, col="grey")
legend(200, 50, legend=c("well-sampled", "poor-sampled"), col=c("orange", "blue"), pch=1, bty="n")

# get the index of each type of method
methods <- as.matrix(read.csv("~/Desktop/APHIS/abundanceModeling/wild-pigs/fitted_mods/fit20190416/method_matrix.csv"))
c("firearms", "fixed wing", "helicopter", "snare", "trap")
pred_harv <- harv[, grep("pred_harv", colnames(harv))]
meth_used <- as.integer(1:5)
ones1 <- ifelse(methods==1, 1, NA)
sum(ones1, na.rm=TRUE)
ones <- matrix(NA, nrow=length(survey_idx), ncol=172) # these are the indices where methods ==1
twos1 <- ifelse(methods==2, 2, NA)
twos <- matrix(NA, nrow=length(survey_idx), ncol=172) 
threes1 <- ifelse(methods==3, 3, NA)
threes <- matrix(NA, nrow=length(survey_idx), ncol=172) 
fours1 <- ifelse(methods==4, 4, NA)
fours <- matrix(NA, nrow=length(survey_idx), ncol=172) 
fives1 <- ifelse(methods==5, 5, NA)
fives <- matrix(NA, nrow=length(survey_idx), ncol=172) 
for(i in survey_idx){
  # ones
  ones2 <- which(ones1[i,] == 1)
  ones3 <- tmp[!is.na(tmp)]
  if(length(ones3>0)){
    ones[i, 1:length(tmp)] <- tmp
  }
  # twos
  twos2 <- which(twos1[i,] == 2)
  twos3 <- tmp[!is.na(tmp)]
  if(length(twos3>0)){
    twos[i, 1:length(tmp)] <- tmp
  }
  # threes
  threes2 <- which(threes1[i,] == 3)
  threes3 <- tmp[!is.na(tmp)]
  if(length(threes3>0)){
    threes[i, 1:length(tmp)] <- tmp
  }
  # fours
  fours2 <- which(fours1[i,] == 4)
  fours3 <- tmp[!is.na(tmp)]
  if(length(fours3>0)){
    fours[i, 1:length(tmp)] <- tmp
  }
  # fives
  fives2 <- which(fives1[i,] == 5)
  fives3 <- tmp[!is.na(tmp)]
  if(length(fives3>0)){
    fives[i, 1:length(tmp)] <- tmp
  }
} # I didn't need to do this because method is in the same index format as dev

# compare deviance across methods
pdf("./fitted_mods/fit20190416/absVal_deviance_byMethod_20190906.pdf")
par(mar=c(5,8,4,3))
boxplot(log(y_sum$dev + 1) ~ stan_d$method, axes=FALSE, 
        xlab="Method of capture",
        ylab="Absolute value of deviance of \nabundance prediction from observed")
axis(1, at=1:5, labels=c("firearms", "fixed wing", "helicopter", "snare", "trap"))
#axis(2, las=2, at=0:5, labels=exp(0:5))
axis(2, las=2, at=log(c(1, 3, 8, 21, 56, 151)), labels=c(0, 2, 7, 20, 55, 150))
dev.off()

# do this with just value (not absolute value)
y_sum$devVal <- y_sum$med -y_sum$y_obs
y_sum2 <- y_sum[y_sum$devVal > -280,]
which(y_sum)
#pdf("./fitted_mods/fit20190416/val_deviance_byMethod_20190906.pdf")
par(mar=c(5,8,4,3))
bp <- boxplot(y_sum$devVal ~ stan_d$method, axes=FALSE, 
              #bp <- boxplot(log(y_sum$devVal[-which.min(y_sum$devVal)]+185) ~ stan_d$method[-which.min(y_sum$devVal)], axes=FALSE, 
              xlab="Method of capture",
              ylab="deviance", range=1000)
axis(1, at=1:5, labels=c("firearms", "fixed wing", "helicopter", "snare", "trap"))
#axis(2, las=2, at=0:5, labels=exp(0:5))
axis(2, las=2)
dev.off()
# firearms and helicopter are overestimating the worst

#---------------------------------------------------------------
# look at the effect of effort on abundance of pigs
#---------------------------------------------------------------
# read in density and abundance data
#y_sum <- read_csv("./fitted_mods/fit20190416/post_predicted_capture.csv")
y_sum <- read_csv("./fitted_mods/fit20190416/y_rep_all_sum.csv") # this one has all of the datapoints, but need to add some stuff from stan
y_sum$y_obs <- c(stan_d$y, stan_d$y_dev)
y_sum$scaled_effort <- c(stan_d$scaled_effort, stan_d$scaled_effort_dev)
y_sum$effort <- c(stan_d$effort, stan_d$effort_dev)
y_sum$property_idx <- c(stan_d$p_property_idx, stan_d$p_property_idx_dev)
y_sum$county_idx <- c(stan_d$p_county_idx, stan_d$p_county_idx_dev)
#st_d$county_index
scaled_effort <- stan_d$scaled_effort
effort <- stan_d$effort
survey_idx <- c(stan_d$survey_idx, stan_d$survey_idx_dev)
# add spatio temporal data to y_sum
y_sum$FIPS <- st_d$FIPS[survey_idx]
y_sum$timestep <- st_d$timestep[survey_idx]
y_sum$agrp_prp_id <- st_d$agrp_prp_id[survey_idx]
y_sum$method <- stan_d$method[survey_idx]
y_sum$state <- st_d$STATE_NAME[survey_idx]
y_sum$property_area_km2 <- st_d$property_area_km2[survey_idx]
y_sum$med_dens <- y_sum$med/y_sum$property_area_km2
y_sum$lo_dens <- y_sum$lo/y_sum$property_area_km2
y_sum$hi_dens <- y_sum$hi/y_sum$property_area_km2
y_sum2 <- merge(x=y_sum, y=timestep_df)
y_sum <- y_sum2
# deal with dates
y_sum$year <- year(y_sum$end_dates)
mo <- month(y_sum$end_dates)
mon <- ifelse(nchar(mo)==2, mo, paste0(0, mo))
y_sum$ym <- paste0(y_sum$year, mon)

# read in density from abundance prediction
dens <- read_csv("./fitted_mods/fit20190416/dens_summary_cnty_20190910.csv")
dens$FIPS <- dens$fips
dens$timestep <- as.integer(dens$timestep)
dens2 <- merge(dens, y_sum, by=c("FIPS", "timestep"))
anyNA(dens2)

# add in modifier of acceptable area
{
  scwds <- read_csv("./data/SCWDS_CNTY_All_Years_1960-2018_FINAL_3Apr2019.csv")
  scwds$fips <- ifelse(nchar(scwds$FIPS) == 5, scwds$FIPS, paste0(0, scwds$FIPS))
  # add years where no scwds data
  # I'm assuming nothing changed in proportion of county occupied between 2004-2007
  a2005 <- scwds[scwds$YEAR==2004,]
  a2005$YEAR <- 2005
  a2006 <- scwds[scwds$YEAR==2004,]
  a2006$YEAR <- 2006
  a2007 <- scwds[scwds$YEAR==2004,]
  a2007$YEAR <- 2007
  scwds2 <- rbind(scwds, a2005, a2006, a2007)
  scwds <- scwds2
  # remove exclude flags
  scwds2 <- scwds[scwds$ExcludeFlag == 0, ]
  prp_cnty <- aggregate(scwds2$PropCNTY_Pig, by=list(scwds2$fips,scwds2$YEAR), FUN=mean)
  colnames(prp_cnty) <- c("FIPS", "year", "prp_cnty_occupied")
}
dens3 <- merge(dens2, prp_cnty)
dens3$med_dens <- dens3$med.x * dens3$prp_cnty_occupied
#write.csv(prp_cnty, "proporition_county_occupied.csv")
# make a function that does this
calc_slopes <- function(max_months, 
                        min_months=1, 
                        trueAbun=TRUE, # using true predicted abundance instead of y_sum
                        compare_meth="slope", 
                        sp_unit = "property",
                        slope_calc="effort", 
                        remove_zeros = FALSE
){
  
  cap_methods <- c("firearms", "fixed wing", "helicopter", "snare", "trap")
  
  
  if(remove_zeros){
    y_sum <- y_sum[y_sum$med_dens != 0, ]
  }
  
  # set up for using different spatial units
  if(sp_unit=="property"){
    props <- unique(y_sum$agrp_prp_id)
    seqer <- seq_along(props)
    reduce_df <- matrix(NA, nrow=length(props), ncol=6)
  } 
  if(sp_unit=="county"){
    cntys <- unique(y_sum$FIPS)
    seqer <- seq_along(cntys)
    reduce_df <- matrix(NA, nrow=length(cntys), ncol=6)
  } 
  if(sp_unit=="state"){
    states <- unique(y_sum$state)
    seqer <- seq_along(states)
    reduce_df <- matrix(NA, nrow=length(states), ncol=6)
  }
  
  for(i in seqer){
    # subset data based on spatial unit.
    # I'm calling them all prop1, regardless of the spatial unit for simplicity below
    if(sp_unit=="property"){
      prop1 <- y_sum[y_sum$agrp_prp_id == props[i],]
    } 
    if(sp_unit=="county"){
      prop1 <- y_sum[y_sum$FIPS == cntys[i],]
    } 
    if(sp_unit=="state"){
      prop1 <- y_sum[y_sum$state== states[i],]
    }
    
    
    # if work was only conducted at one timestep, ignore this property
    if(length(unique(prop1$timestep)) < 2){ 
      reduce_df[i, ] <- rep(NA, ncol(reduce_df))
    } else{
      # calculate the time differences from one timestep to the next
      diffs <- ave(prop1$timestep, FUN=function(x) c(NA,diff(x)))
      # get the indeces for time differences that are of interest
      short_idxs <- which(diffs <= max_months & diffs >= min_months)
      
      # if there are not 2 instances where the time difference between removal is approporiate, ignore this st unit
      if(length(short_idxs) < 2 ){ 
        reduce_df[i, ] <- rep(NA, ncol(reduce_df))
      } else{
        # I cannot calculate the difference if it is the last year of work being done at the st unit
        if(max(short_idxs) == nrow(prop1)){ 
          short_idxs <- short_idxs[-(max(short_idxs))]
        } else{
          # if there are not 2 instances where the time difference between removal is approporiate, ignore this st unit
          if(length(short_idxs) < 2) { 
            reduce_df[i, ] <- rep(NA, ncol(reduce_df))
          } else{
            prop1$diff <- rep(NA, nrow(prop1))
            # the equation for calculating difference
            difference <- (prop1$med_dens[short_idxs]-prop1$med_dens[short_idxs+1])/(prop1$med_dens[short_idxs]+1) # relative difference
            #difference <- (prop1$med_dens[short_idxs]-prop1$med_dens[short_idxs+1]) #absolute difference
            #difference <- (prop1$hi_dens[short_idxs]-prop1$hi_dens[short_idxs+1])
            prop1$diff[short_idxs] <- difference
            
            # compare difference to effort by slope or correlation
            if(compare_meth == "slope"){
              if(slope_calc == "effort"){
                # slope for all methods combined
                #lm1 <- lm(prop1$diff ~ prop1$scaled_effort)
                #lm1 <- lm(prop1$med_dens ~ prop1$scaled_effort)
                lm1 <- lm(prop1$diff ~ prop1$effort, na.action="na.exclude")
                #lm1 <- lm(prop1$y_obs ~ prop1$effort)
              }
              if(slope_calc == "num_events"){
                prop2 <- prop1[short_idxs, ]
                ag1 <- aggregate(prop2$method ~ prop2$ym, FUN=length)
                ag2 <- aggregate(prop2$med_dens ~ prop2$ym, FUN=median)
                #ag2 <- aggregate(prop2$med ~ prop2$ym, FUN=median)
                lm1 <- lm(ag2[,2] ~ ag1[,2])
              }
              
              reduce_df[i,1] <- lm1$coefficients[2]
            }
            if(compare_meth == "cor"){
              reduce_df[i,1] <- cor(prop1$scaled_effort, prop1$diff, use="complete.obs",
                                    method="spearman")
            }
            
            # now look at the slope for each method
            for(j in seq_along(cap_methods)){
              prop1_meth <- prop1[prop1$method==j, ]
              prop1_meth <- prop1_meth[!is.na(prop1_meth$diff), ]
              if(nrow(prop1_meth) < 2){
                reduce_df[i,(j+1)] <- NA
              } else{
                if(compare_meth == "slope"){
                  #lm2 <- lm(prop1_meth$diff ~ prop1_meth$scaled_effort)
                  #lm2 <- lm(prop1_meth$med_dens ~ prop1_meth$scaled_effort)
                  lm2 <- lm(prop1_meth$diff ~ prop1_meth$effort, na.action="na.exclude")
                  #lm2 <- lm(prop1_meth$y_obs ~ prop1_meth$effort)
                  reduce_df[i,(j+1)] <- lm2$coefficients[2]
                }
                if(compare_meth == "cor"){
                  reduce_df[i,(j+1)] <- cor(prop1_meth$scaled_effort, prop1_meth$diff, use="complete.obs",
                                            method="spearman")
                }
              }
            } 
          }
        }
        
      }
      
    }
  }
  colnames(reduce_df) <- c("all_methods", "firearms", "fixed_wing", "helicopter", "snare", "trap")
  assign(paste0("quant_", min_months, "_",max_months),  apply(reduce_df, 2, function(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975), na.rm=TRUE)), envir = .GlobalEnv)
  print(paste0("number of slopes/correlations calculated = ", sum(!is.na(reduce_df))))
}

# make a plot of slope by method

# now run this for several versions of months after removal effort at the property level
calc_slopes(1,1,"slope")
calc_slopes(3,1,"slope")
calc_slopes(6,3,"slope")
calc_slopes(9,6,"slope")
calc_slopes(12, 9, "slope")
# using number of events as a measure of effort
calc_slopes(1,1,"slope", sp_unit = "property", slope_calc = "num_events")
calc_slopes(3,1,"slope", sp_unit = "property", slope_calc = "num_events")
calc_slopes(6,3,"slope", sp_unit = "property", slope_calc = "num_events")
calc_slopes(9,6,"slope", sp_unit = "property", slope_calc = "num_events")
calc_slopes(12, 9, "slope", sp_unit = "property", slope_calc = "num_events" )

#- at the county level
#*** this next line makes it so that we are doing density from the abundance estimates
y_sum <- dens3
calc_slopes(1,1,"slope", sp_unit = "county")
calc_slopes(3,1,"slope", sp_unit = "county")
calc_slopes(6,3,"slope", sp_unit = "county")
calc_slopes(9,6,"slope", sp_unit = "county")
calc_slopes(12, 9, "slope", sp_unit = "county")
# using number of events as a measure of effort
calc_slopes(1,1,"slope", sp_unit = "county", slope_calc = "num_events", remove_zeros = TRUE)
calc_slopes(3,1,"slope", sp_unit = "county", slope_calc = "num_events", remove_zeros = TRUE)
calc_slopes(6,3,"slope", sp_unit = "county", slope_calc = "num_events", remove_zeros = TRUE)
calc_slopes(9,6,"slope", sp_unit = "county", slope_calc = "num_events", remove_zeros = TRUE)
calc_slopes(12, 9, "slope", sp_unit = "county", slope_calc = "num_events" , remove_zeros = TRUE)

# at the state level
calc_slopes(1,1,"slope", sp_unit = "state", slope_calc = "num_events", remove_zeros = TRUE)
calc_slopes(3,1,"slope", sp_unit = "state", slope_calc = "num_events", remove_zeros = TRUE)
calc_slopes(6,3,"slope", sp_unit = "state", slope_calc = "num_events", remove_zeros = TRUE)
calc_slopes(9,6,"slope", sp_unit = "state", slope_calc = "num_events", remove_zeros = TRUE)
calc_slopes(12, 9, "slope", sp_unit = "state", slope_calc = "num_events", remove_zeros = TRUE)

# plotting parameters
months_used <- c(1,3,6,9,12)
var <- 1/(length(months_used)+1)
j=1 # for the first one
lwd50 <- 3
lwd95 <- 1
cex.pch <- 1

# make plot
#pdf("./fitted_mods/fit20190416/method_reduce_density_over_time2.pdf")
par(mar=c(5,8,4,3))
plot(0, 0,
     #xlim=c(min(eval(parse(text=paste0("quant_1_", months_used[j])))), max(eval(parse(text=paste0("quant_1_", months_used[j]))))),
     xlim=c(-3,2.5),
     #ylim=c(1, 7), #xlim=c(-2,1),
     ylim=c(1, 2),
     axes=FALSE, 
     xlab="", 
     ylab="", cex=cex.pch, pch=16)
abline(v=0, lty=2, col="grey", lwd=1.5)
for(j in seq_along(months_used)){
  # points(eval(parse(text=paste0("quant_1_", months_used[j], "[3,", j, "]"))), j, cex=cex.pch, pch=16)
  
  for(i in 1:6){
    segments(x0=eval(parse(text=paste0("quant_1_", months_used[j], "[2,", i, "]"))),
             x1=eval(parse(text=paste0("quant_1_", months_used[j], "[4,", i, "]"))),
             y0=(i + (j)*var), lwd=lwd50, pch=j, col=j)
    segments(x0=eval(parse(text=paste0("quant_1_", months_used[j], "[1,", i, "]"))),
             x1=eval(parse(text=paste0("quant_1_", months_used[j], "[5,", i, "]"))),
             y0=(i + (j)*var), lwd=lwd95, pch=j, col=j)
    points(eval(parse(text=paste0("quant_1_", months_used[j], "[3,", i, "]"))), 
           (i + (j)*var), cex=cex.pch, pch=j, col=j)
    
  }
}
axis(1)
axis(2, at=(1:6)+.5, tick=FALSE,
     labels=c("all methods", "firearms", "fixed wing", "helicopter", "snare", "trap"), las=2)

mtext(paste0("reduction of wild pig density at a property\n by next removal effort within x", " months"), 1,line=3)
#mtext(paste0("correlation of wild pig density at a property\n and next removal effort within x", " months"), 1,line=3)
legend("topleft", pch=rev(1:length(months_used)), col=rev(1:length(months_used)),
       legend=paste0(rev(months_used), " months since removal"), bty="n")
dev.off()

#---------------------------------------------------------------
# evaluate the effect of effort on abundance by 
# looking at residuals from slope of a moving window
#---------------------------------------------------------------
# read in density and abundance data
#y_sum <- read_csv("./fitted_mods/fit20190416/post_predicted_capture.csv")
y_sum <- read_csv("./fitted_mods/fit20190416/y_rep_all_sum.csv") # this one has all of the datapoints, but need to add some stuff from stan
y_sum$y_obs <- c(stan_d$y, stan_d$y_dev)
y_sum$scaled_effort <- c(stan_d$scaled_effort, stan_d$scaled_effort_dev)
y_sum$effort <- c(stan_d$effort, stan_d$effort_dev)
y_sum$property_idx <- c(stan_d$p_property_idx, stan_d$p_property_idx_dev)
y_sum$county_idx <- c(stan_d$p_county_idx, stan_d$p_county_idx_dev)
#st_d$county_index
scaled_effort <- stan_d$scaled_effort
effort <- stan_d$effort
survey_idx <- c(stan_d$survey_idx, stan_d$survey_idx_dev)
# add spatio temporal data to y_sum
y_sum$FIPS <- st_d$FIPS[survey_idx]
y_sum$timestep <- st_d$timestep[survey_idx]
y_sum$agrp_prp_id <- st_d$agrp_prp_id[survey_idx]
y_sum$method <- stan_d$method[survey_idx]
y_sum$state <- st_d$STATE_NAME[survey_idx]
y_sum$property_area_km2 <- st_d$property_area_km2[survey_idx]
#y_sum$med_dens <- y_sum$med/y_sum$property_area_km2
#y_sum$lo_dens <- y_sum$lo/y_sum$property_area_km2
#y_sum$hi_dens <- y_sum$hi/y_sum$property_area_km2
y_sum2 <- merge(x=y_sum, y=timestep_df)
y_sum <- y_sum2
# deal with dates
y_sum$year <- year(y_sum$end_dates)
mo <- month(y_sum$end_dates)
mon <- ifelse(nchar(mo)==2, mo, paste0(0, mo))
y_sum$ym <- paste0(y_sum$year, mon)

# deal with prp_cnty occupied
prp_cnty2 <- prp_cnty[prp_cnty$year > 1998,]

# read in density from abundance prediction
dens <- data.frame(read_csv("./fitted_mods/fit20190416/dens_summary_cnty_20190912.csv"))
dens$year <- year(dens$end_dates)
dens$FIPS <- dens$fips
dens$timestep <- as.integer(dens$timestep)
#dens2 <- merge(dens, prp_cnty2, by=c("FIPS", "year"), all.x=TRUE, all.y=FALSE)
# colnames(dens2)[colSums(is.na(dens2)) > 0]
# tmp <- dens2[is.na(dens2$prp_cnty_occupied),]
# unique(tmp$FIPS)[which(unique(tmp$FIPS) %in% unique(prp_cnty2$FIPS))]
# unique(tmp$FIPS)[6]
# prp_cnty2[prp_cnty2$FIPS=="04001",] # proportion of counties occupied is missing data:: stop using this!
# dens2 <- dens2[!is.na(dens2$med),]
# dens2$med_dens <- dens2$med*dens2$prp_cnty_occupied
dens$FIPS2 <- as.integer(dens$FIPS)
y_sum$FIPS2 <- as.integer(y_sum$FIPS)
dens3 <- merge(dens, y_sum, by=c("timestep", "FIPS2"), all.x=TRUE, all.y=FALSE)
u_fips <- unique(dens3$FIPS.x)
# make a table of just states for adding
all_state <- dens3[!duplicated(dens3$FIPS.x),c(3,9)]


# get the slope of each county in a five year window
# then calculate the residuals for each data point where there 
# is a method (these are the times when culling happened)
#u_fips <- unique(y_sum$FIPS)
t_year <- 2 # = number of timesteps: 66 in 5 years. I might want to make this shorter
timesteps <- 1:174
n_windows <- length(timesteps)-t_year # number of moving windows 
out <- out1 <- out2 <- out3 <- out4 <- out5 <- matrix(NA, nrow=length(u_fips), ncol=n_windows)
#out <- array(NA, dim=c(length(u_fips), (n_windows+1),5))

# remove values we know are too high
#dens3 <- dens3[dens3$med.x < 100, ]
# see if there is an effort value for every method
#methT <- ifelse((!is.na(dens3$method)), 1, 0)
#effT <- ifelse((!is.na(dens3$effort)), 1, 0)
#all.equal(methT, effT)

# do this for all states
# separate it out by method
u_fips <- unique(dens3$FIPS.x)
out_mn <- out_md <- out_direction <-  matrix(NA, nrow=length(u_fips), ncol=n_windows)
out1 <- out2 <- out3 <- out4 <- out5 <- matrix(NA, nrow=length(u_fips), ncol=n_windows)
cor1 <- cor2 <- cor3 <- cor4 <- cor5 <- matrix(NA, nrow=length(u_fips), ncol=n_windows)
dif1 <- dif2 <- dif3 <- dif4 <- dif5 <- matrix(NA, nrow=length(u_fips), ncol=n_windows)
slp <- slp1 <- slp2 <- slp3 <- slp4 <- slp5 <- matrix(NA, nrow=length(u_fips), ncol=n_windows)


for(i in seq_along(u_fips)){
  cntyD <- dens3[dens3$FIPS.x == u_fips[i], ]
  #cntyD <- sub_method[sub_method$FIPS.x == u_fips[i], ]
  for(j in 1:n_windows){
    t_use <- timesteps[j:(j+t_year)]
    cntyDT2 <- cntyD[which(cntyD$timestep %in% t_use),]
    
    # calculate residuals
    if(sum(!is.na(cntyDT2$med.x))>1){
      if(sum(is.na(cntyDT2$method)) != nrow(cntyDT2)){
        lm1 <- lm(cntyDT2$med.x ~ cntyDT2$timestep)
        work_idxs <- which(!is.na(cntyDT2$method)) 
        #direction <- ifelse(lm1$residuals>0, 1, -1)
        #out_mn[i, j] <- mean(lm1$residuals[(work_idxs+1)], na.rm=TRUE) # na.rm to deal with situations when work was done in the last timestep
        #out_md[i, j] <- median(lm1$residuals[work_idxs+1], na.rm=TRUE)
        # out_direction[i, j] <- mean(direction[work_idxs+1])
        
        # separate the work indices by method
        work1 <- which(cntyDT2$method == 1)
        work2 <- which(cntyDT2$method == 2)
        work3 <- which(cntyDT2$method == 3)
        work4 <- which(cntyDT2$method == 4)
        work5 <- which(cntyDT2$method == 5)
        # then put these into out tables
        out1[i,j] <- median(lm1$residuals[work1+1], na.rm=TRUE)
        out2[i,j] <- median(lm1$residuals[work2+1], na.rm=TRUE)
        out3[i,j] <- median(lm1$residuals[work3+1], na.rm=TRUE)
        out4[i,j] <- median(lm1$residuals[work4+1], na.rm=TRUE)
        out5[i,j] <- median(lm1$residuals[work5+1], na.rm=TRUE)
        # removing the plus 1 from the index so I am just looking at the timestep
        # where removal occurred instead of the timestep after removal occurred
        out1[i,j] <- median(lm1$residuals[work1], na.rm=TRUE)
        out2[i,j] <- median(lm1$residuals[work2], na.rm=TRUE)
        out3[i,j] <- median(lm1$residuals[work3], na.rm=TRUE)
        out4[i,j] <- median(lm1$residuals[work4], na.rm=TRUE)
        out5[i,j] <- median(lm1$residuals[work5], na.rm=TRUE)
        
        # just collect the slope from each county in each window
        slp[i,j] <- lm1$coefficients[1]
        
        
        # only calculate the difference from previous step to work step
        dif1[i, j] <- median(cntyDT2$med.x[(work1-1)] - cntyDT2$med.x[work1], na.rm=TRUE)
        dif2[i, j] <- median(cntyDT2$med.x[(work2-1)] - cntyDT2$med.x[work2], na.rm=TRUE)
        dif3[i, j] <- median( cntyDT2$med.x[(work3-1)] - cntyDT2$med.x[work3], na.rm=TRUE)
        dif4[i, j] <- median(cntyDT2$med.x[(work4-1)] - cntyDT2$med.x[work4] , na.rm=TRUE)
        dif5[i, j] <- median(cntyDT2$med.x[(work5-1)] -cntyDT2$med.x[work5] , na.rm=TRUE)
        # subtract timestep with work from timestep after work was done
        # dif1[i, j] <- median(cntyDT2$med.x[work1] - cntyDT2$med.x[(work1+1)], na.rm=TRUE)
        # dif2[i, j] <- median(cntyDT2$med.x[work2] - cntyDT2$med.x[(work2+1)], na.rm=TRUE)
        # dif3[i, j] <- median(cntyDT2$med.x[work3] - cntyDT2$med.x[(work3+1)], na.rm=TRUE)
        # dif4[i, j] <- median(cntyDT2$med.x[work4] - cntyDT2$med.x[(work4+1)], na.rm=TRUE)
        # dif5[i, j] <- median(cntyDT2$med.x[work5] - cntyDT2$med.x[(work5+1)], na.rm=TRUE)
        # calculate the correlation between residuals and effort
        # if(length(work1) > 1){
        #   cor1[i,j] <- cor(lm1$residuals[work1+1], cntyDT2$effort[work1], use="complete.obs")
        # }
        # if(length(work2) >1){
        #   cor2[i,j] <- cor(lm1$residuals[work2+1], cntyDT2$effort[work2], use="complete.obs")
        # }
        # if(length(work3) >1){
        # cor3[i,j] <- cor(lm1$residuals[work3+1], cntyDT2$effort[work3], use="complete.obs")
        # }
        # if(length(work4) >1 ){
        # cor4[i,j] <- cor(lm1$residuals[work4+1], cntyDT2$effort[work4], use="complete.obs")
        # }
        # if(length(work5) >1){
        # cor5[i,j] <- cor(lm1$residuals[work5+1], cntyDT2$effort[work5], use="complete.obs")
        # }
      }
    }
  }
  if(i %% 100 == 0){
    print(paste0(i," of ", length(u_fips), " counties complete"))
  }
} # close i


# post process function
cap_methods <- c("firearms", "fixed wing", "helicopter", "snare", "trap")
ag_out <- function(
  out_md = out1, method=1, median_method="state"
){
  out_md2 <- data.frame(out_md, u_fips)
  out_md3 <- merge(out_md2, all_state, by.x="u_fips", by.y="fips")
  states <- unique(out_md3$state.x)
  states <- states[!is.na(states)]
  if(median_method=="state"){ # taking the median for each state
    out_md_sum <- matrix(NA, nrow=length(states), ncol=3)
    for(i in seq_along(states)){
      state_dat <- out_md3[out_md3$state.x == states[i],]
      state_dat2 <- unlist(state_dat[,2:(ncol(state_dat)-1)])
      out_md_sum[i, ] <- quantile(state_dat2, c(0.025, 0.5, 0.975), na.rm=TRUE)
    }
    out_md_sum2 <- data.frame(state=states, out_md_sum, method=cap_methods[method])
  }
  else{
    #ag_on <- rowMeans(out_md3[,2:(ncol(out_md3)-1)], na.rm=TRUE)
    ag_on <- apply(out_md3, 1, function(x) median(x,na.rm=TRUE))
    out_md_sum <- aggregate(ag_on ~ out_md3$state.x, 
                            FUN=function(x) quantile(x, c(0.025, 0.5, 0.975)))
    colnames(out_md_sum) <- c("state", "residual")
    out_md_sum2 <- data.frame(out_md_sum[,1], out_md_sum[,2], method=cap_methods[method])
    
  }
  colnames(out_md_sum2) <- c("state", "lo", "med", "hi", "method")
  out_md_sum2 <- out_md_sum2[!(is.na(out_md_sum2$med)),]
  #hist(out_md_sum2[,3], breaks=100, main=cap_methods[method])
  hist(as.numeric(unlist(out_md3[,2:(ncol(out_md3)-1)])), breaks=100, main=cap_methods[method])
  return(list(
    full=out_md3,
    summary=out_md_sum2)
  )
  cat(paste0(colMeans(out_md_sum2[2:4])))
}
out1_sum <- ag_out(out1, method=1)
out2_sum <- ag_out(out2,2)
out3_sum <- ag_out(out3, 3)
out4_sum <- ag_out(out4, 4)
out5_sum <- ag_out(out5, 5)

out_66timesteps <- rbind(out1_sum, out2_sum, out3_sum, out4_sum, out5_sum)
#write_csv(out_55timesteps, "./fitted_mods/fit20190416/residuals_55timesteps_20190912.csv")

cor1_sum <- ag_out(cor1, 1)
cor2_sum <- ag_out(cor2, 2)
cor3_sum <- ag_out(cor3, 3)
cor4_sum <- ag_out(cor4, 4)
cor5_sum <- ag_out(cor5, 5)

# look at difference between work timestep and previous timestep
dif1_sum <- ag_out(dif1, 1)
dif2_sum <- ag_out(dif2, 2)
dif3_sum <- ag_out(dif3, 3)
dif4_sum <- ag_out(dif4, 4)
dif5_sum <- ag_out(dif5, 5)

# make plot of residual over time for each county
{
  out1v <- out1_sum$full[rowSums(is.na(out1_sum$full)) != (ncol(out1_sum$full)-2),]
  out1v <- data.frame(out1v[,1], scale(out1v[,2:(ncol(out1v)-1)]), out1v[, ncol(out1v)])
  max(out1v[,2:(ncol(out1v)-1)], na.rm=TRUE)
  min(out1v[,2:(ncol(out1v)-1)], na.rm=TRUE)
  plot(t_year+(1:(ncol(out1v)-2)), out1v[1,2:(ncol(out1v)-1)], type="l", ylim=c(-10,10), axes=FALSE,
       xlab="Timestep at end of moving window", ylab="Residuals at timestep after work was completed")
  abline(h=0, lty=2, col="grey")
  for(i in 1:nrow(out1v)){
    lines(t_year+(1:(ncol(out1v)-2)), out1v[i,2:(ncol(out1v)-1)], lty=1)
  }
  axis(1)
  axis(2, las=2)
  mtext(paste0("Capture method = ", cap_methods[1]), line=1)
  
  # break out this method across all states -- for method 1
  states1 <- unique(out1v[,ncol(out1v)])
  mat1 <- matrix(1:30, nrow=6, ncol=5)
  
  pdf(paste0("./fitted_mods/fit20190416/residual_removal_plots/firearm.pdf"))
  layout(mat1)
  par(mar=c(3,2,0.1,0.1), oma=c(1,3.5,2,.6))
  for(j in seq_along(states1)){
    stD <- out1v[out1v[,ncol(out1v)] == states1[j],]
    up <- max(stD[, 2:(ncol(stD)-1)], na.rm=TRUE)
    down <- min(stD[, 2:(ncol(stD)-1)], na.rm=TRUE)
    plot(t_year+(1:(ncol(stD)-2)), stD[1,2:(ncol(stD)-1)], type="l", axes=FALSE,
         ylim=c(down,up),
         xlab="", ylab="")
    abline(h=0, lty=2, col="grey")
    for(i in 1:nrow(stD)){
      lines(t_year+(1:(ncol(stD)-2)), stD[i,2:(ncol(stD)-1)])
    }
    axis(2, las=2)
    mtext(states1[j])
    if(j %% nrow(mat1) ==0 | j == length(states1)){
      axis(1)
    } 
    else {
      axis(1, label=FALSE)
    }
    
  }
  mtext("Residuals at timestep after work was completed",
        side=2, outer=TRUE, 
        cex=1.2, line=1)
  mtext("Timestep at end of moving window", side=1, outer=TRUE, line=-.4, cex=1.2)
  dev.off()
  
  
  # method 2
  out2v <- out2_sum$full[rowSums(is.na(out2_sum$full)) != (ncol(out2_sum$full)-2),]
  out2v <- data.frame(out2v[,1], scale(out2v[,2:(ncol(out2v)-1)]), out2v[, ncol(out2v)])
  max(out2v[,2:(ncol(out2v)-1)], na.rm=TRUE)
  min(out2v[,2:(ncol(out2v)-1)], na.rm=TRUE)
  plot(t_year+(1:(ncol(out2v)-2)), out2v[1,2:(ncol(out2v)-1)], type="l", ylim=c(-10,10), axes=FALSE,
       xlab="Timestep at end of moving window", ylab="Residuals at timestep after work was completed")
  abline(h=0, lty=2, col="grey")
  for(i in 1:nrow(out2v)){
    lines(t_year+(1:(ncol(out2v)-2)), out2v[i,2:(ncol(out2v)-1)], lty=1)
  }
  axis(1)
  axis(2, las=2)
  mtext(paste0("Capture method = ", cap_methods[2]), line=1)
  
  
  # break out this method across all states -- for method 2
  states2 <- unique(out2v[,ncol(out2v)])
  mat2 <- matrix(1:6, nrow=3, ncol=2)
  
  pdf(paste0("./fitted_mods/fit20190416/residual_removal_plots/fixed_wing.pdf"))
  layout(mat2)
  par(mar=c(3,2,0.1,0.1), oma=c(1,3.5,2,.6))
  for(j in seq_along(states2)){
    stD <- out2v[out2v[,ncol(out2v)] == states2[j],]
    up <- max(stD[, 2:(ncol(stD)-1)], na.rm=TRUE)
    down <- min(stD[, 2:(ncol(stD)-1)], na.rm=TRUE)
    plot(t_year+(1:(ncol(stD)-2)), stD[1,2:(ncol(stD)-1)], type="l", axes=FALSE,
         ylim=c(down,up),
         xlab="", ylab="")
    abline(h=0, lty=2, col="grey")
    for(i in 1:nrow(stD)){
      lines(t_year+(1:(ncol(stD)-2)), stD[i,2:(ncol(stD)-1)])
    }
    axis(2, las=2)
    mtext(states2[j])
    if(j %% nrow(mat2) ==0 | j == length(states2)){
      axis(1)
    } 
    else {
      axis(1, label=FALSE)
    }
    
  }
  mtext("Residuals at timestep after work was completed",
        side=2, outer=TRUE, 
        cex=1.2, line=1)
  mtext("Timestep at end of moving window", side=1, outer=TRUE, line=-.4, cex=1.2)
  dev.off()
  
  # make plot of residual over time for each county
  out3v <- out3_sum$full[rowSums(is.na(out3_sum$full)) != (ncol(out3_sum$full)-2),]
  #out3v <- out3_sum$full[rowSums(is.na(out3_sum$full)) < 50,]
  out3v <- data.frame(out3v[,1], scale(out3v[,2:(ncol(out3v)-1)]), out3v[, ncol(out3v)])
  max(out3v[,2:(ncol(out3v)-1)], na.rm=TRUE)
  min(out3v[,2:(ncol(out3v)-1)], na.rm=TRUE)
  plot(t_year+(1:(ncol(out3v)-2)), out3v[1,2:(ncol(out3v)-1)], type="l", ylim=c(-10,10), axes=FALSE,
       xlab="Timestep at end of moving window", ylab="Residuals at timestep after work was completed")
  abline(h=0, lty=2, col="grey")
  for(i in 1:nrow(out3v)){
    lines(t_year+(1:(ncol(out3v)-2)), out3v[i,2:(ncol(out3v)-1)], lty=1)
  }
  axis(1)
  axis(2, las=2)
  mtext(paste0("Capture method = ", cap_methods[3]), line=1)
  
  # break out this method across all states -- for method 3
  states3 <- unique(out3v[,ncol(out3v)])
  mat3 <- matrix(1:8, nrow=4, ncol=2)
  
  pdf(paste0("./fitted_mods/fit20190416/residual_removal_plots/helicopter.pdf"))
  layout(mat3)
  par(mar=c(3,2,0.1,0.1), oma=c(1,3.5,2,.6))
  for(j in seq_along(states3)){
    stD <- out3v[out3v[,ncol(out3v)] == states3[j],]
    up <- max(stD[, 2:(ncol(stD)-1)], na.rm=TRUE)
    down <- min(stD[, 2:(ncol(stD)-1)], na.rm=TRUE)
    plot(t_year+(1:(ncol(stD)-2)), stD[1,2:(ncol(stD)-1)], type="l", axes=FALSE,
         ylim=c(down,up),
         xlab="", ylab="")
    abline(h=0, lty=2, col="grey")
    for(i in 1:nrow(stD)){
      lines(t_year+(1:(ncol(stD)-2)), stD[i,2:(ncol(stD)-1)])
    }
    axis(2, las=2)
    mtext(states3[j])
    if(j %% 4 ==0 | j == length(states3)){
      axis(1)
    } 
    else {
      axis(1, label=FALSE)
    }
    
  }
  mtext("Residuals at timestep after work was completed",
        side=2, outer=TRUE, 
        cex=1.2, line=1)
  mtext("Timestep at end of moving window", side=1, outer=TRUE, line=-.4, cex=1.2)
  dev.off()
  
  
  # make plot of residual over time for each county
  out4v <- out4_sum$full[rowSums(is.na(out4_sum$full)) != (ncol(out4_sum$full)-2),]
  out4v <- data.frame(out4v[,1], scale(out4v[,2:(ncol(out4v)-1)]), out4v[, ncol(out4v)])
  max(out4v[,2:(ncol(out4v)-1)], na.rm=TRUE)
  min(out4v[,2:(ncol(out4v)-1)], na.rm=TRUE)
  plot(t_year+(1:(ncol(out4v)-2)), out4v[1,2:(ncol(out4v)-1)], type="l", ylim=c(-10,10), axes=FALSE,
       xlab="Timestep at end of moving window", ylab="Residuals at timestep after work was completed")
  abline(h=0, lty=2, col="grey")
  for(i in 1:nrow(out4v)){
    lines(t_year+(1:(ncol(out4v)-2)), out4v[i,2:(ncol(out4v)-1)], lty=1)
  }
  axis(1)
  axis(2, las=2)
  mtext(paste0("Capture method = ", cap_methods[4]), line=1)
  
  # break out this method across all states -- for method 4
  states4 <- unique(out4v[,ncol(out4v)])
  mat4 <- matrix(1:30, nrow=6, ncol=5)
  
  pdf(paste0("./fitted_mods/fit20190416/residual_removal_plots/snare.pdf"))
  layout(mat4)
  par(mar=c(3,2,0.1,0.1), oma=c(1,3.5,2,.6))
  for(j in seq_along(states4)){
    stD <- out4v[out4v[,ncol(out4v)] == states4[j],]
    up <- max(stD[, 2:(ncol(stD)-1)], na.rm=TRUE)
    down <- min(stD[, 2:(ncol(stD)-1)], na.rm=TRUE)
    plot(t_year+(1:(ncol(stD)-2)), stD[1,2:(ncol(stD)-1)], type="l", axes=FALSE,
         ylim=c(down,up),
         xlab="", ylab="")
    abline(h=0, lty=2, col="grey")
    for(i in 1:nrow(stD)){
      lines(t_year+(1:(ncol(stD)-2)), stD[i,2:(ncol(stD)-1)])
    }
    axis(2, las=2)
    mtext(states4[j])
    if(j %% nrow(mat4) ==0 | j == length(states4)){
      axis(1)
    } 
    else {
      axis(1, label=FALSE)
    }
    
  }
  mtext("Residuals at timestep after work was completed",
        side=2, outer=TRUE, 
        cex=1.2, line=1)
  mtext("Timestep at end of moving window", side=1, outer=TRUE, line=-.4, cex=1.2)
  dev.off()
  
  # make plot of residual over time for each county
  out5v <- out5_sum$full[rowSums(is.na(out5_sum$full)) != (ncol(out5_sum$full)-2),]
  out5v <- data.frame(out5v[,1], scale(out5v[,2:(ncol(out5v)-1)]), out5v[, ncol(out5v)])
  max(out5v[,2:(ncol(out5v)-1)], na.rm=TRUE)
  min(out5v[,2:(ncol(out5v)-1)], na.rm=TRUE)
  plot(t_year+(1:(ncol(out5v)-2)), out5v[1,2:(ncol(out5v)-1)], type="l", ylim=c(-10,10), axes=FALSE,
       xlab="Timestep at end of moving window", ylab="Residuals at timestep after work was completed")
  abline(h=0, lty=2, col="grey")
  for(i in 1:nrow(out5v)){
    lines(t_year+(1:(ncol(out5v)-2)), out5v[i,2:(ncol(out5v)-1)], lty=1)
  }
  axis(1)
  axis(2, las=2)
  mtext(paste0("Capture method = ", cap_methods[5]), line=1)
  
  # break out this method across all states
  states5 <- unique(out5v[,ncol(out5v)])
  mat5 <- matrix(1:36, nrow=6, ncol=6)
  
  pdf(paste0("./fitted_mods/fit20190416/residual_removal_plots/trap.pdf"))
  layout(mat5)
  par(mar=c(3,2,0.1,0.1), oma=c(1,3.5,2,.6))
  for(j in seq_along(states5)){
    stD <- out5v[out5v[,ncol(out5v)] == states5[j],]
    up <- max(stD[, 2:(ncol(stD)-1)], na.rm=TRUE)
    down <- min(stD[, 2:(ncol(stD)-1)], na.rm=TRUE)
    plot(t_year+(1:(ncol(stD)-2)), stD[1,2:(ncol(stD)-1)], type="l", axes=FALSE,
         ylim=c(down,up),
         xlab="", ylab="")
    abline(h=0, lty=2, col="grey")
    for(i in 1:nrow(stD)){
      lines(t_year+(1:(ncol(stD)-2)), stD[i,2:(ncol(stD)-1)])
    }
    axis(2, las=2)
    mtext(states5[j])
    if(j %% 6 ==0 | j == length(states5)){
      axis(1)
    } 
    else {
      axis(1, label=FALSE)
    }
    
  }
  mtext("Residuals at timestep after work was completed",
        side=2, outer=TRUE, 
        cex=1.2, line=1)
  mtext("Timestep at end of moving window", side=1, outer=TRUE, line=-.4, cex=1.2)
  dev.off()
}

#----- look at slope by county and across windows
t_year <- 33 # = number of timesteps: 66 in 5 years. I might want to make this shorter
n_windows <- length(timesteps)-t_year
slp <- slp1 <- slp2 <- slp3 <- slp4 <- slp5 <- matrix(NA, nrow=length(u_fips), ncol=n_windows)
for(i in seq_along(u_fips)){
  cntyD <- dens[dens$FIPS == u_fips[i], ]
  #cntyD <- sub_method[sub_method$FIPS.x == u_fips[i], ]
  for(j in 1:n_windows){
    t_use <- timesteps[j:(j+t_year)]
    cntyDT2 <- cntyD[which(cntyD$timestep %in% t_use),]
    
    # calculate slope
    if(sum(!is.na(cntyDT2$med))>1){
      if(sum(is.na(cntyDT2$method)) != nrow(cntyDT2)){
        lm1 <- lm(cntyDT2$med ~ cntyDT2$timestep)
        slp[i,j] <- lm1$coefficients[2]
      }
    }
  }
  if(i %% 100 == 0){
    print(paste0(i," of ", length(u_fips), " counties complete"))
  }
}


# plot slopes
densW <- dens3[!is.na(dens3$med.y),]
states <- unique(densW$state.x)[order(unique(densW$state.x))]
mat1 <- matrix(1:36, nrow=6, ncol=6)
rug_cols <- c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e')

#pdf(paste0("./fitted_mods/fit20190416/slopes_33timesteps.pdf"))
layout(mat1)
par(mar=c(3,2,0.1,0.1), oma=c(1,3.5,2,.6))

slp2 <- data.frame(slp, u_fips)
slp3 <- merge(slp2, all_state, by.x="u_fips", by.y="fips")
for(i in 1:length(states)){
  stD <- slp3[slp3$state.x == states[i],]
  std2 <- stD[rowSums(is.na(stD)) != (ncol(stD)-2), ]
  if(nrow(std2) > 0){
    plot((t_year+1):173, unlist(std2[1,(2:(ncol(std2)-2))]), type="l", 
         ylim=c(min(std2[,(2:(ncol(std2)-2))], na.rm=TRUE), max(std2[,(2:(ncol(std2)-2))], na.rm=TRUE)),
         axes=FALSE)
    for(j in 2:nrow(std2)){
      lines((t_year+1):173, unlist(std2[j,(2:(ncol(std2)-2))]))
    }
    axis(1)
    mtext(states[i])
    abline(h=0, col="dark grey", lty=3, lwd=1.5)
  }
  
  # adding a rug for when work was done by each method
  # I can use out1_md if I run this using [work] (instead of [work+1])
  for(k in 5:1){
    assign(paste0("workD", k), eval(parse(text=paste0("out", k, "_sum$full[out", k, "_sum$full$state.x == states[i],]"))))
    assign(paste0("workD_2_", k), eval(parse(text=paste0("workD", k, "[rowSums(is.na(workD", k,")) != (ncol(workD", k,")-2), (2:(ncol(workD", k,")-2))]"))))
    work_cols <- which(!is.na(eval(parse(text=paste0("workD_2_", k)))), arr.ind=TRUE)
    if(ncol(work_cols)>1){
      rug(work_cols[,2], col=rug_cols[k], ticksize=.05, lwd=.5)
    }else{
      rug(work_cols, col=rug_cols[k], ticksize=.05, lwd=.5)
    }
    
  }
}
mtext(paste0("slope of abundance over time for ", t_year, " timesteps"), 2, outer=TRUE, cex=1.5)
mtext("timestep", 1, outer=TRUE, cex=1.5)
dev.off()

#---------
# looking at just heavy control work

# for each county at each timestep, get the amount of effort
ag_eff <- aggregate(y_sum$effort, by=list(y_sum$FIPS2, y_sum$timestep, y_sum$method),
                    FUN=sum)
colnames(ag_eff) <- c("FIPS2", "timestep", "method", "effort")
hist(ag_eff$effort, breaks=100)
quantile(ag_eff$effort)
# sepaprate by heavy control opeartions
heavy_q <- .25 # quantile of amount of control efforts deemed heavy
heavy <- quantile(ag_eff$effort, heavy_q) # control work >= to this amount will be used
ag_eff2 <- ag_eff[ag_eff$effort >= heavy,]
dens_eff <- merge(dens, ag_eff2, by=c("FIPS2", "timestep"))

# set up some plot parameters
color.gradient <- function(x, colors=c("red","yellow","green"), colsteps=100) {
  return(colorRampPalette(colors,interpolate = c("spline"))(colsteps)[ findInterval(x, seq(min(x),max(x), length.out=colsteps)) ] 
  )
}
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
num.colors<-20
blue.end.pnts <- c("#c6dbef","#08306b")
red.end.pnts <- c("#ffeda0","#730022")

#----- look at slope by county and across windows
slopeBeforeAfterControl <- function(
  t_year=3,
  heavy_q =.1,
  percent_change = TRUE, # using percent change for before/after heavy control?
  how_long_negative = FALSE,
  output = "diff_cnty"
){
  # separate data for heavy control work
  heavy <- quantile(ag_eff$effort, heavy_q) # control work >= to this amount will be used
  ag_eff2 <- ag_eff[ag_eff$effort >= heavy,]
  dens_eff <- merge(dens, ag_eff2, by=c("FIPS2", "timestep"))
  
  timesteps <- min(timestep_df$timestep):max(timestep_df$timestep)
  n_windows <- length(timesteps)-t_year
  u_fips <- unique(dens_eff$FIPS) # this is different here because it is only places with high effort
  slp_bef <- slp_aft <- time_negative <- matrix(NA, nrow=length(u_fips), ncol=200) # ncol is max number of instances of heavy control in a county
  for(i in seq_along(u_fips)){
    cntyD_eff <- dens_eff[dens_eff$FIPS == u_fips[i], ]
    cntyD <- dens[dens$FIPS == u_fips[i], ]
    for(j in 1:nrow(cntyD_eff)){
      time1 <- dens_eff$timestep[j]
      if(how_long_negative == FALSE){
        t_before <- (time1 - t_year) : (time1 - 1)
        t_after <- (time1 +1):(time1 + t_year)
        if((min(t_before) >=1) & max(t_after) <= max(timesteps)){
          cntyD_bef <- cntyD[cntyD$timestep %in% t_before,]
          lm_bef <- lm(cntyD_bef$med ~ cntyD_bef$timestep)
          slp_bef[i,j] <- lm_bef$coefficients[2]
          
          cntyD_aft <- cntyD[cntyD$timestep %in% t_after,]
          lm_aft <- lm(cntyD_aft$med ~ cntyD_aft$timestep)
          slp_aft[i,j] <- lm_aft$coefficients[2]
          
        }
      }
      # if calculating how long a slope stays negative in a county after control
      # set t_year == 2 for this
      if(how_long_negative){
        slope_after <- matrix(0, nrow=nrow(cntyD_eff), ncol=max(timesteps)) # this is a temporary matrix with the slopes
        t_after <- (time1 +1):(time1 + t_year)
        if(max(t_after) > max(timesteps)){
          slope_after[j,] <- rep(1, max(timesteps)) # I should not need to do this
        } else {
          cntyD_aft <- cntyD[cntyD$timestep %in% t_after,]
          lm_aft <- lm(cntyD_aft$med ~ cntyD_aft$timestep)
          slope_after[j,time1] <- lm_aft$coefficients[2]
          if(slope_after[j,time1] < 0){
            for(k in (t_year+1):(max(timesteps)-time1)){
              if(slope_after[j,(time1+k-t_year-1)] <0 ){
                t_after <- (time1 +1):(time1+k)
                cntyD_aft <- cntyD[cntyD$timestep %in% t_after,]
                lm_aft <- lm(cntyD_aft$med ~ cntyD_aft$timestep)
                slope_after[j,(time1+k-t_year)] <- lm_aft$coefficients[2] # this contains the negative slopes for each heavy control work incident in a county, but I need to summarize this by county
              } else {
                slope_after[j,(time1+k-t_year):max(timesteps)] <- 1 # once it goes positive, just fill in the rest with 1s to complete the loop
              }
            } # close k
          }
          num_neg_times <- sum(slope_after[j, ] < 0, na.rm=TRUE) + t_year -1 # number of timesteps where it stayed negative
          time_negative[i,j] <- ifelse(num_neg_times == t_year-1, 0, num_neg_times)
        } 
        
      }# close if how long negative
      
    } # close j
    
    # print progress
    if(i %% 100 == 0){
      print(paste0(i," of ", length(u_fips), " counties complete"))
    }
  }
  
  # control output by options
  if(how_long_negative){
    time_neg <- data.frame(time_negative, u_fips)
    time_neg2 <- merge(time_neg, all_state, by.x="u_fips", by.y="fips")
    return(time_neg2)
  }else{
    # calculate difference before and after control
    if(percent_change){
      slp_df <- (slp_bef - slp_aft)/slp_bef * 100
    } else{
      slp_df <- slp_bef-slp_aft
    }
    
    slp_bef2 <- data.frame(slp_bef, u_fips)
    slp_bef3 <- merge(slp_bef2, all_state, by.x="u_fips", by.y="fips")
    slp_aft2 <- data.frame(slp_aft, u_fips)
    slp_aft3 <- merge(slp_aft2, all_state, by.x="u_fips", by.y="fips")
    slp_df <- slp_bef-slp_aft
    slp_df2 <- data.frame(slp_df, u_fips)
    slp_df3 <-  merge(slp_df2, all_state, by.x="u_fips", by.y="fips")
    
    
    
    
    # compare slope before control with slope after
    # plot(slp_bef3[,3], slp_aft3[,3])
    # quantile(slp_bef, c(0.025, 0.5, 0.975), na.rm=TRUE)
    # quantile(slp_aft, c(0.025, 0.5, 0.975),na.rm=TRUE)
    # quantile(slp_df, c(0.025, 0.5, 0.975),na.rm=TRUE)
    # 
    # compare slope before control with slope after
    states_BA <- unique(slp_bef3$state.x)
    BA_ag <- matrix(NA, nrow=length(states_BA), ncol=9)
    for(i in seq_along(states_BA)){
      st_bef <- slp_bef3[slp_bef3$state.x==states_BA[i], 2:(ncol(slp_bef3)-2)]
      before <- quantile(st_bef, c(0.025,0.5,0.975), na.rm=TRUE)
      st_aft <- slp_aft3[slp_aft3$state.x==states_BA[i], 2:(ncol(slp_aft3)-2)]
      after <- quantile(st_aft, c(0.025,0.5,0.975), na.rm=TRUE)
      st_df <- slp_df3[slp_df3$state.x==states_BA[i], 2:(ncol(slp_df3)-2)]
      diff <- quantile(st_df, c(0.025,0.5,0.975), na.rm=TRUE)
      BA_ag[i,] <- c(before, after, diff)
    }
    qs <- c("lo", "med", "hi")
    colnames(BA_ag) <- c(paste0("before_", qs), paste0("after_", qs), paste0("diff_", qs))
    BA_ag2 <- data.frame(states_BA, BA_ag)
    data.frame(BA_ag2$states_BA, BA_ag2$before_med, BA_ag2$after_med, BA_ag2$diff_med)
    #write_csv(BA_ag2, "./fitted_mods/fit20190416/effects_of_control/slopeBeforeAfterControl_3timesteps_25quant.csv")
    write_csv(BA_ag2, paste0("./fitted_mods/fit20190416/effects_of_control/slopeBeforeAfterControl_", t_year, "timesteps_", (heavy_q * 100), "quant.csv"))
    print(data.frame(BA_ag2$states_BA, BA_ag2$before_med, BA_ag2$after_med, BA_ag2$diff_med))
    
    # look at difference in slope in individual counties
    slp_df_cnty <- apply(slp_df, 1, FUN=function(x) median(x, na.rm=TRUE))
    dif_cnty <- data.frame(u_fips, slp_df_cnty, slp_df3[,ncol(slp_df3)])
    colnames(dif_cnty) <- c("FIPS", "med.slp.dif", "state")
    
    # if(output == "diff_cnty"){
    #   return(dif_cnty)
    # }
    
    #-- make map plot
    # make a map of change before and after control
    data(county.fips)
    y <- dif_cnty$med.slp.dif
    df_col_cnty <- data.frame(FIPS=dif_cnty$FIPS,
                              value=dif_cnty$med.slp.dif,
                              color=grey(range01(y)),
                              fips = as.integer(as.numeric(as.character(dif_cnty$FIPS))))
    counties <- merge(county.fips, df_col_cnty, by="fips", all.x=TRUE)
    
    # set up pallete
    # determine zero point
    lwr.mid <- round(mean(y<0)*num.colors)
    upper.mid <- num.colors-lwr.mid
    # Blue Palette (cut point at .60 so 0 to 12)
    blue.palette <- color.gradient(1:lwr.mid,colors=blue.end.pnts)
    # Red Palette (.6 to 1)
    red.palette <- color.gradient(1:upper.mid,colors=red.end.pnts)
    # Combine Palettes into one
    red.blue.palette <- c(rev(blue.palette),red.palette)
    
    # apply colors to data
    h1 <- c(seq(min(dif_cnty$med.slp.dif), -0.000001, length.out=lwr.mid),
            seq( 0.000001, max(dif_cnty$med.slp.dif),length.out=upper.mid))
    ints <- findInterval(counties$value, h1)
    newcol <- ifelse(is.na(ints), rgb(0,0,0,alpha=0),
                     red.blue.palette[ints])
    
    if(output == "diff_cnty"){
      return(list(dif_cnty=dif_cnty,
                  newcol=newcol
      ))
    }
  }
  
  
  
  
  # pdf(paste0("./fitted_mods/fit20190416/effects_of_control/slopeBeforeAfterControl_", t_year, "timesteps_", (heavy_q * 100), "quant.pdf"))
  # plot1 <- map('state',lwd=2,col="black") # states
  # map('county', fill=TRUE, col=newcol,lty=0, add=TRUE) # data
  # title(main="Change in wild pig abundance after control", line=1)
  # print(plot1)
  #  dev.off()
  
}
dif_cnty <- slopeBeforeAfterControl(heavy_q =.95)

# make plots of different quantiles for heavy work and different amounts of time pre and post control
# put this into a loop because I can't make a map within a fucntion
qs <- c(0, .1, .25, .5, .75, .9)
times <- 2:4
for(i in seq_along(qs)){
  for(j in seq_along(times)){
    dif_cnty <- slopeBeforeAfterControl(t_year=times[j],
                                        heavy_q =qs[i],
                                        percent_change = TRUE)
    
    pdf(paste0("./fitted_mods/fit20190416/effects_of_control/slopeBeforeAfterControl_", times[j], "timesteps_", (qs[i] * 100), "quant.pdf"))
    map('state',lwd=2,col="black") # states
    map('county', fill=TRUE, col=dif_cnty$newcol,lty=0, add=TRUE) # data
    title(main="Change in wild pig abundance after control", line=1)
    #legend("bottomright",legend=c("percent change after control"))
    #plot(rep(1,num.colors),1:num.colors, col=red.blue.palette, pch=15, cex=5, axes=FALSE,xlab="",ylab="")
    dev.off()
    
  }
}

# look at how long slope stays negative
t_neg <- slopeBeforeAfterControl(t_year=2, heavy_q =0, how_long_negative = TRUE)
apply(t_neg[,2:(ncol(t_neg)-2)], 1, FUN=function(x) median(x, na.rm=TRUE))

# make map of this
data(county.fips)
t_neg_df <- data.frame(FIPS=t_neg$u_fips,
                       value=apply(t_neg[,2:(ncol(t_neg)-2)], 1, FUN=function(x) max(x, na.rm=TRUE)),
                       fips = as.integer(as.numeric(as.character(t_neg$u_fips))))
counties <- merge(county.fips, t_neg_df, by="fips", all.x=TRUE)
green.palette <- c("grey", color.gradient(1:20, c('#e5f5f9','#2ca25f')))
green.palette <- c("grey", color.gradient(1:20, c('#74c476','#00441b')))
purple.pallette <- color.gradient(1:20, c('#9e9ac8','#3f007d'))
green.purple.pallette <- c("grey", (green.palette), purple.pallette)
plot( rep(1,21),1:21, col=green.palette, pch=15, cex=5, axes=FALSE,xlab="",ylab="")
h1 <- seq(0, max(t_neg_df$value), length.out = 20)
ints <- findInterval(counties$value, h1)
newcol <- ifelse(is.na(ints), rgb(0,0,0,alpha=0),
                 green.palette[ints])

# map plot
pdf(paste0("./fitted_mods/fit20190416/effects_of_control/monthsNegSlope_max.pdf"))
mat1 <- matrix(c(1,1,2,2,1,1,2,2), byrow = TRUE, nrow=2, ncol=4)
layout(mat1)
maps::map('state',lwd=2,col="black") 
maps::map('county', fill=TRUE, col=newcol,lty=0, add=TRUE) # data
title(main="Maximum number of months that abundance\n decreases following control", line=1)
#legend("bottomright", legend=c(0, seq(1,21, by=2)), col=green.palette[c(TRUE, FALSE)], horiz = FALSE, fill = newcol)
plot( rep(1,21),1:21, col=c("grey", green.palette), pch=15, cex=5, axes=FALSE,xlab="",ylab="")
axis(2, at=c(1, 6, 11, 16, 21), labels = c(0, 42, 83, 124, 170), las=2, line=-10, tick=FALSE)
dev.off()

library(english)
library(toOrdinal)

# make map plot for each of the first X control efforts in a county
pdf(paste0("./fitted_mods/fit20190416/effects_of_control/monthsNegSlope_controlEffortNum", i,".pdf"))
for(i in 1:10){
  # set up data for this control effort
  t_neg_df <- data.frame(FIPS=t_neg$u_fips,
                         value=t_neg[,(i+1)],
                         fips = as.integer(as.numeric(as.character(t_neg$u_fips))))
  counties <- merge(county.fips, t_neg_df, by="fips", all.x=TRUE)
  h1 <- seq(0, max(t_neg_df$value, na.rm=TRUE), length.out = 20)
  ints <- findInterval(counties$value, h1)
  newcol <- ifelse(is.na(ints), rgb(0,0,0,alpha=0),
                   green.palette[ints])
  # plot
  pdf(paste0("./fitted_mods/fit20190416/effects_of_control/monthsNegSlope_controlEffortNum", i,".pdf"))
  mat1 <- matrix(c(1,1,2,2,1,1,2,2), byrow = TRUE, nrow=2, ncol=4)
  layout(mat1)
  maps::map('state',lwd=2,col="black") 
  maps::map('county', fill=TRUE, col=newcol,lty=0, add=TRUE) # data
  title(main=paste0("Number of months that abundance decreases\n following the ", toOrdinal(i), " control operation in county"), line=1)
  #legend("bottomright", legend=c(0, seq(1,21, by=2)), col=green.palette[c(TRUE, FALSE)], horiz = FALSE, fill = newcol)
  plot( rep(1,21),1:21, col=c("grey", green.palette), pch=15, cex=5, axes=FALSE,xlab="",ylab="")
  axis(2, at=c(1, 6, 11, 16, 21), labels = c(0, 42, 83, 124, 170), las=2, line=-10, tick=FALSE)
  dev.off()
}

# make mulipanel plot for effect over time

pdf(paste0("./fitted_mods/fit20190416/effects_of_control/monthsNegSlope_afterControl.pdf"))
mat2 <- matrix(c(1,2,3,4,5,
                 6,7,8,9,10), byrow=T, nrow=2)
layout(mat2)
for(i in 1:9){
  # set up data for this control effort
  t_neg_df <- data.frame(FIPS=t_neg$u_fips,
                         value=t_neg[,(i+1)],
                         fips = as.integer(as.numeric(as.character(t_neg$u_fips))))
  counties <- merge(county.fips, t_neg_df, by="fips", all.x=TRUE)
  h1 <- seq(0, max(t_neg_df$value, na.rm=TRUE), length.out = 20)
  ints <- findInterval(counties$value, h1)
  newcol <- ifelse(is.na(ints), rgb(0,0,0,alpha=0),
                   green.palette[ints])
  # plot
  #pdf(paste0("./fitted_mods/fit20190416/effects_of_control/monthsNegSlope_controlEffortNum", i,".pdf"))
  #mat1 <- matrix(c(1,1,2,2,1,1,2,2), byrow = TRUE, nrow=2, ncol=4)
  #layout(mat1)
  maps::map('state',lwd=2,col="black") 
  maps::map('county', fill=TRUE, col=newcol,lty=0, add=TRUE) # data
  title(main=paste0("Number of months that abundance decreases\n following the ", toOrdinal(i), " control operation in county"), line=1)
}
plot( rep(1,21),1:21, col=c("grey", green.palette), pch=15, cex=5, axes=FALSE,xlab="",ylab="")
axis(2, at=c(1, 6, 11, 16, 21), labels = c(0, 42, 83, 124, 170), las=2, line=-10, tick=FALSE)
dev.off()

#-------------------------------------------
# slope/percent change after control until next control event
perC <- function(x,y){
  (x-y)/x
}

#dens4 <- dens3[!(is.na(dens3$effort)), ]
timesteps <- timestep_df$timestep
firstPositive <- matrix(NA, nrow=length(u_fips), ncol=length(timesteps))
btwnControls <- function(
  
  slope=FALSE
){
  for(i in seq_along(u_fips)){
    cntyD <- dens3[dens3$fips == u_fips[i], ]
    cntyD2 <- cntyD[!is.na(cntyD$effort), ]
    if(nrow(cntyD2) >1){ # if more than one control effort
      controlTimes <- unique(cntyD2$timestep)
      if(max(controlTimes) == max(timesteps)){ # if there is control in the last timestep, can't measure change after this.
        controlTimes <- controlTimes[controlTimes != max(controlTimes)]
      }
      for(j in seq_along(controlTimes)){ # j is current timestep
        if(j != max(seq_along(controlTimes))){
          first <- controlTimes[j]
          last <- controlTimes[j+1]
          
          perChange <- rep(NA, (last-first))
          for(k in (first+1):last){
            perChange[which((first+1):last %in% k)] <- 
              perC(x=mean(cntyD[cntyD$timestep==first, "med.x"]),
                   y=mean(cntyD[cntyD$timestep==k, "med.x"]))
          
          } 
          firstPositive[i,j] <- which(perChange > 0)[1] # the first positive percent change from that control time
          
        } else{
          first <- controlTimes[j]
          last <- max(timesteps)
          perChange <- rep(NA, (last-first))
          for(k in (first+1):last){
            perChange[which((first+1):last %in% k)] <- 
              perC(x=mean(cntyD[cntyD$timestep==first, "med.x"]),
                   y=mean(cntyD[cntyD$timestep==k, "med.x"]))
            
          }
          firstPositive[i,j] <- which(perChange > 0)[1] # the first positive percent change from that control time
          
        }
        
      } # close j
    }
    if(i %% 100 == 0){
      print(paste0(i," of ", length(u_fips), " counties complete"))
    }
  } # close i
  return(firstPositive)
}
firstPositive <- btwnControls()
fp2 <- firstPositive[rowSums(is.na(firstPositive)) != ncol(firstPositive), ]
timeReduced <- data.frame(firstPositive - 1, u_fips)

# add state data
tR2 <- merge(timeReduced, all_state, by.x="u_fips", by.y="fips")
CA <- tR2[tR2$state.x=="California",]

# make map of mean time reduced
green.palette <- c("grey", color.gradient(1:20, c('#e5f5f9','#2ca25f')))
green.palette <- c("grey", color.gradient(1:20, c('#74c476','#00441b')))

# set up data for this control effort
t_red_df <- data.frame(FIPS=u_fips,
                       value=rowMeans(timeReduced[,1:length(timesteps)], na.rm=TRUE),
                       #value=apply(timeReduced, 1, function(x) median(x, na.rm=TRUE)),
                       fips = u_fips,
                       state = tR2$state.x)
#counties <- merge(county.fips, t_red_df, by="fips", all.x=TRUE)
h1 <- seq(0, max(t_red_df$value, na.rm=TRUE), length.out = 20)
h1 <- seq(0, 5, #max(t_red_df$value, na.rm=TRUE), 
          length.out = 20)
h1[length(h1)] <- max(t_red_df$value, na.rm=TRUE)
#h1 <- c(seq(0, 5, by=.5), max(t_red_df$value, na.rm=TRUE))
ints <- findInterval(counties$value, h1)
ints <- findInterval(t_red_df$value, h1)
newcol <- ifelse(is.na(ints), rgb(0,0,0,alpha=0),
                 green.palette[ints])
hist(ints, breaks=200)
CA <- t_red_df[t_red_df$state=="California",]

# plot
mat <- matrix(c(1,2), nrow=1)
#pdf("./fitted_mods/fit20190416/effects_of_control/mean_months_neg_slope.pdf")
maps::map('state',lwd=2,col="black", main="Mean number of months with \nnegative slope after control") 
maps::map('county', fill=TRUE, col=newcol,lty=0, add=TRUE) # data
mtext("Mean number of consecutive months with \ndecreasing population size after control", cex=2)
#legend("bottomright", legend=0:20, fill=green.palette)
dev.off()

pdf("./fitted_mods/fit20190416/effects_of_control/mean_months_neg_slope_legend.pdf")
par(font=2)
plot(rep(1,21),1:21, col=c(green.palette), pch=15, cex=5, axes=FALSE,xlab="",ylab="")
axis(4, at=c(.8, 5, 8.5, 12.5, 16.5, 21.5), labels = c(0, 1,2, 3, 4, expression("">=5)), 
     las=2, line=-10, tick=FALSE, pos=1.05, cex.axis=2)
dev.off()
  #---------------------------------------------------------------
  # archive
  #---------------------------------------------------------------
  # find max number of surveys at at st unit
  num_st <- rep(NA, n_st)
  for(i in 1:n_st){
    num_st[i] <- length(st_index[[i]])
  }
  max(num_st) #==22
  
  # how to get the indeces in a table
  set.seed(100)
  n <- 1000
  samps <- sample(1:10, n, replace=TRUE)
  tbl1 <- table(samps)
  which(samps %in% tbl1)
  order(samps)[!duplicated(sort(samps))]
  split(seq_along(samps), samps)
  
  # make sure subtraction works how I want 
  mat <- matrix(rpois(10, 2), nrow=10, ncol=20)
  vec <- t(replicate(nrow(mat), 1:20))
  (mat-vec)/vec
  apply(mat, 2, function(x){x-vec})
  sweep(mat[,1:ncol(mat)], 2, function(x) {x-vec})
  
  # test how to extract data I want for each survey
  n <- 1000
  y_obs <- rpois(n, 3)
  n_surveys <- 200
  survey_num <- 1:n_surveys
  survey_idx <- sample(1:n_surveys, size=n, replace=TRUE)
  n_surv <- table(survey_idx)
  y_surv <- list()
  o_idx <- which(survey_idx %in% survey_num[i])
  surv_list <- list()
  for(i in 1:n_surveys){
    surv_list[[i]] <- which(survey_idx %in% survey_num[i])
  }
  
  # check that adding spatiotemporal data to y_sum worked
  unique(st_d$county_index)
  data.frame(y_sum[y_sum$county_idx== 6,])
  data.frame(st_d[st_d$county_index == 6, ])
  
  # I need to break this out by timestep. Otherwise, there are multiple values per county
  df_map <- data.frame(county=as.integer(d2009$FIPS), value=d2009$medDens, 
                       color=grey(d2009$medDens/max(d2009$medDens)))
  counties <- merge(county.fips, df_map, by.x="fips", by.y="county", all.x=TRUE)
  cc <- counties[!is.na(counties$value), ]
  dim(cc)
  #postscript("./fig/densMap_20190410.eps")
  #map("county", fill=TRUE, col=counties$color)
  #dev.off()
  
  # plot density
  par(mar=c(5,8,4,3))
  plot(quant[3,1], 1, ylim=c(.5, 6.5),xlim=c(min(quant), max(quant)),
       axes=FALSE, 
       xlab="", 
       ylab="", cex=1, pch=16)
  segments(x0=quant[2,1], x1=quant[4,1], y0=1, lwd=lwd50)
  segments(x0=quant[1,1], x1=quant[5,1], y0=1, lwd=lwd95)
  points(quant[3,1], 1, cex=cex.pch, pch=16)
  for(i in 2:6){
    segments(x0=quant[2,i], x1=quant[4,i], y0=i, lwd=lwd50)
    segments(x0=quant[1,i], x1=quant[5,i], y0=i, lwd=lwd95)
    points(quant[3,i], i, cex=cex.pch, pch=16)
  }
  axis(1)
  axis(2, at=1:6, labels=c("all methods", "firearms", "fixed wing", "helicopter", "snare", "trap"), las=2)
  mtext(paste0("reduction of wild pig density at a property\n by next removal effort within ", max_months, " months"), 1,
        line=3)
  
  # look at a specific county and abundance vs effort 
  max_months <- 3 # the max number of months in which to assess the effect of removal
  u_fips <- unique(y_sum$FIPS)
  for(i in u_fips){
    tmp <- y_sum[y_sum$FIPS == i,]
    if(nrow(tmp) > 200){
      print(paste0(i, "--" ,nrow(tmp)))
    }
  }
  
  # do this for all properties
  max_months <- 1 # the max number of months in which to assess the effect of removal
  min_months <- 3 # the minimum number of months in which to assess the effect of removal
  reduce_df <- matrix(NA, nrow=length(props), ncol=6)
  props <- unique(y_sum$agrp_prp_id)
  cap_methods <- c("firearms", "fixed wing", "helicopter", "snare", "trap")
  for(i in seq_along(props)){
    prop1 <- y_sum[y_sum$agrp_prp_id == props[i],]
    if(length(unique(prop1$timestep)) < 2){ # if work was only conducted at one timestep, ignore this property
      reduce_df[i, ] <- rep(NA, ncol(reduce_df))
    } else{
      diffs <- ave(prop1$timestep, FUN=function(x) c(NA,diff(x)))
      short_idxs <- which(diffs <= max_months & diffs >= min_months)
      if(length(short_idxs) < 2 ){ # if there are not 2 instances where the time difference between removal is approporiate, ignore this property
        reduce_df[i, ] <- rep(NA, ncol(reduce_df))
      } else{
        if(max(short_idxs) == nrow(prop1)){ # I cannot calculate the difference if it is the last year of work being done at the property
          short_idxs <- short_idxs[-(max(short_idxs))]
        } else{
          if(length(short_idxs) < 2) { # if there are not 2 instances where the time difference between removal is approporiate, ignore this property
            reduce_df[i, ] <- rep(NA, ncol(reduce_df))
          } else{
            prop1$diff <- rep(NA, nrow(prop1))
            # the equation for calculating difference
            difference <- (prop1$med_dens[short_idxs]-prop1$med_dens[short_idxs+1])/(prop1$med_dens[short_idxs]+1)
            prop1$diff[short_idxs] <- difference
            # slope for all methods combined
            lm1 <- lm(prop1$diff ~ prop1$scaled_effort)
            reduce_df[i,1] <- lm1$coefficients[2]
            
            # now look at the slope for each method
            for(j in seq_along(cap_methods)){
              prop1_meth <- prop1[prop1$method==j, ]
              prop1_meth <- prop1_meth[!is.na(prop1_meth$diff), ]
              if(nrow(prop1_meth) < 2){
                reduce_df[i,(j+1)] <- NA
              } else{
                lm2 <- lm(prop1_meth$diff ~ prop1_meth$scaled_effort)
                reduce_df[i,(j+1)] <- lm2$coefficients[2]
              }
            } 
          }
        }
        
      }
      
    }
  }
  
  
  # look at slope and residuals in just Texas
  densTX <- dens3[dens3$state.x=="Texas",]
  u_fips <- unique(densTX$FIPS)[!is.na(unique(densTX$FIPS))]
  u_fips <- unique(densTX$FIPS.x)[!is.na(unique(densTX$FIPS.x))]
  out_mn <- out_md <- out_direction <-  matrix(NA, nrow=length(u_fips), ncol=n_windows)
  for(i in seq_along(u_fips)){
    cntyD <- densTX[densTX$FIPS.x == u_fips[i], ]
    for(j in 1:n_windows){
      t_use <- timesteps[j:(j+t_year)]
      cntyDT2 <- cntyD[which(cntyD$timestep %in% t_use),]
      if(sum(!is.na(cntyDT2$med.x))>1){
        if(sum(is.na(cntyDT2$method)) != nrow(cntyDT2)){
          lm1 <- lm(cntyDT2$med.x ~ cntyDT2$timestep)
          work_idxs <- which(!is.na(cntyDT2$method)) 
          direction <- ifelse(lm1$residuals>0, 1, -1)
          out_mn[i, j] <- mean(lm1$residuals[(work_idxs+1)], na.rm=TRUE) # na.rm to deal with situations when work was done in the last timestep
          out_md[i, j] <- median(lm1$residuals[work_idxs+1], na.rm=TRUE)
          # out_direction[i, j] <- mean(direction[work_idxs+1])
        }
      }
    }
  }
  # do this for all states
  u_fips <- unique(dens3$FIPS.x)
  out_mn <- out_md <- out_direction <-  matrix(NA, nrow=length(u_fips), ncol=n_windows)
  for(i in seq_along(u_fips)){
    cntyD <- dens3[dens3$FIPS.x == u_fips[i], ]
    for(j in 1:n_windows){
      t_use <- timesteps[j:(j+t_year)]
      cntyDT2 <- cntyD[which(cntyD$timestep %in% t_use),]
      if(sum(!is.na(cntyDT2$med.x))>1){
        if(sum(is.na(cntyDT2$method)) != nrow(cntyDT2)){
          lm1 <- lm(cntyDT2$med.x ~ cntyDT2$timestep)
          work_idxs <- which(!is.na(cntyDT2$method)) 
          direction <- ifelse(lm1$residuals>0, 1, -1)
          #out_mn[i, j] <- mean(lm1$residuals[(work_idxs+1)], na.rm=TRUE) # na.rm to deal with situations when work was done in the last timestep
          out_md[i, j] <- median(lm1$residuals[work_idxs+1], na.rm=TRUE)
          # out_direction[i, j] <- mean(direction[work_idxs+1])
        }
      }
    }
  }
  
