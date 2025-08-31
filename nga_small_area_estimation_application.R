##---Bayesian Hierarchical Population modelling for NIGERIA Using
##--Author: Dr Chris Nnanatu, WorldPop, University of Southampton
#--Date: 5th January 2025
rm(list=ls()) #----Clear the workspace
#packages <- c("raster", "haven", "sf","sp", "tmap","tidyverse","rgdal",
#       "lattice", "gridExtra", "devtools", "rlang", "viridis", "spdep",
#     "car")
#if(length(setdiff(packages, rownames(installed.packages())), type="binary") > 0) {
# install.packages(setdiff(packages, rownames(installed.packages()))) }

#Installing INLA!!
#if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
#  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#}
library(INLA)
library(tmap)
library(tmaptools)
library(sf)
library(sp)
library(tmap)
library(tmaptools)
library(viridis)
library(tidyverse)
library(gridExtra)
library(lattice)
library(MASS)
library(car)
library(raster)


# Set working directory for input and output files
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/NGA/CHRIS_N"
data_path <- paste0(path, "/data/input/pop")#---paths for survey dat_sca:
cov_path <- paste0(path, "/data/input/covariates")#---paths for covariates:
results_path <- paste0(path , "/data/new_output_29.08.2025")

##-------------------------------------------------------------------------------------------------------------
#dim(shp <- st_read(paste0(data_path, "/NGA_NMEP_data_train4.shp"))) # NMEP shapefile
dim(shp <- st_read(paste0(data_path, "/input_data.gpkg"))) # NMEP shapefile


plot(shp["pop"], main="Observed  wards(raw)")
plot(shp["pop_imp"], main="Observed  wards(imputed)")
plot(shp["pop_drp"], main="Observed  wards(dropped)")


# Convert the shapefile to data frame
dim(dat <- data.frame(shp))
dim(dat <- dat %>% dplyr::select(-geom))



#Rename key variables
dat$lga_name <- dat$GEONAME #LGA
dat$wardname <- dat$wrdnm_x # Ward
dat$state_name <- dat$stat_nm # state


# summarise observed counts by state
(state_df <- data.frame(dat %>% group_by(state_name)%>%
                          summarise(Base = sum(pop, na.rm=T), # Base
                                    Imputed = sum(pop_imp, na.rm=T), # Imputed
                                    Dropped = sum(pop_drp, na.rm=T), # Dropped
                                    bld = sum(bld_count, na.rm=T)))) # Building counts

# Summarise by states as shapefile
(state_shp <- shp %>% group_by(stat_nm)%>%
                          summarise(Base = sum(pop, na.rm=T), # Base
                                    Imputed = sum(pop_imp, na.rm=T), # Imputed
                                    Dropped = sum(pop_drp, na.rm=T), # Dropped
                                    bld = sum(bld_count, na.rm=T))) # Building counts
state_shp$Base[state_shp$Base==0]=state_shp$Imputed[state_shp$Imputed==0]=state_shp$Dropped[state_shp$Dropped==0]=NA
plot(state_shp["Base"], main="Observed (Base)")
plot(state_shp["Imputed"], main="Observed (imputed)")
plot(state_shp["Dropped"], main="Observed (dropped)")


# Summarise by LGA as shapefile
shp$ID_lga <- as.numeric(as.factor(paste(shp$stat_nm, shp$GEONAME)))
length(unique(shp$ID_lga))
(lga_shp <- shp %>% group_by(ID_lga, GEONAME, stat_nm)%>%
    summarise(Base = sum(pop, na.rm=T), # Base
              Imputed = sum(pop_imp, na.rm=T), # Imputed
              Dropped = sum(pop_drp, na.rm=T), # Dropped
              bld = sum(bld_count, na.rm=T))) # Building counts


lga_shp$Base[lga_shp$Base==0]=lga_shp$Imputed[lga_shp$Imputed==0]=lga_shp$Dropped[lga_shp$Dropped==0]=NA
plot(lga_shp["Base"], main="Observed (Base)")
plot(lga_shp["Imputed"], main="Observed (imputed)")
plot(lga_shp["Dropped"], main="Observed (dropped)")




# Extract and add the coordinates to the data
shp2 <- as(st_geometry(shp), "Spatial") # 9266 wards
plot(shp2)

dat$lon <- coordinates(shp2)[,1] # longitude
dat$lat <- coordinates(shp2)[,2]  # longitudes
plot(dat$lon, dat$lat)


#------------------------------------------------------------------------------------
# Variable Selection-----------------------------------------------------------

# Covariates scaling function
stdize2 <- function(x)
{
  
  stdz <- (x - mean(x[!is.na(x)]))/sd(x[!is.na(x)])
  return(stdz)
}


# define the population density variable
par(mfrow=c(1,1))

boxplot(dat$dens <- dat$pop/dat$bld_count) # based on the raw data

# more than 30 people per building is considered as outlier
out = 30
boxplot(dat$dens[dat$dens <  out])# Use 16 across data to make it comparable
dat$dens[dat$dens > out] = NA # set identified outliers to NA
boxplot(dat$dens)

dat$dens_imp <- dat$pop_imp/dat$bld_count # based on the imputed data
boxplot(dat$dens_imp[dat$dens_imp <  out])
dat$dens_imp[dat$dens_imp > out] = NA

dat$dens_drp <- dat$pop_drp/dat$bld_count # based on the dropped data
boxplot(dat$dens_drp[dat$dens_drp <  out])
dat$dens_drp[dat$dens_drp > out] = NA
par(mfrow=c(1,1))

# rename the geospatial covariates for ease of use
names <- colnames(dat)
vars <- paste0("x", 1:55)
#vars <- vars[-c(1,14,15,55)] 
dim(dat[,vars]  <- apply(dat[,vars], 2, stdize2)) #apply scaling


covs_dens <- c("x28", "x30","x45","x48","x53", "x54") # besf fit covariates identified from stepwise regression


#  Build Mesh for geostatistical modelling
coords2 <- cbind(dat$lon, dat$lat) # mesh
non_convex_bdry <- inla.nonconvex.hull(coords2, -0.035, -0.05, resolution = c(100, 100))
mesh <- inla.mesh.2d(boundary = non_convex_bdry, max.edge=c(0.2,1),
                     offset = c(0.2, 0.7),
                     cutoff = 0.2)
# visualise
par(mfrow=c(1,1))
plot(mesh)
plot(shp2, add=T)
points(dat$lon, dat$lat,
       col="red", pch=15, cex=0.3)

mesh$n # 1698 mesh nodes


###---Build projector matrix A
A <-inla.spde.make.A(mesh=mesh,loc=as.matrix(coords2));dim(A)

##---Create the SPDE
spde <- inla.spde2.matern(mesh, alpha=2)

##----specify the observation indices for estimation
iset <- inla.spde.make.index(name = "s", spde$n.spde)

class(dat$bld_count <- as.integer(round(dat$bld_count )))



# prepare data for initial step modelling
dat$eps <- 1:nrow(dat)
dat$ward_f <- factor(as.numeric(as.factor(dat$wardname)))
dat$set_typ <- factor(as.numeric(as.factor(dat$set_type)))
dat$lga_f <- factor(as.numeric(as.factor(dat$GEONAME)))
dat$state_f <- factor(as.numeric(as.factor(dat$state_name)))
class(dat$pop <- as.integer(dat$pop))
sum(dat$pop, na.rm=T)

# Prepare covariate for INLA stack
covars <- dat[, c(covs_dens,"set_typ","state_f","ward_f", "lga_f", "eps")]

# Raw data stack
dat$dens[is.infinite(dat$dens)] = NA
dat$dens[is.nan(dat$dens)] = NA
dat$dens[dat$dens==0] = NA
stk_nmep <- inla.stack(data=list(y=dat$dens), #the response
                       A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                       effects=list(c(list(Intercept=1), #the Intercept
                                      iset),
                                    list(covars)
                       ),
                       #this is a quick name so you can call upon easily
                       tag='est')


# Imputed Wards stack
dat$dens_imp[is.infinite(dat$dens_imp)] = NA
dat$dens_imp[is.nan(dat$dens_imp)] = NA
dat$dens_imp[dat$dens_imp==0] = NA
stk_nmep_imp <- inla.stack(data=list(y=dat$dens_imp), #the response
                           A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                           effects=list(c(list(Intercept=1), #the Intercept
                                          iset),
                                        list(covars)
                           ),
                           #this is a quick name so you can call upon easily
                           tag='est')

# Dropped wards stack
dat$dens_drp[is.infinite(dat$dens_drp)] = NA
dat$dens_drp[is.nan(dat$dens_drp)] = NA
dat$dens_drp[dat$dens_drp==0] = NA
stk_nmep_drp <- inla.stack(data=list(y=dat$dens_drp), #the response
                           A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                           effects=list(c(list(Intercept=1), #the Intercept
                                          iset),
                                        list(covars)
                           ),
                           #this is a quick name so you can call upon easily
                           tag='est')


# Specify the PC prior
pcprior <- list(prec = list(prior = "pc.prec", param = c(1, 0.01)))  # PC prior


#-----------------------
###  Base  -----
#--------------------------
# model 11  - with state, lga and wards random effects
f_nmep301 <-  y ~ -1 + Intercept +x28 + x30 + x45 + x48  + x53 + x54+#
  f(s, model=spde) + f(state_f, model="iid", hyper = pcprior) 
mod_nmep301 <- inla(f_nmep301,
                   family="gamma",
                   control.compute = list(dic=T, waic=T, cpo=T, config=T),
                   data = inla.stack.data(stk_nmep, spde=spde),
                   control.predictor = list(A = inla.stack.A(stk_nmep), compute=T))

# model 12 -- with state and wards random effects
f_nmep302 <-  y ~ -1 + Intercept +x28 + x30 + x45 + x48  + x53 + x54+#
   f(state_f, model="iid", hyper = pcprior)
mod_nmep302 <- inla(f_nmep302,
                   family="gamma",
                   control.compute = list(dic=T, waic=T, cpo=T, config=T),
                   data = inla.stack.data(stk_nmep, spde=spde),
                   control.predictor = list(A = inla.stack.A(stk_nmep), compute=T))


# model 13  with wars random effects only
f_nmep30 <-  y ~ -1 + Intercept +x28 + x30 + x45 + x48  + x53 + x54+#
  f(s, model=spde) + f(eps, model="iid", hyper = pcprior)
mod_nmep30 <- inla(f_nmep30,
                   family="gamma",
                   control.compute = list(dic=T, waic=T, cpo=T, config=T),
                   data = inla.stack.data(stk_nmep, spde=spde),
                   control.predictor = list(A = inla.stack.A(stk_nmep), compute=T))

# Model fit checks
(DIC <- t(c( mod301 = mod_nmep301$dic$dic, # for raw data
             mod302 = mod_nmep302$dic$dic,
             mod303 = mod_nmep30$dic$dic)))
#    mod301  mod302   mod303
#[1,] 9992.497 10566.4 9924.652


summary(mod_nmep30)
index_nmep30 <-inla.stack.index(stk_nmep, "est")$data #---extract the data location indices
pred_nmep30<- exp(mod_nmep30$summary.linear.predictor[index_nmep30 ,"mean"]) #--predicted mean count
pred_nmep30L<- exp(mod_nmep30$summary.linear.predictor[index_nmep30 ,"0.025quant"]) #--lower bound
pred_nmep30U<- exp(mod_nmep30$summary.linear.predictor[index_nmep30 ,"0.975quant"]) #--Upper bound 
pred_nmep30Uncert<- (pred_nmep30U - pred_nmep30L)/pred_nmep30 #-Uncertainty


sum(fit_nmep30 <- pred_nmep30*dat$bld_count)
fit_nmep30L <- pred_nmep30L*dat$bld_count
fit_nmep30U <- pred_nmep30U*dat$bld_count
fit_nmep30Uncert<- (fit_nmep30U - fit_nmep30L)/fit_nmep30 #-Uncertainty
boxplot(pred_nmep30)
#plot(dat$dens, pred_nmep30)
#abline(lm(pred_nmep30 ~ dat$dens), col=2)

boxplot(fit_nmep30)
#-----------------------
# Imputed  -----
#------------------------

# model 21
f_nmep30b1 <-  y ~ -1 + Intercept +x28 + x30 + x45 + x48  + x53 + x54+#
  f(s, model=spde) + f(state_f, model="iid", hyper = pcprior) 
mod_nmep30b1 <- inla(f_nmep30b1,
                    family="gamma",
                    control.compute = list(dic=T, waic=T, cpo=T, config=T),
                    data = inla.stack.data(stk_nmep_imp, spde=spde),
                    control.predictor = list(A = inla.stack.A(stk_nmep_imp), compute=T))

# model 22
f_nmep30b2 <-   y ~ -1 + Intercept +x28 + x30 + x45 + x48  + x53 + x54+#
  f(state_f, model="iid", hyper = pcprior)
mod_nmep30b2 <- inla(f_nmep30b2,
                    family="gamma",
                    control.compute = list(dic=T, waic=T, cpo=T, config=T),
                    data = inla.stack.data(stk_nmep_imp, spde=spde),
                    control.predictor = list(A = inla.stack.A(stk_nmep_imp), compute=T))

# model 23
f_nmep30b <-  y ~ -1 + Intercept +x28 + x30 + x45 + x48  + x53 + x54+#
  f(s, model=spde) + f(eps, model="iid", hyper = pcprior)
mod_nmep30b <- inla(f_nmep30b,
                    family="gamma",
                    control.compute = list(dic=T, waic=T, cpo=T, config=T),
                    data = inla.stack.data(stk_nmep_imp, spde=spde),
                    control.predictor = list(A = inla.stack.A(stk_nmep_imp), compute=T))


# Model fit checks
(DIC <- t(c( mod301b = mod_nmep30b1$dic$dic, # for raw data
             mod302b = mod_nmep30b2$dic$dic,
             mod303b = mod_nmep30b$dic$dic)))
#mod301b  mod302b  mod303b
#[1,] 9287.32 9924.142 9287.088


summary(mod_nmep30b)
index_nmep30b <-inla.stack.index(stk_nmep_imp, "est")$data #---extract the data location indices
pred_nmep30b<- exp(mod_nmep30b$summary.linear.predictor[index_nmep30b ,"mean"]) #--predicted mean count
pred_nmep30bL<- exp(mod_nmep30b$summary.linear.predictor[index_nmep30b ,"0.025quant"]) #--lower bound
pred_nmep30bU<- exp(mod_nmep30b$summary.linear.predictor[index_nmep30b ,"0.975quant"]) #--Upper bound 
pred_nmep30bUncert<- (pred_nmep30bU - pred_nmep30bL)/pred_nmep30b #-Uncertainty


sum(fit_nmep30b <- pred_nmep30b*dat$bld_count)
fit_nmep30bL <- pred_nmep30bL*dat$bld_count
fit_nmep30bU <- pred_nmep30bU*dat$bld_count
fit_nmep30bUncert<- (fit_nmep30bU - fit_nmep30bL)/fit_nmep30b #-Uncertainty


boxplot(pred_nmep30b)
#plot(dat$dens_imp, pred_nmep30b)
#abline(lm(pred_nmep30b ~ dat$dens_imp), col=2)
boxplot(fit_nmep30b)


#-----------------------
# Dropped  -----
#------------------------
# model 31  
f_nmep301c <-  y ~ -1 + Intercept +x28 + x30 + x45 + x48  + x53 + x54+#
  f(s, model=spde) + f(state_f, model="iid", hyper = pcprior) 
mod_nmep301c <- inla(f_nmep301c,
                    family="gamma",
                    control.compute = list(dic=T, waic=T, cpo=T, config=T),
                    data = inla.stack.data(stk_nmep_drp, spde=spde),
                    control.predictor = list(A = inla.stack.A(stk_nmep_drp), compute=T))

# model 32 -- with state and wards random effects
f_nmep302c <-   y ~ -1 + Intercept +x28 + x30 + x45 + x48  + x53 + x54+#
  f(state_f, model="iid", hyper = pcprior)
mod_nmep302c <- inla(f_nmep302c,
                    family="gamma",
                    control.compute = list(dic=T, waic=T, cpo=T, config=T),
                    data = inla.stack.data(stk_nmep_drp, spde=spde),
                    control.predictor = list(A = inla.stack.A(stk_nmep_drp), compute=T))


# model 33  with wars random effects only
f_nmep30c <- y ~ -1 + Intercept +x28 + x30 + x45 + x48  + x53 + x54+#
  f(s, model=spde) + f(eps, model="iid", hyper = pcprior)
mod_nmep30c <- inla(f_nmep30c,
                   family="gamma",
                   control.compute = list(dic=T, waic=T, cpo=T, config=T),
                   data = inla.stack.data(stk_nmep_drp, spde=spde),
                   control.predictor = list(A = inla.stack.A(stk_nmep_drp), compute=T))

# Model fit checks
(DICc<- t(c( mod301c = mod_nmep301c$dic$dic, # for raw data
             mod302c = mod_nmep302c$dic$dic,
             mod303c = mod_nmep30c$dic$dic)))
#mod301c  mod302c  mod303c
#[1,] 6874.722 7493.952 2590.333




f_nmep30c <-  y ~ -1 + Intercept +x28 + x30 + x45 + x48  + x53 + x54+
  f(s, model=spde) +
  f(eps, model="iid", hyper = pcprior)
mod_nmep30c <- inla(f_nmep30c,
                    family="gamma",
                    control.compute = list(dic=T, waic=T, cpo=T, config=T),
                    data = inla.stack.data(stk_nmep_drp, spde=spde),
                    control.predictor = list(A = inla.stack.A(stk_nmep_drp), compute=T))

summary(mod_nmep30c)
index_nmep30c <-inla.stack.index(stk_nmep_drp, "est")$data #---extract the data location indices
pred_nmep30c<- exp(mod_nmep30c$summary.linear.predictor[index_nmep30c ,"mean"]) #--predicted mean count
pred_nmep30cL<- exp(mod_nmep30c$summary.linear.predictor[index_nmep30c ,"0.025quant"]) #--lower bound
pred_nmep30cU<- exp(mod_nmep30c$summary.linear.predictor[index_nmep30c ,"0.975quant"]) #--Upper bound 
pred_nmep30cUncert<- (pred_nmep30cU - pred_nmep30cL)/pred_nmep30c #-Uncertainty


sum(fit_nmep30c <- pred_nmep30c*dat$bld_count)
fit_nmep30cL <- pred_nmep30cL*dat$bld_count
fit_nmep30cU <- pred_nmep30cU*dat$bld_count
fit_nmep30cUncert<- (fit_nmep30cU - fit_nmep30cL)/fit_nmep30c #-Uncertainty

boxplot(pred_nmep30c)

boxplot(fit_nmep30c)
# Initial model fit checks
(DIC <- t(c( mod30 = mod_nmep30$dic$dic, # for raw data
             mod30b = mod_nmep30b$dic$dic,
             mod30c = mod_nmep30c$dic$dic)))

#        mod30   mod30b   mod30c
#[1,] 9826.869 9155.738 4898.707
# Add to the datasets
# raw
dat$pred30 <- fit_nmep30
dat$pred30uncert <- fit_nmep30Uncert
sum(dat$pred30uncert, na.rm=T)

# imputed
dat$pred30b <- fit_nmep30b
dat$pred30b_uncert <- fit_nmep30bUncert
sum(dat$pred30b_uncert, na.rm=T)

# dropped data
dat$pred30c <- fit_nmep30c
dat$pred30c_uncert <- fit_nmep30cUncert
sum(dat$pred30c_uncert, na.rm=T)

# summarise by state
data.frame(state_total <- dat %>% group_by(state_name) %>%
             summarise(obs = round(sum(pop, na.rm=T)),
                       pred30a = round(sum(pred30, na.rm=T)),
                       pred30b = round(sum(pred30b, na.rm=T)),
                       pred30c = round(sum(pred30c, na.rm=T))))

apply(state_total[,2:5], 2, sum, na.rm=T)

#  Model Validation
mod_metrics2 <- function(obs, pred)
{
  residual = pred - obs
  INACCURACY = mean(abs(residual), na.rm=T)#MAE
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])
  
  output <- list(MAE  = INACCURACY ,
                 RMSE = RMSE,
                 #BIAS = abs(BIAS),
                 corr = corr)
  return(output)
}


####   Model fit checks
met3a <- unlist(mod_metrics2(dat$pop, dat$pred30))
met3b <- unlist(mod_metrics2(dat$pop_imp, dat$pred30b))
met3c <- unlist(mod_metrics2(dat$pop_drp, dat$pred30c))
rbind(met3a, met3b, met3c)


#          MAE     RMSE      corr
#met3a 8660.184 19267.36 0.8086793
#met3b 8361.501 18663.77 0.8055178
#met3c 5212.302  9701.98 0.9131219


###----------------------------------------------------------------------------------
#    Model cross validation
####---------------------------------------------------------------------------------

#----Extract settement type effects
set_t <- function(dat, st)
{
  uniq <- unique(dat$set_typ)
  uniq[1]
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:3)
    {
      if(dat$set_typ[i]==uniq[j]) dat$set_typ2[i] = st[j]
    }
    
  }
  dat$set_typ2
}



#'-----------------------------
# Make scatter plots
#------------------------------
# Base
pred_dt1 <- cross_val1$pred_dat[cross_val1$pred_dat$data == "outsample",]
pred_dt1$method <- rep("Base", nrow(pred_dt1))
# Imputed
pred_dt2 <- cross_val2$pred_dat[cross_val2$pred_dat$data == "outsample",]
pred_dt2$method <- rep("Imputed", nrow(pred_dt2))
# Dropped
pred_dt3 <- cross_val3$pred_dat[cross_val3$pred_dat$data == "outsample",]
pred_dt3$method <- rep("Dropped", nrow(pred_dt3))

# Join all
pred_dt <- rbind(pred_dt1, pred_dt2, pred_dt3)
pred_dt$method <- factor(pred_dt$method,
                         levels = c("Base", "Imputed", "Dropped"))

library(ggpubr)

pred_dtb <- pred_dt %>% filter(pred < 250000)
plot_cval <-ggplot(pred_dtb, aes(x=obs, y=pred))+
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw()+
  theme(strip.text = element_text(size=20))+
  facet_wrap(~method)

rcval <-  ggpar(plot_cval, xlab="Observed counts", ylab="Predicted Counts",
                legend = "top",
                legend.title = "Fold (k{=5}-fold)",size=18,
                font.legend=c(18),
                font.label = list(size = 18, face = "bold", color ="red"),
                font.x = c(18),
                font.y = c(18),
                font.main=c(18),
                font.xtickslab =c(18),
                font.ytickslab =c(18),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rcval

#--------------------------------------------------------------------
## *******************Run Grid Cell Predictions******************
pred_covs2 <-  readRDS(paste0(cov_path, # read on the prediction stack
                              "/pred_covs/NGA_NMEP_predict_cov4.rds"))
dim(pred_covs2 <- as.data.frame(pred_covs2))
names(pred_covs2)


# select variables to include to minimize memory use
vars2put <- c("gridcell_id","x", "y", "bld_count", "wardname","set_typ",
              "GEONAME3","state_name",covs_dens)
pred_covs2 <- pred_covs2[, vars2put]

names(pred_covs2)
# checks for NAs
sapply(pred_covs2[,covs_dens], function(x) length(x[is.na(x)]))

# drop grid cells with no state names
dim(pred_covs2 <- pred_covs2 %>% drop_na(state_name))

# apply scaling
pred_covs2[,covs_dens] <- apply(pred_covs2[,covs_dens],2, stdize2)


# summarise covariates
apply(pred_covs2[,covs_dens], 2, summary)


# Build prediction projection matrix
Apred <- inla.spde.make.A(mesh = mesh,
                          loc = cbind(pred_covs2$x, pred_covs2$y));dim(Apred)


dim(Apred)

# Checks
table(pred_covs2$set_typ <- as.factor(as.numeric(pred_covs2$set_typ)))
length(pred_covs2$set_typ[is.na(pred_covs2$set_typ)])

#pred_covs2$set_typ[is.na(pred_covs2$set_typ)] =  0 # this is dummy


# function for extraction of random effects
str_ranef <- function(dat, strat, st)
{
  uniq <- unique(strat[strat!=0])
  ranef <- rep(0, length(uniq))
  for(i in  1:nrow(dat))
  {
    
    for(j in 1:length(uniq))
    {
      if(strat[i]==uniq[j]) ranef[i] = st[j]
    }
    
  }
  ranef
}


# Set all standardized covariates with NA values to zero
pred_covs2$x28[is.na(pred_covs2$x28)] = 0
pred_covs2$x30[is.na(pred_covs2$x30)] = 0
pred_covs2$x37[is.na(pred_covs2$x37)] = 0
pred_covs2$x43[is.na(pred_covs2$x43)] = 0
pred_covs2$x45[is.na(pred_covs2$x45)] = 0
pred_covs2$x46[is.na(pred_covs2$x46)] = 0
pred_covs2$x47[is.na(pred_covs2$x47)] = 0
pred_covs2$x48[is.na(pred_covs2$x48)] = 0
pred_covs2$x53[is.na(pred_covs2$x53)] = 0
pred_covs2$x54[is.na(pred_covs2$x54)] = 0


# Run posterior simulation and grid prediction simultaneously
simCom <- function(model1, dat, Aprediction1,run)
{
  fixedeff1  <- dens_hat <- pop_hat  <-  matrix(0,nrow=nrow(dat), ncol = run)
  #inla.seed = as.integer(runif(1)*.Machine$integer.max)
  inla.seed =   1912025672
  set.seed(inla.seed)
  print(inla.seed)
  
  # Obtain posterior samples of the Betas
  m1.samp1 <- inla.posterior.sample(run,
                                    model1, seed = inla.seed ,
                                    selection=list(x28=1,
                                                   x30=1,
                                                   x45=1,
                                                   x48=1,
                                                   x53=1,
                                                   x54=1
                                    ),
                                    num.threads="1:1")

  #Extract the spatial random effects at the mesh nodes
  sfield_nodes_mean1 <- model1$summary.random$s['mean']
  field_mean1 <- (Aprediction1%*% as.data.frame(sfield_nodes_mean1)[, 1])
  
  for(i in 1:run)
  {
    fixedeff1[,i] <-
      model1$summary.fixed['Intercept', 'mean'] +
      m1.samp1[[i]]$latent[1,] * dat[,'x28'] +
      m1.samp1[[i]]$latent[2,] * dat[,'x30'] +
      m1.samp1[[i]]$latent[3,] * dat[,'x45'] +
      m1.samp1[[i]]$latent[4,] * dat[,'x48'] +
      m1.samp1[[i]]$latent[5,] * dat[,'x53'] +
      m1.samp1[[i]]$latent[6,] * dat[,'x54'] +
      
      rnorm(nrow(dat), 0, 1/model1$summary.hyperpar$mean[4]) + # ward
      field_mean1[,1]
    
    dens_hat[,i]<- exp(fixedeff1[,i]) # population density
    pop_hat[,i] <- dens_hat[,i]*dat$bld_count # population count
  }
  
  dat$mean_pop_hat <- apply(pop_hat, 1, mean, na.rm=T) #
  dat$lower_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.025), na.rm=T) #
  dat$upper_pop_hat <- apply(pop_hat, 1, quantile, probs=c(0.975), na.rm=T) #
  dat$sd_pop_hat <- apply(pop_hat, 1, sd, na.rm=T) #
  dat$cv_pop_hat <- dat$sd_pop_hat/dat$mean_pop_hat
  
  output <- list(pop_hat = pop_hat,
                 est_data = dat)
}

run = 30

#
dim(pred_covs2 <- as.data.frame(pred_covs2))
system.time(str(sim.nmep <-  simCom(mod_nmep30,pred_covs2, Apred,run))) #run the draws and predict takes about 34 minutes
sum(sim.nmep$est_data$mean_pop_hat, na.rm=T)

system.time(str(sim.nmepb <-  simCom(mod_nmep30b,pred_covs2, Apred,run))) #run the draws and predict takes about 34 minutes
sum(sim.nmepb$est_data$mean_pop_hat, na.rm=T)

system.time(str(sim.nmepc <-  simCom(mod_nmep30c,pred_covs2, Apred,run))) #run the draws and predict takes about 34 minutes
sum(sim.nmepc$est_data$mean_pop_hat, na.rm=T)


min(sim.nmepc$est_data$mean_pop_hat, na.rm=T)
max(sim.nmepc$est_data$mean_pop_hat, na.rm=T)
max(sim.nmep$est_data$mean_pop_hat, na.rm=T)
max(sim.nmepb$est_data$mean_pop_hat, na.rm=T)


draw <- 1 :run # draws


# extract posterior estimates for zonal statistcs calculation
dim(data.nmep<- data.frame(cbind(pred_covs2[, c("wardname","GEONAME3", "state_name")], sim.nmep$pop_hat)))
names(data.nmep)
dim(data.nmepb<- data.frame(cbind(pred_covs2[, c("wardname","GEONAME3", "state_name")], sim.nmepb$pop_hat)))
names(data.nmepb)
dim(data.nmepc<- data.frame(cbind(pred_covs2[, c("wardname","GEONAME3", "state_name")], sim.nmepc$pop_hat)))
names(data.nmepc)


# National total
nat_total<- function(dat, run)
{
  p_hat <- dat[,(draw+3)]
  tots <- apply(p_hat,2, sum, na.rm=T)
  tot_sd  <- sd(tots, na.rm=T)
  
  tot_mean  <- mean(tots, na.rm=T)
  
  tot_lower <- quantile(tots, probs=c(0.025))
  
  tot_median <- quantile(tots, probs=c(0.5))
  tot_upper <- quantile(tots, probs=c(0.975))
  
  return(estimates <- data.frame(estimates=unlist(list(total=tot_mean,
                                                       lower=tot_lower, median=tot_median, upper=tot_upper))))
}


# Base
(nat_nmep <- nat_total(data.nmep,run))
(nat_nmep <- data.frame(total= nat_nmep[1,],
                        lower = nat_nmep[2,],
                        median=nat_nmep[3,],
                        upper=nat_nmep[4,]))

# Imputed
(nat_nmepb <- nat_total(data.nmepb,run))
(nat_nmepb <- data.frame(total= nat_nmepb[1,],
                         lower = nat_nmepb[2,],
                         median=nat_nmepb[3,],
                         upper=nat_nmepb[4,]))

# Dropped
(nat_nmepc <- nat_total(data.nmepc,run))
(nat_nmepc <- data.frame(total= nat_nmepc[1,],
                         lower = nat_nmepc[2,],
                         median=nat_nmepc[3,],
                         upper=nat_nmepc[4,]))

(nat_est <- data.frame(rbind(nat_nmep, nat_nmepb, nat_nmepc)))
nat_est$method <- c("Base", "Imputed", "Dropped")
nat_est
# save
write.csv(nat_est, paste0(results_path, "/national_estimates.csv"), row.names = F)


##---State level estimates---------------------------------
state_est <- function(datr,run)
{
  uniR <-as.numeric(as.factor(unique(datr$state_name)))
  names <-unique(datr$state_name)
  outR <- matrix(0, nrow=length(uniR), ncol=5)
  for(j in 1:length(uniR))
  {
    reg <- datr[datr$state_name==names[j],]
    rtots <- apply(reg[,(draw+3)], 2, sum, na.rm=T)
    rtot_mean  <- mean(rtots, na.rm=T)
    rtot_sd <- sd(rtots, na.rm=T)
    rtot_lower <- quantile(rtots, probs=c(0.025))
    rtot_median <- quantile(rtots, probs=c(0.5))
    rtot_upper <- quantile(rtots, probs=c(0.975))
    rtot_uncert <- (rtot_upper - rtot_lower)/rtot_mean
    
    restimates <- round(c(rtot_mean, rtot_lower, rtot_median,rtot_upper, rtot_uncert),4)
    outR[j,] <- restimates
  }
  outR <- data.frame(outR)
  return(reg_est <- data.frame(ID = uniR,
                               state = names,
                               total = outR[,1],
                               lower = outR[,2],
                               median = outR[,3],
                               upper = outR[,4],
                               uncertainty = outR[,5]))
}
#
# Base
(state.nmep <- state_est(data.nmep, run))
(state.nmep <- state.nmep[order(state.nmep$state),])

# Imputed
(state.nmepb <- state_est(data.nmepb, run))
(state.nmepb <- state.nmepb[order(state.nmepb$state),])

# Dropped
(state.nmepc <- state_est(data.nmepc, run))
(state.nmepc <- state.nmepc[order(state.nmepc$state),])
write.csv(state.nmepc, paste0(results_path, "/state_estimates.csv"), row.names = F)

#--------------------------------------------------------------
## lga estimates -------------------------------------------

lga_est <- function(datr, run)
{
  uniR <-as.numeric(as.factor(unique(datr$lga_name)))
  names <-unique(datr$lga_name)
  outR <- matrix(0, nrow=length(uniR), ncol=5)
  state <- rep(0, length(uniR))
  for(j in 1:length(uniR))
  {
    reg <- datr[datr$lga_name==names[j],]
    state[j] <- unique(datr$state_name[datr$lga_name==names[j]])
    rtots <- apply(reg[,(draw+3)], 2, sum, na.rm=T)
    rtot_mean  <- mean(rtots, na.rm=T)
    rtot_sd <- sd(rtots, na.rm=T)
    rtot_lower <- quantile(rtots, probs=c(0.025))
    rtot_median <- quantile(rtots, probs=c(0.5))
    rtot_upper <- quantile(rtots, probs=c(0.975))
    rtot_cv<- rtot_sd/rtot_mean
    
    restimates <- round(c(rtot_mean, rtot_lower, rtot_median,rtot_upper, rtot_cv),4)
    outR[j,] <- restimates
  }
  outR <- data.frame(outR)
  return(reg_est <- data.frame(ID = uniR,
                               lga = names,
                               state = state,
                               total = outR[,1],
                               lower = outR[,2],
                               median = outR[,3],
                               upper = outR[,4],
                               cv = outR[,5]))
}

#data.nmep$lga_name <- data.nmep$GEONAME3
#(lga.nmep <- lga_est(data.nmep, run))


# Base
data.nmep$lga_name <- data.nmep$GEONAME3
(lga.nmep <- lga_est(data.nmep, run))
(lga.nmep <- lga.nmep[order(lga.nmep$ID),])

# Imputed
data.nmepb$lga_name <- data.nmepb$GEONAME3
(lga.nmepb <- lga_est(data.nmepb, run))
(lga.nmepb <- lga.nmepb[order(lga.nmepb$ID),])

# Dropped
data.nmepc$lga_name <- data.nmepc$GEONAME3
(lga.nmepc <- lga_est(data.nmepc, run))
(lga.nmepc <- lga.nmepc[order(lga.nmepc$ID),])
write.csv(lga.nmepc, paste0(results_path, "/lga_estimates.csv"), row.names = F)
#--------------------------------------------------------------
## WARD estimates -------------------------------------------

wd_est <- function(datr, run)
{
  
 # datr <- data.nmepc
  uniR <- as.numeric(as.factor(paste(datr$wardname, 
                                     datr$GEONAME3,  # lga names
                                     datr$state_name)))
  
  datr$uniR <- uniR
  unique_rows <- datr %>% distinct(uniR, .keep_all = TRUE)
  names <- unique_rows$wardname
  outR <- matrix(0, nrow=length(unique(uniR)), ncol=5)
  state <- rep(0, length(unique(uniR)))
  lga <- rep(0, length(unique(uniR)))
  for(j in 1:length(unique(uniR)))
  {
    print(j)
    reg <- datr[datr$uniR==unique(uniR)[j],]
    lga[j] <- unique(datr$GEONAME3[datr$uniR==unique(uniR)[j]])
    state[j] <- unique(datr$state_name[datr$uniR==unique(uniR)[j]])
    rtots <- apply(reg[,(draw+3)], 2, sum, na.rm=T)
    rtot_mean  <- mean(rtots, na.rm=T)
    rtot_sd <- sd(rtots, na.rm=T)
    rtot_lower <- quantile(rtots, probs=c(0.025))
    rtot_median <- quantile(rtots, probs=c(0.5))
    rtot_upper <- quantile(rtots, probs=c(0.975))
    rtot_cv<- rtot_sd/rtot_mean
    
    restimates <- round(c(rtot_mean, rtot_lower, rtot_median,rtot_upper, rtot_cv),4)
    outR[j,] <- restimates
  }
  outR <- data.frame(outR)
  return(reg_est <- data.frame(ID = unique(uniR),
                               ward = names,
                               lga = lga,
                               state = state,
                               total = outR[,1],
                               lower = outR[,2],
                               median = outR[,3],
                               upper = outR[,4],
                               cv = outR[,5]))
}

#data.nmep$wd_name <- data.nmep$GEONAME3
#(wd.nmep <- wd_est(data.nmep, run))


# Base
data.nmep$wd_name <- data.nmep$wardname
(wd.nmep <- wd_est(data.nmep, run))
(wd.nmep <- wd.nmep[order(wd.nmep$ID),])

# Imputed
data.nmepb$wd_name <- data.nmepb$wardname
(wd.nmepb <- wd_est(data.nmepb, run))
(wd.nmepb <- wd.nmepb[order(wd.nmepb$ID),])

# Dropped

data.nmepc$wd_name <- data.nmepc$wardname
(wd.nmepc <- wd_est(data.nmepc, run))
(wd.nmepc <- wd.nmepc[order(wd.nmepc$ID),])


sum(wd.nmepc$total); sum(lga.nmepc$total); sum(state.nmepc$total)
write.csv(wd.nmepc, paste0(results_path, "/ward_estimates.csv"), row.names = F)
cross_validate <- function(dat, n.folds, mod, form, A, seed)
{
  #--------------------------------------------
  # dat: the input survey data containing the all the variables
  # n.folds: number of test (k) folds to use
  # mod: the best model of the full or reference data
  # A: the projection  matrix used in training the full data model
  # seed: a random sample seed to make results reproducible
  #--------------------------------------------
  seed = 13235
  set.seed(seed)
  #dat <- dat2 # survey data
  N <- nrow(dat)

  
  table(ind_train <- factor(sample(x = rep(1:n.folds,c(rep(926,9), 932)),  # Sample IDs for training data
                                   size = N))) 
  
  table(as.numeric(ind_train))
  dat$k_fold <- as.numeric(ind_train)
  coords = cbind(dat$lon, dat$lat)
  
  
  k_uniq <-sort(unique(dat$k_fold))
  
  
  #---------------------------------------------------------------
  #                   in-sample
  #---------------------------------------------------------------
  
  met_list_in <- list()
  pred_list_in <- list()
  for(i in 1:length(k_uniq))
  {
    
    print(paste0("in-sample cross-validation using fold ", i, sep=""))
    test_ind <- which(dat$k_fold==k_uniq[i])
    dim(test <- dat[test_ind, ]) #---test set for fold i
    
    
    train_coords <- coords
    test_coords <- coords[test_ind,]
    
    
    # spatial random effects based on the full data best model
    sfield_nodes_mean <- mod$summary.random$s['mean']
    field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
    
    
    ##--------
    fixed <-
      mod$summary.fixed['Intercept', 'mean'] +
      mod$summary.fixed['x28', 'mean'] * test[,'x28'] +
      mod$summary.fixed['x30', 'mean'] * test[,'x30'] +
      mod$summary.fixed['x45', 'mean'] * test[,'x45'] +
      mod$summary.fixed['x48', 'mean'] * test[,'x48'] +
      mod$summary.fixed['x53', 'mean'] * test[,'x53'] +
      mod$summary.fixed['x54', 'mean'] * test[,'x54'] +
      
      mod$summary.random$eps['mean'][test_ind,1] + #--uncorrelated spatial random effects
      
      field_mean[test_ind,1]
    
    dens_ht <- exp(fixed)
    sum(pop_ht <- dens_ht*test$bld)

    
    # visualise samples
   par(mfrow =c(1,1))
    plot(test$obs, pop_ht, xlab = "Observed",
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5)
    
    
    
    # calculate fit metrics
    met_in <- mod_metrics2(test$pop,
                           pop_ht)
    
    met_list_in[[i]]<- unlist(met_in)
    pred_list_in[[i]] <- data.frame(obs = test$obs, pred = pop_ht,
                                    fold = rep(i, length(test$obs)),
                                    data = rep("insample", length(test$obs)))
  }
  met_list_in_dat <- do.call(rbind,met_list_in)
  metrics_in <- apply(met_list_in_dat, 2, mean)
  pred_list_in_dat <- do.call(rbind,pred_list_in)
  
  #-----------------------------------------------------------
  #               out - of -sample
  #-----------------------------------------------------------
  met_list_out <- list()
  pred_list_out <- list()
  for(i in 1:length(k_uniq))
  {
    
    print(paste0("out-of-sample cross-validation using fold ", i, sep=""))
    train_ind <- which(dat$k_fold!=k_uniq[i])
    test_ind <- which(dat$k_fold==k_uniq[i])
    dim(train <- dat[train_ind, ])#---train set for fold i
    dim(test <- dat[test_ind, ]) #---test set for fold i
    
    
    train_coords <- coords[train_ind,]
    test_coords <- coords[test_ind,]
    
    
    ###---Create projection matrices for training and testing datasets
    Ae<-inla.spde.make.A(mesh=mesh,loc=as.matrix(train_coords));dim(Ae) #training
    
    #####------------------------
    covars_train <- train[,c("x28", "x30", "x45", "x48", "x53", "x54", "eps")]; dim(covars_train)
    
    
    stk_train <- inla.stack(data=list(y=train$dens), #the response
                            
                            A=list(Ae,1),  #the A matrix; the 1 is included to make the list(covariates)
                            
                            effects=list(c(list(Intercept=1), #the Intercept
                                           iset),  #the spatial index
                                         #the covariates
                                         list(covars_train)
                            ),
                            tag='train')
    
    
    ###---Rerun INLA for model test prediction
    model <-inla(form, #the formula
                 data=inla.stack.data(stk_train,spde=spde),  #the data stack
                 family= 'gamma',   #which family the data comes from
                 control.predictor=list(A=inla.stack.A(stk_train),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE) #can include verbose=TRUE to see the log of the model runs
    summary(model)
    
    
    # Extract Spatial random effects from the full data best model
    
    sfield_nodes_mean <- mod$summary.random$s['mean']
    field_mean <- (A%*% as.data.frame(sfield_nodes_mean)[, 1])
    
    
    ##--------
    fixed <-
      model$summary.fixed['Intercept', 'mean'] +
      model$summary.fixed['x28', 'mean'] * test[,'x28'] +
      model$summary.fixed['x30', 'mean'] * test[,'x30'] +
      model$summary.fixed['x45', 'mean'] * test[,'x45'] +
      model$summary.fixed['x48', 'mean'] * test[,'x48'] +
      model$summary.fixed['x53', 'mean'] * test[,'x53'] +
      model$summary.fixed['x54', 'mean'] * test[,'x54'] +
      
      rnorm(nrow(test), 0, 1/model$summary.hyperpar$mean[4]) +
      field_mean[test_ind,1]
    dens_ht <- exp(fixed)
    sum(pop_ht <- dens_ht*test$bld)
    
    
    # visualise samples
 
    par(mfrow =c(1,1))
    plot(test$obs, pop_ht, xlab = "Observed",
         ylab = "Predicted", col=c('blue','orange'),
         pch=c(16,16), cex.axis=1.5)
    abline(0,1)
    legend("topleft", c("Observed", "Predicted"), col=c("blue", "orange"), pch=c(16,16),
           bty="n", cex=1.5)
    
    
    
    # calculate fit metrics
    met_out <- mod_metrics2(test$obs,
                            pop_ht)
    
    met_list_out[[i]]<- unlist(met_out)
    pred_list_out[[i]] <- data.frame(obs = test$obs, pred = pop_ht,
                                     fold = rep(i, length(test$obs)),
                                     data = rep("outsample", length(test$obs)))
  }
  met_list_out_dat <- do.call(rbind,met_list_out)
  metrics_out <- apply(met_list_out_dat, 2, mean) # fit metrics
  
  pred_list_out_dat <- do.call(rbind,pred_list_out)# predictions
  
  cv_mets <- rbind(metrics_in, metrics_out)
  output <- list( met_list_in_dat = met_list_in_dat,
                  met_list_out_dat = met_list_out_dat,
                  pred_dat = rbind(pred_list_in_dat, pred_list_out_dat),
                  cv_metrics = rbind(metrics_in, metrics_out))
}



##  Raw (Base) -----------------------------
dat1 <- dat
dat1$pop <- dat1$pop
dat1$dens <- dat1$dens
dat1$bld <- dat1$bld_count
dat1$obs <- dat1$pop # observed household size
######


(cross_val1 <- cross_validate(dat1, n.folds = 10,
                              mod = mod_nmep30,
                              form = f_nmep30,
                              A,
                              seed = 135))

cross_val1$met_list_out_dat  # out-of-sample metrics per fold

val_dat1 <- data.frame(cross_val1$met_list_out_dat)
val_dat1$method <- rep("Base", nrow(val_dat1))
val_dat1$DIC <- rep(mod_nmep30$dic$dic, nrow(val_dat1))


dim(cross_val1$pred_dat1 <- cross_val1$pred_dat %>% filter(fold==1))
pred_out <- cross_val1$pred_dat %>% filter(data == "outsample")
pred_out %>% group_by(fold)%>%
  summarise(obs = sum(obs, na.rm=T),
            pred = sum(pred, na.rm=T))

# Imputed ----------------
dat2 <- dat
dat2$pop <- dat2$pop_imp
dat2$dens <- dat2$dens_imp
dat2$bld <- dat2$bld_count
dat2$obs <- dat2$pop_imp # observed household siz

(cross_val2 <- cross_validate(dat2, n.folds = 10,
                              mod = mod_nmep30b,
                              form = f_nmep30b,
                              A,
                              seed = 1357))


cross_val2$met_list_out_dat  # out-of-sample metrics per fold
apply(cross_val2$met_list_out_dat, 2, mean)
val_dat2 <- data.frame(cross_val2$met_list_out_dat)
val_dat2$method <- rep("Imputed", nrow(val_dat2))
val_dat2$DIC <- rep(mod_nmep30b$dic$dic, nrow(val_dat2))



# Dropped ---------------------------------------
dat3 <- dat
dat3$pop <- dat3$pop_drp
dat3$dens <- dat3$dens_drp
dat3$bld <- dat3$bld_count
dat3$obs <- dat3$pop_drp # observed household siz
(cross_val3 <- cross_validate(dat3, n.folds = 10,
                              mod = mod_nmep30c,
                              form = f_nmep30c,
                              A,
                              seed = 138))

cross_val3$met_list_out_dat  # out-of-sample metrics per fold
val_dat3 <- data.frame(cross_val3$met_list_out_dat)
val_dat3$method <- rep("Dropped", nrow(val_dat3))
val_dat3$DIC <- rep(mod_nmep30c$dic$dic, nrow(val_dat3))


## join all validation data
val_dat <- rbind(val_dat1, val_dat2, val_dat3)

write.csv(val_dat, paste0(results_path, "/validation_data.csv"))


rbind(cross_val1$met_list_out_dat,
      cross_val2$met_list_out_dat,
      cross_val3$met_list_out_dat)
# means
(mets_all <- rbind(apply(cross_val1$met_list_out_dat, 2, mean),
                   apply(cross_val2$met_list_out_dat, 2, mean),
                   apply(cross_val3$met_list_out_dat, 2, mean)))

write.csv(mets_all, paste0(results_path, "/cross_validated_metrics.csv"), row.names=F)
#      MAE     RMSE      corr
#[1,] 9245.941 17309.74 0.8014975
#[2,] 7890.979 14828.00 0.8300593
#[3,] 7557.046 10623.00 0.8458744

# Reductions in relative MAE
(1- (mets_all[2,1]/mets_all[1,1]))*100 # 14.66% for imputed
(1- (mets_all[3,1]/mets_all[1,1]))*100 # 18.27% for dropped

# Reductions in relative MAE
(1- (mets_all[2,2]/mets_all[1,2]))*100 # 14.34% for imputed
(1- (mets_all[3,2]/mets_all[1,2]))*100 # 27.08% for dropped


# extract fixed effects
fixeffs1 <- round(mod_nmep30$summary.fixed[,c(1:3,5)],4)
fixeffs2 <- round(mod_nmep30b$summary.fixed[,c(1:3,5)],4)
fixeffs3 <- round(mod_nmep30c$summary.fixed[,c(1:3,5)],4)

fixed_effects <- rbind(fixeffs1, fixeffs2, fixeffs3)
write.csv(fixed_effects, paste0(results_path, "/fixed_effects.csv"), row.names=F)



#-------------------------------------------------------------
# cross-validation results plots
#------------------------------------------------------------
# Boxplots for fit metrics
metrics_all <- val_dat
metrics_all$Method <- factor(metrics_all$method,
                             levels = c("Base", "Imputed", "Dropped"))
# MAE
plot_mae <- metrics_all %>%
  ggplot(aes(x=Method, y=MAE, color=Method)) +
 # geom_violin(width=1.2) +
  geom_boxplot(width=0.8, 
               alpha=0.8,
                notchwidth = 0.2) +
  
  theme_bw()+
  theme(strip.text = element_text(size = 20),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))


rmae <-  ggpar(plot_mae, xlab="Method", ylab="MAE",
               legend = "none",
               legend.title = "Method:",size=20,,
               font.legend=c(20),
               font.label = list(size = 20, face = "bold", color ="red"),
               palette = "lancet",
               font.y = c(20),
               font.x = c(20),
               font.main=c(20),
               font.xtickslab =c(18),
               font.ytickslab =c(18),
               ytickslab.rt = 45)
rmae


# RMSE
plot_rmse <- metrics_all %>%
  ggplot( aes(x=Method, y=RMSE, color=Method)) +
  #geom_violin(width=1.2) +
   geom_boxplot(width=0.8, 
   alpha=0.8,
    notchwidth = 0.2) +
  
  theme_bw()+
  theme(strip.text = element_text(size = 20),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))


rrmse <-  ggpar(plot_rmse, xlab="Method", ylab="RMSE",
               legend = "none",
               legend.title = "Method:",size=20,,
               font.legend=c(20),
               font.label = list(size = 15, face = "bold", color ="red"),
               palette = "lancet",
               font.y = c(20),
               font.x = c(20),
               font.main=c(20),
               font.xtickslab =c(18),
               font.ytickslab =c(18),
               ytickslab.rt = 45)
rrmse



# CC
plot_CC <- metrics_all %>%
  ggplot(aes(x=Method, y=corr, color=Method)) +
  #geom_violin(width=1.2) +
  geom_boxplot(width=0.8, 
               alpha=0.8,
               notchwidth = 0.2) +
  
  theme_bw()+
  theme(strip.text = element_text(size = 20),
        axis.text.x=element_text(size=15),
        axis.text.y=element_text(size=15),
        legend.title=element_text(size=15),
        legend.text=element_text(size=14))


rCC  <-  ggpar(plot_CC , xlab="Method", ylab="CC ",
                legend = "none",
                legend.title = "Method:",size=20,,
                font.legend=c(20),
                font.label = list(size = 15, face = "bold", color ="red"),
                palette = "lancet",
                font.y = c(20),
                font.x = c(20),
                font.main=c(20),
                font.xtickslab =c(18),
                font.ytickslab =c(18),
                ytickslab.rt = 45)
rCC 


ggarrange(rmae, rrmse, rCC,
          nrow=1, ncol=3)

#'-----------------------------
# Make scatter plots
#------------------------------
# Base
pred_dt1 <- cross_val1$pred_dat[cross_val1$pred_dat$data == "outsample",]
pred_dt1$method <- rep("Base", nrow(pred_dt1))
# Imputed
pred_dt2 <- cross_val2$pred_dat[cross_val2$pred_dat$data == "outsample",]
pred_dt2$method <- rep("Imputed", nrow(pred_dt2))
# Dropped
pred_dt3 <- cross_val3$pred_dat[cross_val3$pred_dat$data == "outsample",]
pred_dt3$method <- rep("Dropped", nrow(pred_dt3))

# Join all
pred_dt <- rbind(pred_dt1, pred_dt2, pred_dt3)
pred_dt$method <- factor(pred_dt$method,
                         levels = c("Base", "Imputed", "Dropped"))

library(ggpubr)

pred_dtb <- pred_dt #%>% filter(pred < 250000)
plot_cval <-ggplot(pred_dtb, aes(x=obs, y=pred))+
  geom_point() +
  geom_smooth(method="lm") +
  theme_bw()+
  theme(strip.text = element_text(size=20))+
  facet_wrap(~method)

rcval <-  ggpar(plot_cval, xlab="Observed counts", ylab="Predicted Counts",
                legend = "top",
                legend.title = "Fold (k{=5}-fold)",size=18,
                font.legend=c(18),
                font.label = list(size = 18, face = "bold", color ="red"),
                font.x = c(18),
                font.y = c(18),
                font.main=c(18),
                font.xtickslab =c(18),
                font.ytickslab =c(18),
                # orientation = "reverse",
                xtickslab.rt = 45, ytickslab.rt = 45)
rcval

#---------------------------------------------------------------------------------------
# Write raster files for the best fit outputs
#--------------------------------------------------------------------------------------
ref_coords <- cbind(pred_covs2$x, pred_covs2$y)
x <- as.matrix(ref_coords)

## Mean
z3a.nmep <- as.matrix(sim.nmepc$est_data$mean_pop_hat)
ng3a.nmep = rasterFromXYZ(cbind(x, z3a.nmep))

# Lower
z3a.nmepL <- as.matrix(sim.nmepc$est_data$lower_pop_hat)
ng3a.nmepL = rasterFromXYZ(cbind(x, z3a.nmepL))

# Upper
z3a.nmepU <- as.matrix(sim.nmepc$est_data$upper_pop_hat)
ng3a.nmepU = rasterFromXYZ(cbind(x, z3a.nmepU))


# Coefficient of variation
z3a.nmepCV <- as.matrix(sim.nmepc$est_data$cv_pop_hat)
ng3a.nmepCV = rasterFromXYZ(cbind(x, z3a.nmepCV))


writeRaster(ng3a.nmep, filename=paste0(results_path, "/gridded_nga_total.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

writeRaster(ng3a.nmepL, filename=paste0(results_path, "/gridded_nga_lower.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

writeRaster(ng3a.nmepU, filename=paste0(results_path, "/gridded_nga_upper.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))

writeRaster(ng3a.nmepCV, filename=paste0(results_path, "/gridded_nga_CV.tif"),
            overwrite=TRUE, options = c('COMPRESS' = 'LZW'))



### Make wd maps
# add wd data to wd shapefile


wd.nmep <- wd.nmep[order(wd.nmep$state, wd.nmep$lga, wd.nmep$ward),]
wd.nmepb <- wd.nmepb[order(wd.nmepb$state, wd.nmepb$lga, wd.nmepb$ward),]
wd.nmepc <- wd.nmepc[order(wd.nmepc$state, wd.nmepc$lga, wd.nmepc$ward),]


# Merge with the shapefile
wd_shp <- shp
wd_shp$ward <- wd_shp$wrdnm_x
wd_base <- merge(wd_shp, wd.nmep, by = "ward")
wd_imp <- merge(wd_shp, wd.nmepb, by = "ward")
wd_drp <- merge(wd_shp, wd.nmepc, by = "ward")
#--------------------------
summary(wd_base$pop)
summary(wd_base$total)
summary(wd_imp$pop_imp)
summary(wd_imp$total)
summary(wd_drp$pop_drp)
summary(wd_drp$total)
breaks = c(0, 1000, 5000, 10000, 50000, 100000, 500000, 1100000)
tmap_options(check.and.fix = TRUE)

# base
#  Observed counts 
wdbase_obs <- tm_shape(wd_base) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("pop",# color="white",
              palette=inferno(100),
              style = "cont",
              title="Population",
              legend.show=F,
              textNA = "No data",
              breaks = c(0, 1000, 5000, 10000, 50000, 100000, 500000, 1100000)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout( frame=F)


#  Predicted counts 
wdbase_prd <- tm_shape(wd_base) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("total",
              palette=inferno(100),
              style = "cont",
              title="Population",
              legend.show=F,
              textNA = "No data",
              breaks = c(0, 1000, 5000, 10000, 50000, 100000, 500000, 1100000)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout( frame=F)


# Uncertainty
summary(wd_base$cv)
wdbase_uncert <- tm_shape(wd_base) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("cv",
              #palette=magma(100),
              palette=plasma(250),
              #palette="Reds",
              style = "cont",
              title="Uncertainty",
              legend.show=F,
              breaks = c(0.01, 0.05, 0.1,0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout(frame=F)


## imputed
#  Observed counts 
wdimp_obs <- tm_shape(wd_imp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("pop_imp",
              palette=inferno(100),
              style = "cont",
              title="Population",
              legend.show=F,
              textNA = "No data",
              breaks = c(0, 1000, 5000, 10000, 50000, 100000, 500000, 1100000)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout( frame=F)


#  Predicted counts 
wdimp_prd <- tm_shape(wd_imp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("total",
              palette=inferno(100),
              style = "cont",
              title="Population",
              legend.show=F,
              textNA = "No data",
              breaks = c(0, 1000, 5000, 10000, 50000, 100000, 500000, 1100000)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout( frame=F)


# Uncertainty
summary(wd_imp$cv)
wdimp_uncert <- tm_shape(wd_imp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("cv",
              #palette = c("#33FF57", "steelblue", "red"),
              palette=plasma(250),
              #palette="Reds",
              style = "cont",
              title="Uncertainty",
              legend.show=F,
              breaks = c(0.01, 0.05, 0.1,0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout(frame=F)




# Dropped
#  Observed counts 
wddrp_obs <- tm_shape(wd_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("pop_drp",
              palette=inferno(100),
              style = "cont",
              title="Population",
              legend.show=F,
              textNA = "No data",
              breaks = c(0, 1000, 5000, 10000, 50000, 100000, 500000, 1100000)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout( frame=F)


#  Predicted counts 
wddrp_prd <- tm_shape(wd_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("total",
              palette=inferno(100),
              style = "cont",
              title="Population",
              legend.show=T,
              textNA = "No data",
              breaks = c(0, 1000, 5000, 10000, 50000, 100000, 500000, 1100000)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout( frame=F)



# Uncertainty
summary(wd_drp$cv)
wddrp_uncert <- tm_shape(wd_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("cv",
              #palette = c("#33FF57", "darkblue", "#FF5733"),
              palette=plasma(250),
              #palette="Reds",
              style = "cont",
              title="Uncertainty",
              legend.show=F,
              breaks = c(0.01, 0.05, 0.1,0.15, 0.2, 0.25, 0.3, 0.35, 0.4)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout(frame=F)


tmap_arrange(wdbase_obs, wdimp_obs, wddrp_obs,
             wdbase_prd, wdimp_prd, wddrp_prd, 
             wdbase_uncert, wdimp_uncert, wddrp_uncert,
             nrow = 3, ncol= 3)


### added
#  Predicted counts 
wd_drp$uncertainty <- round((wd_drp$upper - wd_drp$lower)/wd_drp$total,4)
summary(wd_drp$uncertainty)
wd_drp_tot <- tm_shape(wd_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("total",
              palette=viridis(100),
              style = "cont",
              title="Population",
              legend.show=T,
              textNA = "No data",
              breaks = c(0, 1000, 5000, 10000, 50000, 100000, 500000, 1100000)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout( frame=T)



# Uncertainty
summary(wd_drp$cv)
wd_drp_uncert <- tm_shape(wd_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("cv",
              #palette = c("#33FF57", "darkblue", "#FF5733"),
              palette="Oranges",
              #palette="Reds",
              style = "cont",
              title="Uncertainty",
              legend.show=F,
              breaks = c(0.01, 0.05, 0.1,0.2, 0.3, 0.4, 0.5, 0.8,1.3)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout(frame=T)


# Scatter plots with 95% credible inmtervals



# Create scatter plot with vertical confidence intervals
### For ward
wd_base$obs <-wd_base$pop
wd_base$Method <- factor(rep("Base", nrow(wd_base)))

wd_imp$obs <- wd_imp$pop_imp
wd_imp$Method <- factor(rep("Imputed", nrow(wd_imp)))

wd_drp$obs <- wd_drp$pop_drp
wd_drp$Method <- factor(rep("Dropped", nrow(wd_drp)))


plot_wdrp <- wd_drp %>% filter(!is.na(obs), total!=0, obs!=0) %>%
  ggplot(aes(x = obs,y = total)) +
  geom_point(size = 3, color = "black") +  # Scatter points
  geom_smooth(method = "lm", se=F)+
  #geom_line(data = lg_dat_lng, aes(x = obs,y = total), color = "darkblue",size=1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width =2,size=0.6, color = "magenta") +  # CI lines
  #geom_pointrange(aes(ymin = lower, ymax = upper), size=0.6, color = "black") +  # CI lines
  theme_bw()+
  theme(strip.text = element_text(size = 18))
pscat <- ggpubr::ggpar(scat2, xlab="Observed population", ylab="Predicted population",
                       legend = "none", legend.title = "",size=22,
                       palette = "lancet",
                       font.label = list(size = 18, face = "bold", color ="red"),
                       font.x = c(18),
                       font.y = c(18),
                       font.main=c(18),
                       font.xtickslab =c(18),
                       font.ytickslab =c(18),
                       xtickslab.rt = 45,
                       ytickslab.rt = 45
)
pscat

# Combine all datasets
wd_dat_lng <- rbind(wd_base, wd_imp, wd_drp)
wd_dat_lng$Method <- factor(wd_dat_lng$Method)


wd_dat_lng2 <- wd_dat_lng %>% filter(!is.na(obs)) %>% filter(obs !=0, total!=0)

scat2 <- ggplot(wd_dat_lng2, aes(x = obs,y = total)) +
  geom_point(size = 3, color = "black") +  # Scatter points
  geom_smooth(method = "lm", se=F)+
  #geom_line(data = lg_dat_lng, aes(x = obs,y = total), color = "darkblue",size=1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width =2,size=0.6, color = "magenta") +  # CI lines
  #geom_pointrange(aes(ymin = lower, ymax = upper), size=0.6, color = "black") +  # CI lines
  theme_bw()+
  theme(strip.text = element_text(size = 18))+
  facet_wrap(~Method)

pscat <- ggpubr::ggpar(scat2, xlab="Observed population", ylab="Predicted population",
                       legend = "none", legend.title = "",size=22,
                       font.legend=c(18),
                       palette = "lancet",
                       font.label = list(size = 18, face = "bold", color ="red"),
                       font.x = c(18),
                       font.y = c(18),
                       font.main=c(18),
                       font.xtickslab =c(18),
                       font.ytickslab =c(18),
                       xtickslab.rt = 45,
                       ytickslab.rt = 45
)
pscat



## grid aggregated ward data
#wd_base2 <- wd_base %>% filter(!is.na(obs)) %>% filter(obs !=0, total!=0)
#wd_imp2 <- wd_imp %>% filter(!is.na(obs)) %>% filter(obs !=0, total!=0)
#wd_drp2 <- wd_drp %>% filter(!is.na(obs)) %>% filter(obs !=0, total!=0)

#met_base<- unlist(mod_metrics2(wd_base2 $obs,
#                               wd_base2$total))
#met_imp<- unlist(mod_metrics2(wd_imp2$obs,
#                              wd_imp2$total))
#met_drp<- unlist(mod_metrics2(wd_drp2$obs,
#                              wd_drp2$total))
#rbind(met_base, met_imp, met_drp)



### Make LGA maps
# add lga data to lga shapefile


lga.nmep <- lga.nmep[order(lga.nmep$state, lga.nmep$lga),]
lga.nmepb <- lga.nmepb[order(lga.nmepb$state, lga.nmepb$lga),]
lga.nmepc <- lga.nmepc[order(lga.nmepc$state, lga.nmepc$lga),]



lga_shp <- lga_shp[order(lga_shp$ID_lga),]
lga_base <- cbind(lga_shp, lga.nmep)
lga_imp <- cbind(lga_shp, lga.nmepb)
lga_drp <- cbind(lga_shp, lga.nmepc)
#--------------------------
summary(lga_base$Base)
summary(lga_base$total)
summary(lga_imp$Imputed)
summary(lga_imp$total)
summary(lga_drp$Dropped)
summary(lga_drp$total)
breaks = c(19000,40000, 70000,0.15e6,0.3e6, 0.6e6, 1e6, 2.5e6,4.3e6)
tmap_options(check.and.fix = TRUE)

# base
#  Observed counts 
lgbase_obs <- tm_shape(lga_base) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("Base",
              palette=magma(100),
              title="Population",
              legend.show=F,
              breaks = c(19000,40000, 70000,0.15e6,0.3e6, 0.6e6, 1e6, 2.5e6,4.3e6)
  )+
  tm_legend(outside = TRUE)+
  tm_layout(frame=F)

#  Predicted counts 
lgbase_prd <- tm_shape(lga_base) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("total",
              palette=magma(100),
              style= "cont",
              title="Population",
              legend.show=T,
              breaks = c(19000,40000, 70000,0.15e6,0.3e6, 0.6e6, 1e6, 2.5e6,4.3e6)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout(frame=F)


# Uncertainty
summary(lga_base$cv)
lgbase_uncert <- tm_shape(lga_base) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("cv",
              palette=magma(100),
              style="cont",
              title="Population",
              legend.show=F,
              breaks = c(0.01,0.025, 0.05, 0.1, 0.15, 0.2,0.25, 0.3,0.35,0.4)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout(frame=F)


## imputed
#  Observed counts 
lgimp_obs <- tm_shape(lga_imp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("Imputed",
              palette=magma(100),
              title="Population",
              legend.show=F,
              breaks = c(19000,40000, 70000,0.15e6,0.3e6, 0.6e6, 1e6, 2.5e6,4.3e6)
  )+
  tm_legend(outside = TRUE)+
  tm_layout(frame=F)

#  Predicted counts 
lgimp_prd <- tm_shape(lga_imp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("total",
              palette=magma(100),
              title="Population",
              legend.show=F,
              breaks = c(19000,40000, 70000,0.15e6,0.3e6, 0.6e6, 1e6, 2.5e6,4.3e6)
  )+
  tm_legend(outside = TRUE)+
  tm_layout(frame=F)


# Uncertainty
summary(lga_imp$cv)
lgimp_uncert <- tm_shape(lga_imp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("cv",
              palette=viridis(100),
              title="Population",
              legend.show=F,
              breaks = c(0.01,0.025, 0.05, 0.1, 0.15, 0.2,0.25, 0.3,0.35,0.4)
  )+
  tm_legend(outside = TRUE)+
  tm_layout(frame=F)




# Dropped
#  Observed counts 
lgdrp_obs <- tm_shape(lga_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("Dropped",
              palette=magma(100),
              title="Population",
              legend.show=F,
              breaks = c(19000,40000, 70000,0.15e6,0.3e6, 0.6e6, 1e6, 2.5e6,4.3e6)
  )+
  tm_legend(outside = TRUE)+
  tm_layout(frame=F)

#  Predicted counts 
lgdrp_prd <- tm_shape(lga_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("total",
              palette=viridis(100),
              style= "cont",
              title="Population",
              legend.show=F,
              breaks = c(19000,40000, 70000,0.15e6,0.3e6, 0.6e6, 1e6, 2.5e6,4.3e6)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout(frame=T)


# Uncertainty
summary(lga_drp$cv)
lgdrp_uncert <- tm_shape(lga_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("uncertainty",
              #palette=viridis(100),
              palette="Oranges",
              style = "cont",
              title="Population",
              legend.show=F,
              breaks = c(0.01, 0.05, 0.1,0.2, 0.3, 0.4, 0.5, 0.8,1.3)
  )+
  tm_legend(outside = TRUE, text.size=1.2)+
  tm_layout(frame=T)


tmap_arrange(lgdrp_prd, 
             lgdrp_uncert,
             nrow = 2, ncol= 1)

tmap_arrange(lgbase_obs, lgimp_obs, lgdrp_obs,
             lgbase_prd, lgimp_prd, lgdrp_prd, 
             lgbase_uncert, lgimp_uncert, lgdrp_uncert,
             nrow = 3, ncol= 3)



# LGA maps by states 
#  Predicted counts 
lga_drp$uncertainty <- round((lga_drp$upper - lga_drp$lower)/lga_drp$total,4)
summary(lga_drp$uncertainty)
lgdrp_prd_all <- tm_shape(lga_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("total",
              palette=viridis(100),
              style= "cont",
              title="pop",
              legend.show=T,
              breaks = c(19000,40000, 70000,0.15e6,0.3e6, 0.6e6, 1e6, 2.5e6,4.3e6)
  )+
  tm_legend(outside = TRUE, text.size=0.7)+
  tm_facets(by="state")+
  tm_layout(panel.label.size = 2.2) 


# Uncertainty
summary(lga_drp$cv)
lgdrp_uncert_all <- tm_shape(lga_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("uncertainty",
              #palette=viridis(100),
              palette="Oranges",
              style = "cont",
              title="Uncertainty",
              legend.show=F,
              breaks = c(0.01, 0.05, 0.1,0.2, 0.3, 0.4, 0.5, 0.8,1.3)
  )+
  tm_legend(outside = TRUE, text.size=0.7)+
  tm_facets(by="state")+
  tm_layout(panel.label.size = 2.2) 

# Scatter plots with 95% credible inmtervals



# Create scatter plot with vertical confidence intervals
### For state
lga_base$obs <-lga_base$Base
lga_base$Method <- factor(rep("Base", nrow(lga_base)))

lga_imp$obs <- lga_imp$Imputed
lga_imp$Method <- factor(rep("Imputed", nrow(lga_imp)))

lga_drp$obs <- lga_drp$Dropped
lga_drp$Method <- factor(rep("Dropped", nrow(lga_drp)))


plot_ldrp <- lga_drp %>% filter(total <1e+6) %>%
ggplot(aes(x = obs,y = total)) +
  geom_point(size = 3, color = "black") +  # Scatter points
  geom_smooth(method = "lm", se=F)+
  #geom_line(data = lg_dat_lng, aes(x = obs,y = total), color = "darkblue",size=1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width =2,size=0.6, color = "magenta") +  # CI lines
  #geom_pointrange(aes(ymin = lower, ymax = upper), size=0.6, color = "black") +  # CI lines
  theme_bw()+
  theme(strip.text = element_text(size = 18))

pldrp <- ggpubr::ggpar(plot_ldrp, xlab="Observed population", ylab="Predicted population",
                       legend = "none", legend.title = "",size=22,
                       font.legend=c(18),
                       palette = "lancet",
                       font.label = list(size = 18, face = "bold", color ="red"),
                       font.x = c(18),
                       font.y = c(18),
                       font.main=c(18),
                       font.xtickslab =c(18),
                       font.ytickslab =c(18),
                       xtickslab.rt = 45,
                       ytickslab.rt = 45
)
pldrp



# Combine all datasets
lg_dat_lng <- rbind(lga_base, lga_imp, lga_drp)
lg_dat_lng$Method <- factor(lg_dat_lng$Method)


lg_dat_lng <- lg_dat_lng %>% filter(!is.na(obs))
scat2 <- ggplot(lg_dat_lng, aes(x = obs,y = total)) +
  geom_point(size = 3, color = "black") +  # Scatter points
  geom_smooth(method = "lm", se=F)+
  #geom_line(data = lg_dat_lng, aes(x = obs,y = total), color = "darkblue",size=1.5) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width =2,size=0.6, color = "magenta") +  # CI lines
  #geom_pointrange(aes(ymin = lower, ymax = upper), size=0.6, color = "black") +  # CI lines
  theme_bw()+
  theme(strip.text = element_text(size = 18))+
  facet_wrap(~Method)

pscat <- ggpubr::ggpar(scat2, xlab="Observed population", ylab="Predicted population",
                       legend = "none", legend.title = "",size=22,
                       font.legend=c(18),
                       palette = "lancet",
                       font.label = list(size = 18, face = "bold", color ="red"),
                       font.x = c(18),
                       font.y = c(18),
                       font.main=c(18),
                       font.xtickslab =c(18),
                       font.ytickslab =c(18),
                       xtickslab.rt = 45,
                       ytickslab.rt = 45
)
pscat


lga_base <- lga_base %>% drop_na(obs)
lga_imp <- lga_imp %>% drop_na(obs)
lga_drp <- lga_drp %>% drop_na(obs)



### State maps
state.nmep <- state.nmep[order(state.nmep$state),]
state.nmepb <- state.nmepb[order(state.nmepb$state),]
state.nmepc <- state.nmepc[order(state.nmepc$state),]



state_shp <- state_shp[order(state_shp$stat_nm ),]
state_base <- cbind(state_shp, state.nmep)
state_imp <- cbind(state_shp, state.nmepb)
state_drp <- cbind(state_shp, state.nmepc)
#--------------------------
summary(state_base$Base)
summary(state_base$total)
summary(state_imp$Imputed)
summary(state_imp$total)
summary(state_drp$Dropped)
summary(state_drp$total)

#  Predicted counts 
stdrp_prd <-tm_shape(state_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("total",
              palette=viridis(100),
              style= "cont",
              title="pop",
              legend.show=F,
              breaks = c(2500000, 3500000, 4500000, 5500000, 7500000,10500000, 14500000, 19000000)
  )+
  tm_legend(outside = TRUE, text.size=0.7)+
  tm_layout(frame=T) 


# Uncertainty
summary(state_drp$uncertainty)
stdrp_uncert <-tm_shape(state_drp) +
  tm_borders(lwd=1, col="black")+
  tm_polygons("uncertainty",
              #palette=viridis(100),
              palette="Oranges",
              style = "cont",
              title="Uncertainty",
              legend.show=F,
              breaks = c(0.01, 0.05, 0.1,0.2, 0.3, 0.4, 0.5, 0.8,1.3)
  )+
  tm_legend(outside = TRUE, text.size=0.7)+
  tm_layout(frame=T) 
  #tm_facets(by="state")+
  #tm_layout(panel.label.size = 2.2) 

tmap_arrange(stdrp_prd, stdrp_uncert,
             lgdrp_prd,lgdrp_uncert,
             wd_drp_tot,wd_drp_uncert,
             nrow = 3, ncol= 2)

save.image(paste0(results_path, "/nga_ward_level_model_07_02_2025.Rdata"))
