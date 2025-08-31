

###---BAYESIAN GEOSTATISTICAL POP MODEL PAPER SIMULATION STUDY CODES
###---NNANATU ET AL.(2025)
### - Simulation Study

rm(list=ls())

library(INLA); library(raster);
library(gtools); library(sp); library(spdep)
library(fields); library(mvtnorm); library(gtools)
library(geoR);library(actuar);library(viridisLite)
require(grid);require(gridExtra);require(lattice);require(tidyverse)
library(dplyr); library(sf); library(tmap); library(tmaptools)

# Set working directory for input and output files
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/NGA/CHRIS_N/paper/simstudy"
out_path <- paste0(path , "/output")
data_path <- paste0(path , "/data")

##-------------------------------------------------------------------------------------------------------------
#dim(dat_covs <- read.csv(paste0(cov_path,"/mod_covs/MICS_LGA_covariates.csv")))# Model fitting covariates
dim(shp <- st_read(paste0(data_path, "/shapefile.gpkg"))) # MICS shapefile

ls()
#path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3/Working/nga/Chris_N/paper1"
#shp_path <- paste0(path, "/application/data/Input_Settlement_Boundaries")
#out_path <- paste0(path, "/sim_study/output_revised1")
# visualise
shp %>%
  ggplot2::ggplot() +
  ggplot2::geom_sf()

# create grid

nga_grid <-
  sf::st_make_grid(shp,
                   # cellsize = 0.1,
                   n = c(100, 100),
                   what = "polygons",
                   square = T)

# Add grid IDs as data.frame
nga_grid <-
  sf::st_sf(id = 1:length(lengths(nga_grid)),
            nga_grid)
dim(nga_grid) # 46320 grid cells (approximately 60 per LGA)

# Visualise
nga_grid %>%
  ggplot2::ggplot() +
  ggplot2::geom_sf()


# Crop to NGA extents
nga_grid_cropped <-
  sf::st_intersection(nga_grid,
                      shp %>% st_make_valid())


dim(nga_grid)# 10000 grid cells
dim(nga_grid_cropped) # 6733 grid cells left after croping to CMR extents

# Visualise
nga_grid_cropped %>%
  ggplot2::ggplot() +
  ggplot2::geom_sf()

# check CRS - coordinates reference system
st_crs(nga_grid_cropped)

# Reproject to lon - lat
nga_shp_tf <- st_transform(nga_grid_cropped,
                           crs = "+proj=longlat +datum=WGS84")
st_crs(nga_shp_tf) # lon-lat

# Convert to spatial object
nga_shp_sp <- as(st_geometry(nga_shp_tf),"Spatial")
coords <- data.frame(coordinates(nga_shp_sp)) # extract the centroids
names(coords)
coords <- coords %>% rename(x = X1, y = X2)
plot(coords)
dim(coords)

# visulise
plot(nga_shp_sp, col="white")
points(coords[,1], coords[,2],
       col=2, cex=0.5, pch="*")



####---Model fit metrics function
model_metrics <- function(obs, pred)
{
  residual = pred - obs
  MAE = mean(abs(residual), na.rm=T)#MAE
  MAPE = (1/length(obs))*sum(abs((obs-pred)/obs))*100#MAPE
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])

  output <- list(MAPE = MAPE,
    MAE  = MAE ,
    #BIAS = abs(BIAS),
    RMSE = RMSE,
    corr = corr)
  return(output)
}

#model_metrics(obs, pred, upper, lower)


#---------------------------------------------------------------------------------------
cover <- seq(0.2,1, by=0.2)

#---Create the data frame
dim(dat.all <- as.data.frame(coords))


####------Build the mesh
# dim(coord)
coord <- cbind(coords$x, coords$y) # must be a matrix object
bnd <- inla.nonconvex.hull(coord, -0.035, -0.04, resolution = c(100, 100))
mesh <- inla.mesh.2d(boundary = bnd, max.edge=c(0.2,1),
                     offset = c(0.2, 0.7),
                     cutoff = 0.2)
par(mfrow=c(1,1))
plot(mesh)
mesh$n # 1133 mesh nodes
#--%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##---Specify SPDE parameters
rho <- 0.37 #--range
nu <- 1 #--smooth parameter

sigma0 <- 0.5#--marginal variance
kappa0 <- sqrt(8*nu)/rho #--scale parameters
(tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)) #--precision

#---SPDE
spde <- inla.spde2.matern(mesh,
                          B.tau = matrix(c(log(tau0), -1, +1), nrow=1, ncol=3),
                          B.kappa = matrix(c(log(kappa0), 0, -1), nrow=1, ncol=3),
                          theta.prior.mean = c(0,0),
                          theta.prior.prec = c(0.1, 0.1))


#--Precision Matrix
Q <- inla.spde2.precision(spde=spde, theta=c(0,0))

#---Simulate the GMRF
sam <- as.vector(inla.qsample(
  n = 1, Q = Q, seed=100))
#length(sam)

###---Build projector matrix A
A <- inla.spde.make.A(mesh=mesh, loc=coord);dim(A)

# Spatial Random effects
sp.raneff <- as.vector(A %*% sam)
length(sp.raneff)
hist(sp.raneff)

##----specify the observation indices for estimation
iset <- inla.spde.make.index(name = "spatial.field", spde$n.spde)


#-- generate the Covariates
nn <- nrow(dat.all)
covs = cbind(1,runif(nn,1,10),rnorm(nn, 2, 2),rpois(nn,4),rnorm(nn, 4, 3), runif(nn))
apply(covs, 2, mean) #---check averages
dim(dat.all)

# Standardization
stdize <- function(x)
{
  #return(stdze <- (x-min(x, na.rm=T))/(max(x, na.rm =T)-min(x, na.rm =T)))
  return(stdze <- (x-mean(x, na.rm=T))/(sd(x, na.rm =T)))
}

covs_std <- cbind(covs[,1],apply(covs[,-1], 2, stdize))

# Add covariates to data
head(dat.all)
rownames(dat.all) <- NULL
dim(ddat <- cbind(dat.all, covs_std))
ddat <- data.frame(ddat)
colnames(ddat) <- c("lon", "lat", "x1", "x2", "x3", "x4", "x5", "x6")
names(ddat)
head(ddat)
###############
#Parameters
#
#dim(zpred <- covs_std)
dim(zpred <- covs)

require(MASS)
#--generate data type random effects


# print(AdmUnit[i]^2)
#####----------ALTERNATIVELY---------########
#---Simulate Pop Count and Building Count separately and calculate pop density
#---see the alternative codes below
###---Simulate building count
sigma_b <- 0.55
epsb <- rnorm(nrow(coord), 0, sigma_b)#--iid random effect for building count
#betaB <- c(3.51, 0.46, 0.85, 0.21, 0.58, 0.235) #---betas- fixed effects
betaB <- c(2.21, 0.06, 0.15, 0.21, 0.18, 0.27)
bld <- lambdaB <- numeric(nn) #
for (i in 1:nn)
{
  lambdaB[i] <- exp(zpred[i,1]*betaB[1] + zpred[i,2]*betaB[2] +
                      zpred[i,3]*betaB[3] + zpred[i,4]*betaB[4] +
                      zpred[i,5]*betaB[5] + zpred[i,6]*betaB[6] + sp.raneff[i] + epsb[i])
  bld[i] <- rpois(1, lambdaB[i])
}
bld
min(bld); max(bld)
hist(bld); hist(log(bld))
mean(bld); var(bld);
summary(bld)


# set any 0 buildings to 1
bld[bld==0]=1

###---Simulate Population count
sigma_p <- 0.015
epsp <- rnorm(nrow(coord), 0, sigma_p)
betaP <- c(3.50, 0.41, 0.08, 0.04, 0.15, 0.22) #--betas - fixed effects

pop <- lambdaP <- numeric(nn)



for (i in 1:nn)
{
  lambdaP[i] <- exp(zpred[i,1]*betaP[1] + zpred[i,2]*betaP[2] +
                      zpred[i,3]*betaP[3] + zpred[i,4]*betaP[4] +
                      zpred[i,5]*betaP[5] + zpred[i,6]*betaP[6]  + sp.raneff[i] +  epsp[i]
  )
  pop[i] <- rpois(1, lambdaP[i])
}
sum(pop)
hist(pop); hist(log(pop))
mean(pop); var(pop)
summary(pop)
# set any 0 population to 1
pop[pop==0]=1

#--------Add to dataset
ddat$bld <- bld
ddat$pop <- pop

# define population density
ddat$dens <- ddat$pop/ddat$bld #----population density
hist(ddat$dens); hist(log(ddat$dens))


datam <- ddat
head(ddat)
names(datam); head(datam)
datam[,paste0("x",1:6)] <- covs
head(datam)



## Visualise the simulated full data
dim(dat_sf <- st_as_sf(datam, coords=c("lon", "lat"), crs=st_crs(shp)))
nga_grid_cropped$pop <- dat_sf$pop

sim_map <- tm_shape(nga_grid_cropped) +
  tm_borders(col="black", lwd=0.5)+
  tm_polygons(col="pop",title="Simulated \n Population  counts",
              palette = viridis(100))+

  #breaks = c(1.23, 1.64, 2.36, 3.23,3.36,3.55))+
  tm_layout(frame=F, legend.outside = T,
            legend.text.size = 1, legend.title.size = 1.2) #+
 # tm_compass(position = c("left","top"),
             #text.size=1.5, size = 1)+
  #tm_scale_bar(position = c("right","bottom"),
              # text.size=1, width = 0.3)



coverp <- c(0.1,  0.3, 0.5, 0.7, 0.9)# percentage missing
coverb <- c(0.1,  0.2, 0.3, 0.4, 0.5) # percentage biased
coverm <- c(0.1, 0.2, 0.3, 0.4, 0.5) # magnitude of bias


for(i in 1:length(coverp))
{
  #i=4
  result_path1 <- paste0(out_path,"/outputs_for_", coverp[i]*100,"%","missingness")
  if (file.exists(result_path1)){
    setwd(file.path(result_path1))
  } else {
    dir.create(file.path(result_path1))
    setwd(file.path(result_path1))
  }

  dat.clust <- datam
  pm <- coverp[i]
#  print(pm*100)#--

  # Apply missing values and inflated values to the cluster data
  dat.clust$popm <-  dat.clust$pop
  ind.obs <- sample(nrow(dat.clust), pm*nrow(dat.clust))
  dat.clust$popm[ind.obs] = NA



  for(j in 1:length(coverb))
  {
   # j=3

     pb = coverb[j]
    result_path2 <- paste0(result_path1,"/", coverb[j]*100, "%","_biased")
    if (file.exists(result_path2)){
      setwd(file.path(result_path2))
    } else {
      dir.create(file.path(result_path2))
      setwd(file.path(result_path2))
    }
    #print(c(pm*100,pb*100))#---

    ind.pm <- which(!is.na(dat.clust$popm))
    indd <- sample(ind.pm, length(ind.pm)*pb)

  for(k in 1:length(coverm))
  {

    # k = 3
    bm = coverm[k]
    result_path3 <- paste0(result_path2,"/", coverm[k]*100, "%","_magnitude_of_bias")
    if (file.exists(result_path3)){
      setwd(file.path(result_path3))
    } else {
      dir.create(file.path(result_path3))
      setwd(file.path(result_path3))
    }

    print(c(paste0(pm*100,"% ", "missing samples:"),
         paste0(pb*100,"% ","biased samples with"),
         paste0(bm*100, "% ", "magnitude of bias")))#---

    dat.clust$popm[indd] = dat.clust$popm[indd] + round(dat.clust$popm[indd]*bm)



    coords.clust <- cbind(dat.clust$lon, dat.clust$lat)

    ###--Cluster level mesh, spde and projection matrix
    bnd_clust <- inla.nonconvex.hull(coords.clust, -0.035, -0.05, resolution = c(100, 100))
    mesh.clst <- inla.mesh.2d(boundary =bnd_clust, max.edge=c(0.1,1),
                              offset = c(0.1, 0.3),
                              cutoff = 0.2)

    par(mfrow=c(1,1))
    plot(mesh.clst)
    with(dat.clust, points(lon, lat, pch=16, col="red"))
    mesh.clst$n

    ###---Build projector matrix A for clusters
    A.clst<-inla.spde.make.A(mesh=mesh.clst,loc=as.matrix(coords.clust));dim(A.clst)

    ##---Create the SPDE for clusters
    spde.clst <- inla.spde2.matern(mesh.clst, alpha=2)

    ##----specify the observation indices for estimation  for clusters
    iset.clst <- inla.spde.make.index(name = "s", spde.clst$n.spde)

    #@@@@@@@@@@------------------------------------------------------------------------------
    library(car) ##--For calculating variance inflation factor (vif)

    dat.clust$popm1 <- dat.clust$popm # leave as it is
    dat.clust$dens1 <- dat.clust$popm1/dat.clust$bld
    dat.clust$dens1[is.infinite(dat.clust$dens1)]=NA
    dat.clust$dens1[is.nan(dat.clust$dens1)]=NA


    dat.clust$popm2 <- dat.clust$popm # mean imputation
    dat.clust$popm2[indd] <- round(mean(dat.clust$popm2[-indd], na.rm=T))
    dat.clust$dens2 <- dat.clust$popm2/dat.clust$bld
    dat.clust$dens2[is.infinite(dat.clust$dens2)]=NA
    dat.clust$dens2[is.nan(dat.clust$dens2)]=NA


    dat.clust$popm3 <- dat.clust$popm # drop or set to NA
    dat.clust$popm3[indd] <- NA
    dat.clust$dens3 <- dat.clust$popm3/dat.clust$bld
    dat.clust$dens3[is.infinite(dat.clust$dens3)]=NA
    dat.clust$dens3[is.nan(dat.clust$dens3)]=NA


    ## Covariates and data stacking for both grid and cluster
    #dat.fit.grd[, c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.grid[, c("x2", "x3", "x4", "x5", "x6")], 2, stdize2)
    dat.clust[, c("x2", "x3", "x4", "x5", "x6")] <- apply(dat.clust[, c("x2",
                                                                            "x3", "x4", "x5", "x6")], 2, stdize)

    #---
    dat.clust$ID <- 1:nrow(dat.clust)
    covars.clust <- dat.clust[,c("ID","x2", "x3", "x4", "x5", "x6", "bld")]; dim(covars.clust) ##---Population density


    #---Build the stacks and fit Bayesian Hierarchical models 
    
    ## Base 
    stk1<- inla.stack(data=list(y=dat.clust$dens1), #the response

                              A=list(A.clst,1),  #the A matrix; the 1 is included to make the list(covariates)

                              effects=list(c(list(Intercept=1), #the Intercept
                                             iset.clst),  #the spatial index
                                           #the covariates
                                           list(covars.clust)
                              ),
                              #this is a quick name so you can call upon easily
                              tag='est')
    
    f1 <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + f(s, model=spde.clst) +
      f(ID, model='iid')
    m1 <-inla(f1, #the formula
              data=inla.stack.data(stk1,spde=spde.clst),  #the data stack
              family= 'gamma',   #which family the data comes from
              control.predictor=list(A=inla.stack.A(stk1),compute=TRUE),  #compute gives you the marginals of the linear predictor
              control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
              verbose = FALSE)
    summary(m1)
    
    ind1 <-inla.stack.index(stk1, "est")$data
    fit1 <- exp(m1$summary.linear.predictor[ind1,"mean"])
    sum(dat.clust$pred1 <- fit1*dat.clust$bld)
    plot(dat.clust$pop, dat.clust$pred1)
    
    


    ## Imputed
    stk2<- inla.stack(data=list(y=dat.clust$dens2), #the response

                      A=list(A.clst,1),  #the A matrix; the 1 is included to make the list(covariates)

                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset.clst),  #the spatial index
                                   #the covariates
                                   list(covars.clust)
                      ),
                      #this is a quick name so you can call upon easily
                      tag='est')
    
    f2 <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + f(s, model=spde.clst) +
      f(ID, model='iid')
    m2 <-inla(f2, #the formula
              data=inla.stack.data(stk2,spde=spde.clst),  #the data stack
              family= 'gamma',   #which family the data comes from
              control.predictor=list(A=inla.stack.A(stk2),compute=TRUE),  #compute gives you the marginals of the linear predictor
              control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
              verbose = FALSE)
    summary(m2)
    
    ind2 <-inla.stack.index(stk2, "est")$data
    fit2 <- exp(m2$summary.linear.predictor[ind2,"mean"])
    sum(dat.clust$pred2 <- fit2*dat.clust$bld)
    plot(dat.clust$pop, dat.clust$pred2)
    
    

    ## Dropped
    stk3<- inla.stack(data=list(y=dat.clust$dens3), #the response

                      A=list(A.clst,1),  #the A matrix; the 1 is included to make the list(covariates)

                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset.clst),  #the spatial index
                                   #the covariates
                                   list(covars.clust)
                      ),
                      #this is a quick name so you can call upon easily
                      tag='est')
    f3 <- y ~ -1 + Intercept + x2 + x3 + x4 + x5 + x6 + f(s, model=spde.clst) +
      f(ID, model='iid')
    m3 <-inla(f3, #the formula
              data=inla.stack.data(stk3,spde=spde.clst),  #the data stack
              family= 'gamma',   #which family the data comes from
              control.predictor=list(A=inla.stack.A(stk3),compute=TRUE),  #compute gives you the marginals of the linear predictor
              control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
              verbose = FALSE)
    summary(m3)
    
    ind3 <-inla.stack.index(stk3, "est")$data
    fit3 <- exp(m3$summary.linear.predictor[ind3,"mean"])
    sum(dat.clust$pred3 <- fit3*dat.clust$bld)
    plot(dat.clust$pop, dat.clust$pred3)


    # Model cross validation and fit metrics
    met1 <- unlist(model_metrics(dat.clust$pop, dat.clust$pred1))
    met2 <- unlist(model_metrics(dat.clust$pop, dat.clust$pred2))
    met3 <- unlist(model_metrics(dat.clust$pop, dat.clust$pred3))

    ## Compare fir metrics
    met_dat <- data.frame(rbind(met1, met2, met3))
    met_dat$method <- c("base", "imputed", "dropped")
    print(met_dat)


    # save metrics data
    write.csv(met_dat, file=paste0(result_path3,"/fit_metrics.csv"))
   }
  }
}


#####---------------------------------------------------------------------------------
# Extract the posterior samples
####----------------------------------------------------------------------------------
path <- "//worldpop.files.soton.ac.uk/Worldpop/Projects/WP517763_GRID3_Phase2/Working/NGA/CHRIS_N/paper/simstudy"
out_path <- paste0(path , "/output")
data_path <- paste0(path , "/data")
coverp <- c(0.1,  0.3, 0.5, 0.7, 0.9)# percentage missing
coverb <- c(0.1,  0.2, 0.3, 0.4, 0.5) # percentage biased
coverm <- c(0.1, 0.2, 0.3, 0.4, 0.5) # magnitude of bias

metric <- c("MAPE","MAE", "RMSE", "corr")
method <- c("base", "imputed", "dropped")

n.pop.miss <- length(coverp)
n.pop.bias <- length(coverb)
n.mag.bias <- length(coverm)
n.metric <- length(metric)
n.method <- length(method)

#dat_met <- matrix(0, nrow=n.size, ncol=5)
dat_met1 <- list()
dat_met2 <- list()
dat_met3 <- list()


for(i in 1:n.pop.bias)
{
  result_path1 <- paste0(out_path,"/outputs_for_", coverp[i]*100,"%","missingness")
  for(j in 1:n.pop.bias)
  {
    result_path2 <- paste0(result_path1,"/", coverb[j]*100, "%","_biased")

  for(k in 1:n.mag.bias)
  {
    result_path3 <- paste0(result_path2,"/", coverm[k]*100, "%","_magnitude_of_bias")
    met0 <- read.csv(paste0(result_path3, "/fit_metrics.csv"))
    met0$mag_bias <- rep(coverm[k]*100, n.method) # magnitude of bias
    met0$prop_bias <- rep(coverb[j]*100, n.method) # percentage of biased population
    met0$pop_miss <- rep(coverp[i]*100, n.method) # proportion of missing observations


    dat_met1[[k]] = met0
  }
  dat_met2[[j]] = dat_met1
  }
  dat_met3[[i]] = dat_met2
}

unnest_list <- unlist(unlist(dat_met3, recursive = FALSE),
                      recursive = FALSE)  #--unnest the list
dim(metrics <- do.call(rbind, unnest_list)) # 225 datasets
metrics <- metrics %>% dplyr::select(-X) # drop the variable 'X'
metrics



setwd(out_path) # set working directory
write.csv(metrics, "combined_fit_metrics.csv", row.names=FALSE) # save the combined fit data

# Convert to long format for plotting
# install.packages("reshape2")
library(reshape2)
dim(met_long <- melt(metrics, id.vars=c("method","prop_bias","mag_bias", "pop_miss"),
                     value.name="estimate", variable.name = "metric"))
met_long$method <- factor(met_long$method)
met_long$method = factor(met_long$method,
                         levels = c("base", "imputed", "dropped"),
                         labels = c("Base", "Imputed", "Dropped"))


# Recode some key variables 
met_long$pop_miss2 = factor(met_long$pop_miss,
                            levels=c("10", "30", "50", "70", "90"),
                            labels=c("Missing prop: 10%",
                                    "Missing prop: 30%",
                                    "Missing prop: 50%",
                                     "Missing prop: 70%",
                                     "Missing prop: 90%"))

met_long$mag_bias2 = factor(met_long$mag_bias,
                            levels=c("10", "20", "30", "40", "50"),
                            labels=c("Size of bias: 10%",
                                     "Size of bias: 20%",
                                     "Size of bias: 30%",
                                     "Size of bias: 40%",
                                     "Size of bias: 50%"))


met_long$prop_bias2 <- factor(met_long$prop_bias)
write.csv(met_long, "combined_fit_metrics_long.csv", row.names=FALSE)

# Subset data for each individual metric
MAE<- met_long[met_long$metric=="MAE",]
RMSE<- met_long[met_long$metric=="RMSE",]
MAPE<- met_long[met_long$metric=="MAPE",]
CC<- met_long[met_long$metric=="corr",]

## MAE
plot_mae <-ggbarplot(MAE, x = "prop_bias2", y = "estimate",
                    error.plot = "estimate",
                    fill = "method",
                    palette = "lancet",
                    point.size=1.2,
                    size=1)+
  theme_bw()+
  theme(strip.text = element_text(size = 15))+
  facet_grid(pop_miss2~ mag_bias2)

rmae <- ggpar(plot_mae,xlab="Proportion of biased samples (%)",
             ylab="Mean absolute error (MAE)",
             legend = "top", legend.title=element_text("Method:"),
             font.legend=c(18),
             font.label = list(size = 18, face = "bold", color ="red"),
             font.x = c(18),
             font.y = c(18),
             font.main=c(18),
             font.xtickslab =c(18),
             font.ytickslab =c(18),
             xtickslab.rt = 45, ytickslab.rt = 45)
rmae




## MAPE
plot_mape <-ggbarplot(MAPE, x = "prop_bias2", y = "estimate",
                      error.plot = "estimate",
                      fill = "method",
                      palette = "jco",
                      point.size=1.2,
                      size=1)+
  theme_bw()+
  theme(strip.text = element_text(size = 14))+
  facet_grid(pop_miss2~ mag_bias2)

rmape <- ggpar(plot_mape,xlab="Proportion of biased samples (%)",
               ylab="Mean absolute percentage error (MAPE)",
               legend = "top", legend.title=element_text("Method:"),
               font.legend=c(18),
               font.label = list(size = 18, face = "bold", color ="red"),
               font.x = c(18),
               font.y = c(18),
               font.main=c(18),
               font.xtickslab =c(18),
               font.ytickslab =c(18),
               xtickslab.rt = 45, ytickslab.rt = 45)
rmape


## CC
plot_cc <-ggdotplot(CC, x = "prop_bias2", y = "estimate",
                      error.plot = "estimate",
                      fill = "method",
                      palette = "lancet",
                      point.size=1.2,
                      size=2)+
  theme_bw()+
  theme(strip.text = element_text(size = 15))+
  facet_grid(pop_miss2 ~ mag_bias2)

rcc <- ggpar(plot_cc,xlab="Proportion of biased samples (%)",
               ylab="Correlation Coefficient (CC)",
               legend = "top", legend.title=element_text("Method:"),
               font.legend=c(18),
               font.label = list(size = 18, face = "bold", color ="red"),
               font.x = c(18),
               font.y = c(18),
               font.main=c(18),
               font.xtickslab =c(18),
               font.ytickslab =c(18),
               xtickslab.rt = 45, ytickslab.rt = 45)
rcc





### Reduction in relative mae
rrmae_dat <- as.data.frame(read.csv(paste0(out_path, "/reduction_in_relative_mae.csv")))
rrmae_dat1 <- rrmae_dat %>% drop_na(rrmae)%>% filter(mag_bias > 10, prop_bias>10, pop_miss < 90) %>% dplyr::select(c(pop_miss,
                                                              prop_bias,
                                                              mag_bias,
                                                              method,
                                                              rrmae))
rrmae_dat1$rrmae <- round(rrmae_dat1$rrmae*100, 2)
rrmae_imp <- rrmae_dat1[rrmae_dat1$method == "Imputed",]
rrmae_drpd <- rrmae_dat1[rrmae_dat1$method == "Dropped",]

min(rrmae_imp$rrmae);max(rrmae_imp$rrmae) # Imputed: -46.93 to 59.54% 
min(rrmae_drpd$rrmae);max(rrmae_drpd$rrmae) # Dropped: 1.49 to 65.74%


# for imputed samples
rrmae_imp$prop_bias<- factor(rrmae_imp$prop_bias)
rrmae_imp$pop_miss <- factor(rrmae_imp$pop_miss)
rrmae_imp$mag_bias <- factor(rrmae_imp$mag_bias)
rrmae_imp$pop_miss2 = factor(rrmae_imp$pop_miss,
                              levels=c("10", "30", "50", "70"),
                              labels=c("Missing prop: 10%",
                                       "Missing prop: 30%",
                                       "Missing prop: 50%",
                                       "Missing prop: 70%"))

rrmae_imp$mag_bias2 = factor(rrmae_imp$mag_bias,
                              levels=c("20", "30", "40", "50"),
                              labels=c("Size of bias: 20%",
                                       "Size of bias: 30%",
                                       "Size of bias: 40%",
                                       "Size of bias: 50%"))

bar_imp2 <- ggbarplot(rrmae_imp, x = "prop_bias", y = "rrmae",
                       error.plot = "rrmae",
                       fill = "pop_miss",
                       palette = "jco",
                       point.size=1.2,
                       size=1)+
  theme_bw()+
  theme(strip.text = element_text(size = 14))+
  facet_grid(pop_miss2~mag_bias2)

rbar_imp2 <-  ggpar(bar_imp2 , xlab="Proportion of biased population(%)",
                    ylab="Reduction in relative MAE(%)",
                    legend = "none", legend.title=element_text("Proportion of \n missing samples(%):"),
                    font.legend=c(18),
                    font.label = list(size = 18, face = "bold", color ="red"),
                    font.x = c(18),
                    font.y = c(18),
                    font.main=c(18),
                    font.xtickslab =c(18),
                    font.ytickslab =c(18),
                    xtickslab.rt = 45, ytickslab.rt = 45)
rbar_imp2



# for dropped samples
rrmae_drpd$prop_bias<- factor(rrmae_drpd$prop_bias)
rrmae_drpd$pop_miss <- factor(rrmae_drpd$pop_miss)
rrmae_drpd$mag_bias <- factor(rrmae_drpd$mag_bias)
rrmae_drpd$pop_miss2 = factor(rrmae_drpd$pop_miss,
                            levels=c("10", "30", "50", "70"),
                            labels=c("Missing prop: 10%",
                                     "Missing prop: 30%",
                                     "Missing prop: 50%",
                                     "Missing prop: 70%"))

rrmae_drpd$mag_bias2 = factor(rrmae_drpd$mag_bias,
                            levels=c("20", "30", "40", "50"),
                            labels=c("Size of bias: 20%",
                                     "Size of bias: 30%",
                                     "Size of bias: 40%",
                                     "Size of bias: 50%"))

bar_drpd2 <- ggbarplot(rrmae_drpd, x = "prop_bias", y = "rrmae",
          error.plot = "rrmae",
          fill = "pop_miss",
          palette = "jco",
          point.size=1.2,
          size=1)+
  theme_bw()+
  theme(strip.text = element_text(size = 14))+
  facet_grid(pop_miss2~mag_bias2)

rbar_drpd2 <-  ggpar(bar_drpd2 , xlab="Proportion of biased population(%)",
                     ylab="Reduction in relative MAE(%)",
                     legend = "none", legend.title=element_text("Proportion of  \n missing samples(%):"),
                     font.legend=c(18),
                     font.label = list(size = 18, face = "bold", color ="red"),
                     font.x = c(18),
                     font.y = c(18),
                     font.main=c(18),
                     font.xtickslab =c(18),
                     font.ytickslab =c(18),
                     xtickslab.rt = 45, ytickslab.rt = 45)
rbar_drpd2

ggarrange(rbar_imp2,rbar_drpd2, nrow = 1, ncol=2)
