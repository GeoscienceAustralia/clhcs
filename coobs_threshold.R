## Algorithm for helping determine the locations of additional samples given the existence of prior sampling.
### This method is based on a count of observations approach.
### We want to determine at every grid cell or pixel, how many observations are similar in terms of the covariate space.

###Notes
# The example below only considers numerical data
# The example below has been created using some supplied data. Will need to be adapted for other datasets

rm(list = ls())

#Libraries
library(raster); library(sp); library(rgdal); library(snow); library(doParallel); library(Rfast); library(Rmpi)

# Given pre-exisiting samples how do we incorporate those and then do a cLHC sample?
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Point obsevation data + intersected covaraites at those locations
# for us this contains the global observation sites + their itnerseteced covaraiate values
observations_df <- read.table("/home/sudipta/repos/clhc_sampling/additional/intersected_covs.txt", header = T, sep = ",")

# covariates in the small/local area we are interested in clhs
covariates_df <- read.table("/home/sudipta/repos/clhc_sampling/additional/covs_full.txt", header = T, sep = ",") # change directory as appropriate

start_pos <- 3
end_pos <- 5

#covariance matrix of the covariates - we have only 4 covariates so 3:6
covariate_cov <- as.matrix(cov(covariates_df[, start_pos:end_pos]))

# filter out rows with any null values
observations_df <- observations_df[complete.cases(observations_df),]

calculate_func <- function(i, covariate_df, observations_df, X_dat, cov_mat) {
  # Pixel distance
  pix <- covariate_df[i,] #pixel values
  # pixDist<- mahalanobis(x = X_covariates , center = as.matrix(pix), cov = covMat) # calculate distance of selected pixel to all other pixels
  pix_dist <- Rfast::mahala(x = observations_df, mu = as.matrix(pix), sigma = cov_mat) # calculate distance of selected pixel to all other pixels
  min_pix <- 0.0   # min(pixDist) # minimum distance (will be 0 always)
  max_pix <- quantile(pix_dist, probs = 0.975) # maximum distance  ##the probs variable could change   ####### Hack to avoid outliers

  # data distance
  # datDist<- mahalanobis(x = X_dat, center = as.matrix(pix), cov = covMat) #calculate distance of observations to all other pixels
  dat_dist <- Rfast::mahala(x = X_dat, mu = as.matrix(pix), sigma = cov_mat) #calculate distance of observations to all other pixels
  dat_ndist <- (dat_dist - min_pix) / (max_pix - min_pix) # standardarise

  dat_ndist[dat_ndist > 1] <- 1 #if the datNdist is greater than 1 that means it datDist is greater than maxDist ##HACK
  dat_ndist <- 1 - dat_ndist  # Higher values mean more similar

  # count how many obs are above a given threshold
  return(sum(dat_ndist >= 0.975))
}


# Begin a parallel cluster and register it with foreach:
cpus <- 4 # The number of nodes/cores to use in the cluster
cl <- makeCluster(spec = cpus, type = 'SOCK')
# Register the cluster with foreach:
registerDoParallel(cl)
X_covariates <- as.matrix(covariates_df[, start_pos:end_pos])
X_dat <- as.matrix(observations_df[, start_pos:end_pos])
coobs_classes <- foreach(i = seq_len(nrow(covariates_df)), .combine = "c") %dopar% {  #para mode
  # for (i in 1:1000){  #sequential mode
  calculate_func(i, X_covariates, X_covariates, X_dat, covariate_cov)
}

# stop cluster
stopCluster(cl)

# prepare grid outputs and export
# tempD$sampleNOS<- oper1
# r1 <- rasterFromXYZ(tempD[,c("X_REF", "Y_REF", "sampleNOS")])
# plot(r1)
# writeRaster(r1, filename="coobs.tif", format="GTiff", overwrite=TRUE)

# library("Rmpi")
# ns <- mpi.universe.size() - 1
# mpi.spawn.Rslaves(nslaves=3)
# Tell all slaves to identify themselves
# mpi.remote.exec(paste("I am",mpi.comm.rank(),"of",mpi.comm.size()))
# oper1.mpi = unlist(mpi.applyLB(1:nrow(tempD), FUN=calculate_func, tempD=tempD[, start_pos:end_pos], X_covariates=X_covariates, X_dat=X_dat, covMat=covMat))
# oper1.mpi[101:nrow(tempD)] <- 0 # hack for testing

# prepare grid outputs and export
# tempD$sampleNOS<- oper1.mpi

covariates_df$sampleNOS <- coobs_classes
r1 <- rasterFromXYZ(covariates_df[, c("X_REF", "Y_REF", "sampleNOS")])
plot(r1)
writeRaster(r1, filename = "./global_local_coobs.tif", format = "GTiff", overwrite = TRUE)

# Tell all slaves to close down, and exit the program
# mpi.close.Rslaves(dellog = FALSE)
# mpi.quit()
