# Algorithm for determining how to allocate additional samples within a study area given some existing samples
# basically it entails:
# 1. From the covariate data generate a matrix of quantiles (n x p) n = number of samples that are needed. p = number of available covariates
# 2. Using the quantile matrix create Hypercube matrices for both the existing point data and covariates.
# 3. Work out the densities of the elements in the hypercubes (count of obs/pixels within each quantile / total number of data or pixels)
# 4. Work out the ratio of sampling and grid densities... This will identify under and over sampling in the hypercube
# 5. To add additional samples:
#         1. rank the ratios from smallest to largest
#         2. workout the number of samples required to equalise quantile density of grids and sample data
#         3. Repeat step 5.2 until total number of additonal samples have been allocated.
## Created: 18/05/17



rm(list = ls())
source("R/utils.R")

# Libraries
library(raster); library(sp); library(rgdal); library(clhs)

#working directory
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#load r session with all the necesary objects
# load("clhs_samp.RData")


#load raster data
# load(file="HV_coobs.rda")
covariates_df <- read.table("/home/sudipta/repos/clhc_sampling/additional/covs_full.txt", header = T, sep = ",") # change directory as appropriate
covariates_df <- covariates_df[complete.cases(covariates_df),]  # filter out any null rows

# insert cellNos column used later for tracking rows
covariates_df <- cbind(
  covariates_df[1:2],
  cellNos = seq(1:nrow(covariates_df)),
  type='covariate',
  covariates_df[, 3:ncol(covariates_df)]
)


# Define some constants for convenience
start_pos <- 5
end_pos <- 7

# Point data
# for us this contains the global observation sites + their itnerseteced covaraiate values
observations_df <- read.table("/home/sudipta/repos/clhc_sampling/additional/intersected_covs.txt", header = T, sep = ",")

# insert cellNos column used later
observations_df <- cbind(observations_df[1:2], cellNos = 0, type='existing', observations_df[, 3:ncol(observations_df)])
observations_df <- observations_df[complete.cases(observations_df),]  # filter out any null rows
observations_df <- unique(observations_df)

num_quantiles <- 25


# quantile matrix (of the covariate data)
covariate_quantiles <- generate_covariate_quantile_matrix(
  covariate_data = covariates_df[, start_pos:end_pos],
  number_of_bins = num_quantiles
)

# covariate data hypercube
# This takes a while to do so only do it once if you can
covariate_hypercube <- generate_hypercube_vec(
  covariate_data = covariates_df[, start_pos:end_pos],
  quantiles = covariate_quantiles
)

# sample data hypercube
observations_hypercube <- generate_hypercube_vec(
  covariate_data = observations_df[, start_pos:end_pos],
  quantiles = covariate_quantiles
)


# this could be a high fraction indicating that the existing samples are not very representative of the local area
print(compute_kl_divergence(covariate_data =  covariate_hypercube, sample_data = observations_hypercube))

old_and_samples_combined <- observations_df

# TODO: in universe_df it will be best to remove the pixels corresponding to the observations_df from covariates_df
universe_df <- rbind(observations_df, covariates_df)
num_old_observations <- nrow(observations_df)
num_current_observations <- nrow(observations_df)
must_include <- 1:num_old_observations


# Number of addtional samples to take
number_additional_samples <- 110
up_samp <- number_additional_samples
pass_no <- 0

while (up_samp != 0) {  # while the number of samples to allocate is greater than 0
  pass_no <- pass_no + 1
  print(paste("============PASS NO ==========>>>>>>>>>>>>: ", pass_no))
  sn <- 10 # 10 at a time
  if (up_samp < sn) { sn <- up_samp } # just so we dont over allocate
  total_samples_for_this_loop <- sn + num_current_observations
  # covariate selection
  while (TRUE) {   # loop untill all training indices are unique
    training <- clhs(
      universe_df[, start_pos:end_pos],
      size = total_samples_for_this_loop,
      progress = T, iter = 1000,
      use.cpp = TRUE,
      must.include = must_include
    )
    if (!any(duplicated(training))) { break }
    print("===========>>>>clhs returned duplicate samples. Retrying ......")
  }
  # original + all new observations
  old_and_samples_combined <- rbind(old_and_samples_combined, universe_df[training[!(training %in% must_include)],])
  old_and_samples_combined[old_and_samples_combined$type == 'covariate', "type"] <- paste("clhc_sample_pass", pass_no, sep = "_")
  # update must_include and num_training_observations
  must_include <- training
  num_current_observations <- total_samples_for_this_loop

  #remove selected samples from remaining_covariates_df so that repeated sampling does not occur (Is this necessary??)
  # remaining_covariates_df <- remaining_covariates_df[!(remaining_covariates_df$cellNos %in% subDat2$cellNos),]

  # Append new data to sampling dataframe
  observations_hypercube_updated <- generate_hypercube_vec(
    covariate_data = old_and_samples_combined[, start_pos:end_pos],
    quantiles = covariate_quantiles
  )

  # adjust the while params
  up_samp <- up_samp - sn
  print(table(old_and_samples_combined['type']))
  print(
    paste(
    "KL Divergence with ",  total_samples_for_this_loop - num_old_observations, " additional samples",
    compute_kl_divergence(covariate_data = covariate_hypercube, sample_data = observations_hypercube_updated)
    )
  )
}

# Specify the different surveys (original and additional)

# Spatial points
coordinates(old_and_samples_combined) <- ~X_REF + Y_REF

# Coordinate reference systems
proj4string(old_and_samples_combined) <- CRS("+init=epsg:3577")  # Australian albers


## Write point data to shapefile
writeOGR(old_and_samples_combined, ".", "wa_test", "ESRI Shapefile", overwrite_layer = T)


# composite quantiles raster output
composite <- composite_from_quantiles(covariates_df[, start_pos: end_pos], quantiles = covariate_quantiles)
raster <- rasterFromXYZ(cbind(covariates_df[, c("X_REF", "Y_REF")], composite))
proj4string(raster) <- CRS("+init=epsg:3577")  # Australian albers
writeRaster(raster, filename = "./compositve_quantiles.tif", format = "GTiff", overwrite = TRUE)

composite_sum <- rowSums(composite)
r1 <- rasterFromXYZ(cbind(covariates_df[, c("X_REF", "Y_REF")], composite_sum=composite_sum))
proj4string(r1) <- CRS("+init=epsg:3577")  # Australian albers
plot(r1)
writeRaster(r1, filename = "./compositve_quantiles_sum.tif", format = "GTiff", overwrite = TRUE)

