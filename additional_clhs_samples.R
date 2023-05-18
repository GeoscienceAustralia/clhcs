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
library(raster); library(sp); library(rgdal); library(clhs); library(dplyr)

#working directory
# setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#load r session with all the necesary objects
# load("clhs_samp.RData")

coord_col_names <-  c("X_REF", "Y_REF")

#load raster data
# load(file="HV_coobs.rda")
covariates_df <- read.table("/home/sudipta/repos/clhc_sampling/additional/covs_full.txt", header = T, sep = ",") # change directory as appropriate
covariates_df <- covariates_df[complete.cases(covariates_df),]  # filter out any null rows

# insert cellNos column used later for tracking rows
covariates_df <- cbind(
  covariates_df[, coord_col_names],
  type='covariate',
  covariates_df[, 3:ncol(covariates_df)]
)


# Define some constants for convenience
start_pos <- 4L
end_pos <- 6L

# Point data
# for us this contains the global observation sites + their itnerseteced covaraiate values
observations_df <- read.table("/home/sudipta/repos/clhc_sampling/additional/intersected_covs.txt", header = T, sep = ",")

# insert cellNos column used later
observations_df <- cbind(observations_df[1:2], type='existing', observations_df[, 3:ncol(observations_df)])
observations_df <- observations_df[complete.cases(observations_df),]  # filter out any null rows
observations_df <- unique(observations_df)

num_quantiles <- 25L


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
# note observations can be outside the quantiles of the covariates - do we cannot use the vectorised version
# of this funcion
observations_hypercube <- generate_hypercube(
  covariate_data = observations_df[, start_pos:end_pos],
  quantiles = covariate_quantiles
)


# this could be a high fraction indicating that the existing samples are not very representative of the local area
print(compute_kl_divergence(covariate_data =  covariate_hypercube, sample_data = observations_hypercube))

# initially just contains the existing samples
existing_and_clhs_samples_df <- observations_df

# TODO: in universe_df it will be best to remove the pixels corresponding to the observations_df from covariates_df
universe_df <- rbind(observations_df, covariates_df)
num_old_observations <- nrow(observations_df)
num_current_observations <- nrow(observations_df)
must_include <- 1:num_old_observations


# Number of addtional samples to take
number_additional_samples <- 110L
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
  existing_and_clhs_samples_df <- rbind(existing_and_clhs_samples_df, universe_df[training[!(training %in% must_include)],])
  existing_and_clhs_samples_df[existing_and_clhs_samples_df$type == 'covariate', "type"] <- paste("clhc_sample_pass", pass_no, sep = "_")
  # update must_include and num_training_observations
  must_include <- training
  num_current_observations <- total_samples_for_this_loop

  #remove selected samples from remaining_covariates_df so that repeated sampling does not occur (Is this necessary??)
  # remaining_covariates_df <- remaining_covariates_df[!(remaining_covariates_df$cellNos %in% subDat2$cellNos),]

  # Append new data to sampling dataframe
  observations_hypercube_updated <- generate_hypercube(
    covariate_data = existing_and_clhs_samples_df[, start_pos:end_pos],
    quantiles = covariate_quantiles
  )

  # adjust the while params
  up_samp <- up_samp - sn
  print(table(existing_and_clhs_samples_df['type']))
  print(
    paste(
    "KL Divergence with ",  total_samples_for_this_loop - num_old_observations, " additional samples",
    compute_kl_divergence(covariate_data = covariate_hypercube, sample_data = observations_hypercube_updated)
    )
  )
}


existing_and_clhs_samples_spatial_df <- existing_and_clhs_samples_df


# Spatial points
coordinates(existing_and_clhs_samples_spatial_df) <- ~X_REF + Y_REF
# Coordinate reference systems
proj4string(existing_and_clhs_samples_spatial_df) <- CRS("+init=epsg:3577")  # Australian albers
## Write point data to shapefile
writeOGR(existing_and_clhs_samples_spatial_df, ".", "wa_test", "ESRI Shapefile", overwrite_layer = T)


# composite quantiles raster output
composite <- composite_from_quantiles(covariates_df[, start_pos: end_pos], quantiles = covariate_quantiles)
raster <- rasterFromXYZ(cbind(covariates_df[, c("X_REF", "Y_REF")], composite))
proj4string(raster) <- CRS("+init=epsg:3577")  # Australian albers
writeRaster(raster, filename = "./composite_quantiles.tif", format = "GTiff", overwrite = TRUE)

composite_sum <- rowSums(composite)
r1 <- rasterFromXYZ(cbind(covariates_df[, c("X_REF", "Y_REF")], composite_sum=composite_sum))
proj4string(r1) <- CRS("+init=epsg:3577")  # Australian albers
plot(r1)
writeRaster(r1, filename = "./composite_quantiles_sum.tif", format = "GTiff", overwrite = TRUE)

# create a shapefile with the same covariate values as those of LHC samples
clhc_samples_df <- existing_and_clhs_samples_df %>% dplyr::filter(type != 'existing')
composite_with_coords <- as.matrix(cbind(covariates_df[, c("X_REF", "Y_REF")], composite))

# 2 coordinates + (end_pos-start_pos + 1) covaraites + 1 for the class
classified_sample_quantiles <- matrix(0, nrow = 0, ncol=(2 + (end_pos-start_pos + 1) + 1))

for (i in 1:nrow(clhc_samples_df)) {
  this_sample_quantiles <- composite_with_coords
  for (j in seq(start_pos:end_pos)) {
    dd <- clhc_samples_df[i, start_pos+j-1]
    q <- findInterval(dd, covariate_quantiles[, j])
    this_sample_quantiles <- this_sample_quantiles[this_sample_quantiles[, j+2]==q, ]  # + 2 for the coords
    this_sample_quantiles <- matrix(this_sample_quantiles, ncol = 2 + (end_pos-start_pos + 1))
  }
  classified_sample_quantiles <- rbind(
    classified_sample_quantiles,
    cbind(this_sample_quantiles, i)
  )
}

classified_sample_quantiles_df <- as.data.frame(
  classified_sample_quantiles
)
colnames(classified_sample_quantiles_df) <- c(coord_col_names, names(covariates_df[, start_pos:end_pos]), "clhs_class")

coordinates(classified_sample_quantiles_df) <- ~X_REF + Y_REF
# Coordinate reference systems
proj4string(classified_sample_quantiles_df) <- CRS("+init=epsg:3577")  # Australian albers
## Write point data to shapefile
writeOGR(classified_sample_quantiles_df, ".", "clhs_samples", "ESRI Shapefile", overwrite_layer = T)