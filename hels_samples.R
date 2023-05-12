# Algorithm for determining how to allocate additional samples within a study area given some existing samples
# The method is based mainly on the Carre et al (2007) HELS and  algorithm
# Hypercube Evaluation of a Legacy Sample (HELS)
#
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



# when we have some samples from outside, and possibly some in the local area, how can we use that prior knowledge, and sample from areas that are not sampled well. 

# Hels algorithm - Spend the sampling budget incrementally starting with the worst areas - monitor improvement (ongoing, almost done)
#   a) find the quantiles of the covariates in the local area
#   b) find the hypercube of the covariates in the local area
#   c) find the hypercube of the covariates in the observations
#   d) find the mismatch of densities between local area hypercube vs observartions hypercube
#   e) we take new samples based on the mismatch of densities and add to our observation set
#   f) the order in which we sample is determined by the mismatch of oservered hypercube desnity to taht of covariate hypercube density

# xxxxxxxxx
# xxxxxxxxx
# xxxxxxxxx
# xxxxxxxxx
# xxxxxxxxx
#
# ooooooooo
# ooooooooo
# ooooooooo
# ooooooooo
# ooooooooo

# number_of_bins = 25
# local area each pixel - this is local
# ----2----------2------3------- 25 of them

# observed data - this is global
# ----4----------0------4------   25 of them

# Ratio observed/local
# ----2----------0------0.75---------  25 of them

rm(list = ls())
source("R/utils.R")

# Libraries
library(raster); library(sp); library(rgdal); library(clhs)

#working directory
# setwd("Z:/Dropbox/2018/cLHC_samplingPAPER/gitHub/muddles/additional")
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#load r session with all the necesary objects
# load("clhs_samp.RData")


# Number of addtional samples to take
number_additional_samples <- 100

#load raster data
# load(file="HV_coobs.rda")
covariates_df <- read.table("/home/sudipta/repos/clhc_sampling/additional/covs_full.txt", header = T, sep = ",") # change directory as appropriate

# insert cellNos column used later for tracking rows
covariates_df <- cbind(covariates_df[1:2], cellNos = seq(1:nrow(covariates_df)), covariates_df[, 3:ncol(covariates_df)])


# Define some constants for convenience

num_quantiles <- number_additional_samples
start_pos = 4
end_pos = 6

# quantile matrix (of the covariate data)
covariate_quantiles <- generate_covariate_quantile_matrix(
  covariate_data = covariates_df[, start_pos:end_pos],
  number_of_bins = num_quantiles
)

# covariate data hypercube
# This takes a while to do so only do it once if you can
covariate_hypercube <- generate_hypercube(
  covariate_data = covariates_df[, start_pos:end_pos],
  quantiles = covariate_quantiles,
  number_of_bins = num_quantiles
)

cov_dens <- covariate_hypercube / nrow(covariates_df)  # density matrix


# Point data
# for us this contains the global observation sites + their itnerseteced covaraiate values
observations_df <- read.table("/home/sudipta/repos/clhc_sampling/additional/intersected_covs.txt", header = T, sep = ",")

# insert cellNos column used later
observations_df <- cbind(observations_df[1:2], cellNos = 0, observations_df[, 3:ncol(observations_df)])
observations_df <- observations_df[complete.cases(observations_df),]  # filter out any null rows


# sample data hypercube
observations_hypercube <- generate_hypercube(
  covariate_data = observations_df[, start_pos:end_pos],
  quantiles = covariate_quantiles,
  number_of_bins = num_quantiles
)

# Calculate the data density
dat_dens <- observations_hypercube / nrow(observations_df)


#### selecting new samples

rat <- dat_dens / cov_dens # ratio of data density and covariate density
or <- order(rat) # rank the index where the biggest discrepancy is

## indexes of quantiles that are not adequately sampled ie where rat is less than 1
l1 <- which(rat < 1, arr.ind = T)
l1 <- cbind(l1, which(rat < 1))

# What is the level of the greatest discrepancy? (This is important for comparing to later on when we have the additional sample)
indy <- which(l1[, 3] == or[1])
rp <- l1[indy, 1]
rc <- l1[indy, 2]
rat[rp, rc]
observations_hypercube[rp, rc]
covariate_hypercube[rp, rc] # number of pixel (biggest discrepancy)
dat_dens[rp, rc] # data density
cov_dens[rp, rc] # covariate density


# start from the highest discrepancy then work our way down
up_samp <- number_additional_samples
rpos <- 1

# this could be a high fraction indicating that the existing samples are not very representative of the local area
print(compute_kl_divergence(covariate_data =  covariate_hypercube, sample_data = observations_hypercube))
observartions_and_samples_combined <- observations_df
remaining_covariates_df <- covariates_df


while (up_samp != 0) {  # while the number of samples to allocate is greater than 0
  indy <- which(l1[, 3] == or[rpos])
  rp <- l1[indy, 1]
  rc <- l1[indy, 2]

  ex <- floor(nrow(observartions_and_samples_combined) * (dat_dens[rp, rc])) # existing count of samples within the selected quantile
  eq <- ceiling(nrow(observartions_and_samples_combined) * (cov_dens[rp, rc])) # number of samples needed to get to equal density between data and covariates
  sn <- eq - ex # number of samples needed
  if (up_samp < sn) { sn <- up_samp } # just so we dont over allocate


  # covariate selection
  covL <- covariate_quantiles[rp, rc]
  covU <- covariate_quantiles[rp + 1, rc]
  subDat <- remaining_covariates_df[
    remaining_covariates_df[, (start_pos_sub_one + rc)] >= covL &
    remaining_covariates_df[, (start_pos_sub_one + rc)] <= covU,
  ] # subset the covariates that meet the standard

  # training <- sample(nrow(subDat), sn) #random number   # hels aglorithm
  # training <- clhs(remaining_covariates_df[, start_pos:end_pos], size = sn, progress = T, iter = 1000)  # clhs
  training <- clhs(subDat[, start_pos:end_pos], size = sn, progress = T, iter = 1000)  # a mix
  subDat2 <- subDat[training,]

  #remove selected samples from tempD so that repeated sampling does not occur (Is this necessary??)
  remaining_covariates_df <- remaining_covariates_df[!(remaining_covariates_df$cellNos %in% subDat2$cellNos),]

  # Append new data to sampling dataframe
  # previous_data_length = nrow(observations_df)
  observartions_and_samples_combined <- rbind(observartions_and_samples_combined, subDat2)
  h.mat.updated <- generate_hypercube(
    covariate_data = observartions_and_samples_combined[, start_pos:end_pos],
    quantiles = covariate_quantiles,
    number_of_bins = num_quantiles
  )
  print("===============================")
  print(compute_kl_divergence(covariate_data = covariate_hypercube, sample_data = h.mat.updated))

  # adjust the while params
  rpos <- rpos + 1
  up_samp <- up_samp - sn
  print(sn)
}


# Check the sampling density with the additional samples added

#density
h.mat.updated[which(h.mat.updated == 0)] <- 0.000001
h.mat.updated
dat_dens <- h.mat.updated / nrow(observations_df)


#### check

rat <- dat_dens / cov_dens # ratio of data density and covariate density
or <- order(rat) # rank the index where the biggest discrepancy is
or

## indexes of quantiles that are not adequately sampled
l1 <- which(rat < 1, arr.ind = T)
l1 <- cbind(l1, which(rat < 1))
l1
length(rat)
nrow(l1)


# What the the level of the greatest discrepancy?
indy <- which(l1[, 3] == or[1])
rp <- l1[indy, 1]
rc <- l1[indy, 2]
rat[rp, rc]
observations_hypercube[rp, rc]
covariate_hypercube[rp, rc]
dat_dens[rp, rc]
cov_dens[rp, rc]


## The following code does not have too much to do with the algorithm

# Specify the different surveys (original and additional)
observartions_and_samples_combined$survey <- NA
str(observartions_and_samples_combined)
observartions_and_samples_combined[observartions_and_samples_combined$cellNos == 0, "survey"] <- "original"
observartions_and_samples_combined[observartions_and_samples_combined$cellNos != 0, "survey"] <- "hels_samples"

## Spatial points
coordinates(observartions_and_samples_combined) <- ~X_REF + Y_REF
str(observartions_and_samples_combined)

## Coordinate reference systems
# proj4string(dat) <- CRS("+init=epsg:32756")
proj4string(observartions_and_samples_combined) <- CRS("+init=epsg:3577")  # Australian albers


## Write point data to shapefile
writeOGR(observartions_and_samples_combined, ".", "wa_test", "ESRI Shapefile", overwrite_layer = T)


# The following raster has been derived by comparing the similarity between the mulivariate values of the grids and observed data
# Low numbers mean that there are not many data points showing similarity to the grid cell location. High mumber means quite a few
# observations are similar

#raster extraction
#covariates
# files<- list.files(path= getwd(), pattern = ".tif$", full.names = TRUE)
# files
# r1<- raster(files[1])
# r1
# plot(r1)

# DSM_data2<- extract(r1, dat, sp= 1, method = "simple")

# write table to file
# write.table(as.data.frame(DSM_data2), "HELS_dat.txt", sep = ",", col.names = T)


# Doing some summary statistics between the raster grid values and the sample sites for both original and additional data
# DSM_data2$sampleNos
# dat1<- DSM_data2[DSM_data2$survey == "original", ]
# sum(dat1$sampleNos >= 0 & dat1$sampleNos <= 5) / nrow(dat1)
# sum(dat1$sampleNos > 5 & dat1$sampleNos <= 10) / nrow(dat1)
# sum(dat1$sampleNos > 10 & dat1$sampleNos <= 20) / nrow(dat1)
# sum(dat1$sampleNos > 20 & dat1$sampleNos <= 40) / nrow(dat1)
# sum(dat1$sampleNos > 40) / nrow(dat1)
#
# dat2<- DSM_data2[DSM_data2$survey !=1, ]
# sum(dat2$sampleNos >= 0 & dat2$sampleNos <= 5) / nrow(dat2)
# sum(dat2$sampleNos > 5 & dat2$sampleNos <= 10) / nrow(dat2)
# sum(dat2$sampleNos > 10 & dat2$sampleNos <= 20) / nrow(dat2)
# sum(dat2$sampleNos > 20 & dat2$sampleNos <= 40) / nrow(dat2)
# sum(dat2$sampleNos > 40) / nrow(dat2)

#save.image("clhs_samp.RData") #save R session
