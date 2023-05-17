# cLHC manuscript
## Sample number optimisation
## Original data set.
## Original data set is data that has already been collected.
## The aim of this script is to compare the collected data to the actual population in terms of covariate coverage

rm(list=ls())

source("R/utils.R")

# Get the current working directory
current_wd <- getwd()
# Set the working directory to the current file location
setwd(current_wd)

## Libraries
library(raster);library(rgdal);library(tripack);library(SDMTools); library(manipulate);library(clhs);library(entropy)

## Point data
# p.dat<- read.table("T1_dat.txt", header = T, sep=",") #change directory as appropriate
observations_df<- read.table("/home/sudipta/repos/clhc_sampling/additional/intersected_covs.txt", header = T, sep=",") #change directory as appropriate
str(observations_df)
# coordinates(p.dat)<- ~ X_REF + Y_REF

## Raster data
# cov.dat <- read.table("T1_covs.txt", header = T, sep=",") # change directory as appropriate
covariates_df <- read.table("/home/sudipta/repos/clhc_sampling/additional/covs_full.txt", header = T, sep = ",")
str(covariates_df)

start_pos <- 3
end_pos <- 5

# p.dat_I <- observations_df[, start_pos:end_pos]

## comparison of population and sample distributions

## Start TEST 1
#Kullback-Leibler (KL) divergence

#Quantiles of the population
# Number of bins
num_bins <- 100
start_pos = 3
end_pos = 5

# quantile matrix (of the covariate data)
quantiles <- generate_covariate_quantile_matrix(
  covariate_data = covariates_df[, start_pos: end_pos],
  number_of_bins = num_bins
)

# Hypercube of population
cov_hypercube <- generate_hypercube(
  covariate_data = covariates_df[, start_pos:end_pos],
  quantiles = quantiles
)

####Compare whole study area covariate space with the selected sample
#sample data hypercube (essentially the same script as for the grid data but just doing it on the sample data)
observations_hypercube <- generate_hypercube(
  covariate_data = observations_df[, start_pos:end_pos],
  quantiles = quantiles
)

#Kullback-Leibler (KL) divergence
kl_divergence_original_samples <- compute_kl_divergence(cov_hypercube, observations_hypercube)
print(kl_divergence_original_samples) # KL divergence if the existing soil sample (N=??)

#### END First Test:


#### Second Test:
##points in polygons routine

#principal component of sample
pca.s = prcomp(p.dat_I[,1:3],scale=TRUE, center=TRUE)
scores_pca1 = as.data.frame(pca.s$x)
# plot the first 2 principal components and convex hull
rand.tr<-tri.mesh(scores_pca1[,1],scores_pca1[,2])
rand.ch<-convex.hull(rand.tr, plot.it=F) #convex hull
pr_poly <- cbind(x=c(rand.ch$x),y=c(rand.ch$y)) # save the convex hull vertices
plot(scores_pca1[,1], scores_pca1[,2], xlab="PCA 1", ylab="PCA 2", xlim=c(min(scores_pca1[,1:2]), max(scores_pca1[,1:2])),ylim=c(min(scores_pca1[,1:2]), max(scores_pca1[,1:2])))
lines(c(rand.ch$x,rand.ch$x[1]), c(rand.ch$y,rand.ch$y[1]),col="red",lwd=1) # draw the convex hull(domain of prediction)


# PCA projection of population onto the principal components
PCA_projection<- predict(pca.s, covariates_df[, start_pos:end_pos]) # Project population onto sample PC
newScores = cbind(x=PCA_projection[,1],y=PCA_projection[,2]) # PC scores of projected population

#plot the polygon and all points to be checked
plot(newScores,xlab="PCA 1", ylab="PCA 2", xlim=c(min(newScores[,1:2]), max(newScores[,1:2])),ylim=c(min(newScores[,1:2]), max(newScores[,1:2])),col='black', main='convex hull of ss')
polygon(pr_poly,col='#99999990')

#create check which points fall within the polygon
specMatch = pnt.in.poly(newScores,pr_poly)
sum(specMatch$pip)/nrow(specMatch)*100 # propertion of new spectra that fall within the convex hull
points(specMatch[which(specMatch$pip==0),1:2],pch='X', col='red')
#### END Second Test:


#### Results needed from this script
# KL statistic
print(klo) # KL divergence if the existing soil sample (N=238)

# points in polygon statistic
sum(specMatch$pip)/nrow(specMatch)*100 # propertion of new spectra the fall within the convex hull

# END SCRIPT
