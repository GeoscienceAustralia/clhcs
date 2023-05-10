# compute KL divergence
############################################################
library(entropy)

compute_kl_divergence <- function(cov.mat, h.mat) {
  # Kullback-Leibler (KL) divergence
  klo <- vector(length = ncol(cov.mat), mode = "numeric")
  for (col in 1:ncol(cov.mat)) {
    klo[col] <- KL.empirical(cov.mat[, col], h.mat[, col])
  }
  return(mean(klo))
}

#| matrix numeric[2]
generate_covariate_quantile_matrix <- function(covariate_data, number_of_bins) {
  # covaraite_data contains on the covaraites excluding the coordiantes or any other auxiliary data
  # Check that the input arguments are valid
  stopifnot(is.data.frame(covariate_data))
  stopifnot(is.numeric(number_of_bins))

  quantiles <- matrix(NA, nrow = (number_of_bins + 1), ncol = ncol(covariate_data))
  j <- 1
  for (i in 1:ncol(covariate_data)) { # note the index start here
    # get a quantile matrix together of the covariates
    cov_range <- max(covariate_data[, i]) - min(covariate_data[, i])
    cov_step <- cov_range / number_of_bins
    quantiles[, j] <- seq(min(covariate_data[, i]), to = max(covariate_data[, i]), by = cov_step)
    j <- j + 1
  }
  return(quantiles)
}


#| matrix numeric[2]
generate_hypercube <- function(covariate_data, number_of_bins) {
  stopifnot(is.data.frame(covariate_data))
  stopifnot(is.numeric(number_of_bins))
  hypercube <- matrix(1, nrow = number_of_bins, ncol = ncol(covariate_data))
  for (i in 1:nrow(covariate_data)) { # the number of pixels
    for (j in 1:ncol(covariate_data)) { # for each column/covariate
      dd <- covariate_data[i, j]
      for (k in 1:number_of_bins) {  # for each quantile
        kl <- covarite_quantiles[k, j]
        ku <- covarite_quantiles[k + 1, j]
        if (dd >= kl & dd <= ku) {
          hypercube[k, j] <- hypercube[k, j] + 1
        }
      }
    }
  }
  return(hypercube)
}