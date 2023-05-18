# compute KL divergence
############################################################
library(entropy)

compute_kl_divergence <- function(covariate_data, sample_data) {
  # Kullback-Leibler (KL) divergence
  klo <- vector(length = ncol(covariate_data), mode = "numeric")
  for (col in 1:ncol(covariate_data)) {
    klo[col] <- KL.empirical(covariate_data[, col], sample_data[, col])
  }
  return(mean(klo))
}

#| matrix numeric[2]
generate_covariate_quantile_matrix <- function(covariate_data, number_of_bins) {
  # covaraite_data contains on the covaraites excluding the coordiantes or any other auxiliary data
  # Check that the input arguments are valid
  stopifnot(is.data.frame(covariate_data))
  stopifnot(is.numeric(number_of_bins))
  tolerance <- 1e-6  # this little bit of tolerance makes sure the covariate points don't fall outside when compared
  quantiles <- matrix(NA, nrow = (number_of_bins + 1), ncol = ncol(covariate_data))
  for (i in 1:ncol(covariate_data)) {
    # get a quantile matrix together of the covariates
    cov_range <- max(covariate_data[, i]) - min(covariate_data[, i]) + 2 * tolerance
    cov_step <- cov_range / number_of_bins
    quantiles[, i] <- seq(min(covariate_data[, i]) - tolerance, to = max(covariate_data[, i]) + tolerance, by = cov_step)
  }
  return(quantiles)
}


#| matrix numeric[2]
generate_hypercube <- function(covariate_data, quantiles) {
  stopifnot(is.data.frame(covariate_data))
  stopifnot(ncol(covariate_data) == ncol(quantiles))
  number_of_bins <- nrow(quantiles) - 1
  hypercube <- matrix(0, nrow = number_of_bins, ncol = ncol(covariate_data))
  for (i in 1:nrow(covariate_data)) { # the number of pixels
    for (j in 1:ncol(covariate_data)) { # for each column/covariate
      dd <- covariate_data[i, j]
      for (k in 1:number_of_bins) {  # for each quantile
        kl <- quantiles[k, j]
        ku <- quantiles[k + 1, j]
        if (dd >= kl & dd <= ku) {
          hypercube[k, j] <- hypercube[k, j] + 1
        }
      }
    }
  }
  # Replace zeros with a small number so we don't have to deal with zeros
  hypercube[which(hypercube == 0)] <- 0.000001
  return(hypercube)
}

generate_hypercube_vec <- function(covariate_data, quantiles) {
  stopifnot(is.data.frame(covariate_data))
  stopifnot(ncol(covariate_data) == ncol(quantiles))
  number_of_bins <- nrow(quantiles) - 1
  hypercube <- matrix(0, nrow = number_of_bins, ncol = ncol(covariate_data))
  for (i in 1:nrow(covariate_data)) { # the number of pixels
    for (j in 1:ncol(covariate_data)) { # for each column/covariate
      dd <- covariate_data[i, j]
      hypercube[findInterval(dd, quantiles[, j]), j] <- hypercube[findInterval(dd, quantiles[, j]), j] + 1
    }
  }

  # Replace zeros with a small number so we don't have to deal with zeros
  hypercube[which(hypercube == 0)] <- 0.000001
  return(hypercube)
}


composite_from_quantiles <- function(covariate_data, quantiles) {
  stopifnot(is.data.frame(covariate_data))
  stopifnot(ncol(covariate_data) == ncol(quantiles))
  number_of_bins <- nrow(quantiles) - 1
  hypercube <- matrix(0, nrow = nrow(covariate_data), ncol = ncol(covariate_data))
  for (i in 1:nrow(covariate_data)) { # the number of pixels
    for (j in 1:ncol(covariate_data)) { # for each column/covariate
      dd <- covariate_data[i, j]
      hypercube[i, j] <- findInterval(dd, quantiles[, j])
    }
  }
  return(hypercube)
}