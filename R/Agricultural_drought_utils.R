################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: calculate agricultural drought indices
################################################################################
##
## Creation date: 2022-08-11
##
################################################################################
################################################################################


#' Utils function to calculate the Vegetation Condition Index
#'
#' @param x Numerical vector.
#'
#' @return Numerical vector with the corresponding to the VCI
#' @export
#'
#' @examples
.vci <- function(x, ...){
  
  # Applying VCI formula
  vci <- ( x - min(x, ...) ) / ( max(x, ...) - min(x, ...) ) * 100
  
  # Replace with NAs over areas with full NA vectors
  if(length(vci) != length(x))
    vci <- rep(NA, length(x))
  
  return(vci)
  
}

#' Utils function to calculate the Temperature Condition Index
#'
#' @param x Numerical vector.
#'
#' @return Numerical vector with the corresponding to the TCI
#' @export
#'
#' @examples
.tci <- function(x, ...){
  
  # Applying VCI formula
  tci <- ( max(x, ...) - x ) / ( max(x, ...) - min(x, ...) ) * 100
  
  # Replace with NAs over areas with full NA vectors
  if(length(tci) != length(x))
    tci <- rep(NA, length(x))
  
  return(tci)
  
}



#' Utils function to calculate the Empirical Standardised Soil Moisture Index (ESSMI)
#'
#' @param x Numerical vector.
#' @param dates Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_essmi' function.
#' @param scale Integer value that represents the time scale at which the ESSMI will be computed.
#' @param ref_start optional value that represents the starting point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.
#' @param ref_end Optional value that represents the ending point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.
#' @param distribution Optional value indicating the name of the distribution function to be used in the Kernel Density Estimation
#'  (one of 'Gaussian', 'Rectangular', 'Triangular', 'Epanechnikov', 'Biweight', 'Cosine' and 'Optcosine'). Defaults to 'Gaussian' for ESSMI.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such that this is the standard deviation of the smoothing kernel. 
#' The default is set to the methods of Sheather & Jones (1991) to select the bandwidth using pilot estimation of derivatives.
#' Please see 'Bandwidth Selectors for Kernel Density Estimation' from 'stats' for more information.
#' @param ... Additional variables that can be used for the 'density' function.
#' @return Numerical vector with the corresponding to the ESSMI
#'
#' @examples
.essmi <- function(x, 
                   dates,
                   scale, 
                   ref_start,
                   ref_end,
                   distribution,
                   bw,
                   ...){  ### pass variables to density
  
  # Checking distribution
  if(!distribution %in% c('Gaussian', 'Rectangular', 'Triangular', 'Epanechnikov', 'Biweight', 'Cosine', 'Optcosine'))
    stop("The accepted distributions are 'Gaussian', 'Rectangular', 'Triangular', 'Epanechnikov', 'Biweight', 'Cosine', or 'Optcosine'")
  
  # Changing the distribution name to lower case
  distribution <- tolower(distribution)
  
  # Applying the accumulation according to scale
  x_zoo <- zoo::zoo(x, dates)
  x_zoo <- zoo::rollapply(x_zoo, scale, fill = NA, FUN = sum, align = "right")
  
  # Subsetting 'x' if 'ref_start' or 'ref_end' are given
  pos_ini <- 1
  pos_fin <- length(dates)
  
  if(!is.null(ref_start))
    pos_ini <- grep(ref_start, dates)
  
  if(!is.null(ref_end))
    pos_fin <- grep(ref_end, dates)

  # Substracting the values for the reference period
  x_ref <- x_zoo[pos_ini:pos_fin]
  
  # Extracting NA values
  x_ref <- x_ref[-which(is.na(x_ref))]
  
  # Conditional: if all values are NAs, return NAs (masked regions or oceans)
  nas <- which(is.na(x_ref))
  if(length(nas) > length(x_ref) - 2){
    
    essmi <- rep(NA, length(dates))
    
  } else {
    
    # Computing the Kernel Density Estimation
    dens <- density(as.numeric(x_ref), bw = bw, kernel = distribution, ...)
    
    # Creating the function to calculate the CDF 
    cdf_fun <- spatstat.core::CDF(dens)
    
    # Applying the function to all SM values
    cdf_vals <- cdf_fun(as.numeric(x_zoo))
    
    # Storing the standardized values
    essmi <- qnorm(cdf_vals)
    
  }
  
  return(essmi)
  
}


  
  
  
  
  
  
  
  
