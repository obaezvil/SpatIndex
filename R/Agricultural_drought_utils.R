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
  
  # Removing negative NDVI values
  x[x < 0] <- NA
  
  # Applying VCI formula
  vci <- ( x - min(x, ...) ) / ( max(x, ...) - min(x, ...) ) * 100
  
  # Replace with NAs over areas with full NA vectors
  if(length(vci) != length(x))
    vci <- rep(NA, length(x))
  
  return(vci)
  
}

#' Utils function to calculate the Vegetation Condition Index per period
#'
#' @param x Numerical vector.
#' @param dates Vector of dates passed from the parent function.
#'
#' @return Numerical vector with the corresponding to the VCI
#' @export
#'
#' @examples
.vci_period <- function(x, dates, ...){
  
  # Removing negative NDVI values
  x[x < 0] <- NA
  
  # Convert the value in a zoo object
  x_zoo <- zoo::zoo(x, dates)
  
  # Getting periods
  period    <- substr(dates, 5, 10)
  u_periods <- unique(period)
  
  # Creating empty vector
  vci <- c()
  
  for(i in 1:length(u_periods)){
    
    # Storing position vector for each period
    pos <- which(period %in% u_periods[i])
    
    # Extracting values for period and applying the VCI
    period_vals <- x_zoo[pos]
    vci[pos]    <- ( period_vals - min(period_vals, ...) ) / ( max(period_vals, ...) - min(period_vals, ...) ) * 100
    
  }

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

#' Utils function to calculate the Temperature Condition Index per period
#'
#' @param x Numerical vector.
#' @param dates Vector of dates passed from the parent function.
#'
#' @return Numerical vector with the corresponding to the TCI
#' @export
#'
#' @examples
.tci_period <- function(x, dates, ...){
  
  # Convert the value in a zoo object
  x_zoo <- zoo::zoo(x, dates)
  
  # Getting periods
  period    <- substr(dates, 5, 10)
  u_periods <- unique(period)
  
  # Creating empty vector
  tci <- c()
  
  for(i in 1:length(u_periods)){
    
    # Storing position vector for each period
    pos <- which(period %in% u_periods[i])
    
    # Extracting values for period and applying the VCI
    period_vals <- x_zoo[pos]
    tci[pos]    <- ( max(period_vals, ...) - period_vals ) / ( max(period_vals, ...) - min(period_vals, ...) ) * 100
    
  }
  
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
#' @param missing_ratio Ratio of missing data that is acceptable for the computation. Set to 0.2 by default (20\%).
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
                   missing_ratio,
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
  
  # Obtaining the resulting unique steps
  steps   <- substr(index(x_ref), 6, 10)
  u_steps <- unique(steps)
  
  # Calculating NAs in period
  nas <- length(which(is.na(x_ref)))
  
  if(length(nas) < 1)
    nas <- 0
  
  missing <- nas / length(x_ref)
  
  # Calculating unique values
  uniq_vals <- length(unique(x_ref)) / length(x_ref)
  
  # Conditional: if all values are NAs, return NAs (masked regions or oceans)
  if(missing_ratio < missing){
    
    essmi <- rep(NA, length(dates))
    
  } else {
  
    # Creating an object to store the resulting data
    essmi <- c()
    
    #######
    ####### Loop to evaluate each one of the time steps
    #######
  
    for(i in 1:length(u_steps)){
      
      # Extracting the data according to step 'i'
      pos_step   <- which(steps %in% u_steps[i])
      x_ref_step <- x_ref[pos_step]
      
      # Substracting NA values
      na_pos <- which(is.na(x_ref_step))
      
      if(length(na_pos) > 0)
        x_ref_step <- x_ref_step[-which(is.na(x_ref_step))]
      
      # Calculating unique values
      uniq_vals <- length(unique(round(x_ref_step, 3)))
      
      # Obtaining all values for that time step
      pos_dates   <- grep(u_steps[i], dates)
      values_step <- x_zoo[pos_dates]
      
      if(uniq_vals > 4){
          
        # Computing the Kernel Density Estimation
        dens <- density(as.numeric(x_ref_step), bw = bw, kernel = distribution, ...)
     
        # Creating the function to calculate the CDF 
        cdf_fun <- spatstat.core::CDF(dens)
        
        # Applying the function to all SM values for time step
        cdf_vals <- cdf_fun(as.numeric(values_step))
        
        # Storing the standardized values
        res <- qnorm(cdf_vals)
        res[which(is.infinite(res))] <- NA
        res[which(is.nan(res))]      <- NA
        essmi[pos_dates] <- res
        
      } else {
        
        essmi[pos_dates] <- NA
        
      } # End if else unique values
      
    } # End for
       
  } # End if else missing data
  
  
  return(essmi)

}

#' Utils function to calculate the z-score
#'
#' @param x Numerical vector.
#' @param dates Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_zscore' function.
#' @return Numerical vector with the corresponding to the zFAPAR
#'
#' @examples
.zscore <- function(x,
                    dates){
  
  # Applying the accumulation according to scale
  x_zoo <- zoo::zoo(x, dates)
  
  # Obtaining the resulting unique steps
  steps   <- substr(dates, 6, 10)
  u_steps <- unique(steps)
  
  # Conditional: if all values are NAs, return NAs (masked regions or oceans)
  if(length(x_zoo) == length(which(is.na(x_zoo)))){
    
    zsc <- rep(NA, length(x_zoo))
    
  } else {
    
    # Creating an object to store the resulting data
    zsc <- c()
    
    #######
    ####### Loop to evaluate each one of the time steps
    #######
    
    for(i in 1:length(u_steps)){
      
      # Extracting the data according to step 'i'
      pos_step   <- which(steps %in% u_steps[i])
      x_step     <- x_zoo[pos_step]
      
      # Calculating components of z-score
      x_mean <- mean(x_step, na.rm = TRUE)
      x_sd   <- sd(x_step, na.rm = TRUE)
      
      # Calculating z-score for step 'i'
      zsc[pos_step] <- (x_step - x_mean) / x_sd
      
    } # end for
    
    
  } # end else
  
  # Returning resulting object
  return(zsc)
  
}
