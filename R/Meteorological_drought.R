################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: calculate meteorological drought indices
################################################################################
##
## Creation date: 2022-06-12
##
################################################################################
################################################################################

#' Wrapper function from the SPEI package to calculate the spatially-distributed Standardized Precipitation Index (SPI)
#'
#' @param P_data 'SpatRaster' object that contains spatially-distributed monthly precipitation data that will be used to calculate the SPI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param scale Integer value that represents the time scale at which the SPI will be computed.
#' @param ref_start optional value that represents the starting point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.
#' @param ref_end Optional value that represents the ending point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.
#' @param distribution Optional value indicating the name of the distribution function to be used for computing the SPI 
#'  (one of 'log-Logistic', 'Gamma' and 'PearsonIII'). Defaults to 'Gamma' for SPI.
#' @param fit Optional value indicating the name of the method used for computing the distribution function parameters 
#'  (one of 'ub-pwm', 'pp-pwm' and 'max-lik'). Defaults to 'ub-pwm'.
#' @param na.rm Should the NA values be removed? Set to TRUE.
#' @param package Either 'SCI' or 'SPEI'. Should the SCI or SPEI package be used in the implementation?
#' @param ... Additional variables that can be used for the 'spi' function of the SPEI package.
#'
#' @return Spatially-distributed SPI values.
#' @export
#'
#' @examples
spatial_spi <- function(P_data,
                        scale,
                        ref_start = NULL,
                        ref_end = NULL,
                        distribution = "Gamma",
                        fit = "ub-pwm",
                        na.rm = TRUE,
                        package = "SCI",
                        ...){
  
  # Check P_data
  if(class(P_data) != "SpatRaster")
    stop("The object 'P_data' must be a SpatRaster")
  
  # Check scale
  if(!class(scale) %in% c("integer", "numeric"))
    stop("The object 'scale' must be a integer that represents the time scale at which the SPI will be computed")
  
  # Check ref_start object
  if(class(ref_start) != "character" & !is.null(ref_start))
    stop("If object 'ref_start' is not set to NULL, it must be a character object that indicates the starting point of the reference period used for computing the SPI. The format should be '%Y-%m'")
  
  # Check ref_end object
  if(class(ref_end) != "character" & !is.null(ref_end))
    stop("If object 'ref_end' is not set to NULL, it must be a character object that indicates the ending point of the reference period used for computing the SPI. The format should be '%Y-%m'")
  
  # Check ref_start and ref_end
  if(is.null(ref_start) & !is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  if(!is.null(ref_start) & is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  # Checking the 'package' object
  if(!package %in% c("SPEI", "SCI"))
    stop("The 'package' object must be either 'SPEI' or 'SCI'")
  
  # Extract dates from object
  dates <- terra::time(P_data)
  
  # Apply SPI
  if(package == "SPEI"){
    
    idx <- terra::app(P_data, .spi.spei, scale = scale, dates = dates, distribution = distribution, fit = fit,
                      ref_start = ref_start, ref_end = ref_end, na.rm = na.rm, ...)
    
  } else {
    
    idx <- terra::app(P_data, .spi.sci, scale = scale, dates = dates, distribution = distribution, fit = fit,
                      ref_start = ref_start, ref_end = ref_end, na.rm = na.rm, ...)
    
  }
  
  
  ## set dates and return
  names(idx)        <- paste0(substr(dates, 1, 7)) 
  terra::time(idx)  <- dates
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  idx[is.infinite(idx)] <- NA
  
  return(idx)
  
}

#' Wrapper function from the SPEI package to calculate the spatially-distributed Standardized Precipitation-Evapotranspiration Index (SPEI)
#'
#' @param P_data 'SpatRaster' object that contains spatially-distributed monthly precipitation data that will be used to calculate the SPEI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param PE_data 'SpatRaster' object that contains spatially-distributed monthly potential evaporation (also referred to as potential evapotranspiation) 
#'  data that will be used to calculate the SPEI. This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. 
#'  They can be set with the function time of the terra package.
#' @param scale Integer value that represents the time scale at which the SPI will be computed.
#' @param ref_start Optional value that represents the starting point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.
#' @param ref_end Optional value that represents the ending point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.
#' @param distribution Optional value indicating the name of the distribution function to be used for computing the SPEI 
#'  (one of 'log-Logistic', 'Gamma' and 'PearsonIII'). Defaults to 'log-Logistic' for SPEI.
#' @param fit Optional value indicating the name of the method used for computing the distribution function parameters 
#'  (one of 'ub-pwm', 'pp-pwm' and 'max-lik'). Defaults to 'ub-pwm'.
#' @param na.rm Should the NA values be removed? Set to TRUE.
#' #' @param package Either 'SCI' or 'SPEI'. Should the SCI or SPEI package be used in the implementation?
#' @param ... Additional variables that can be used for the 'spi' function of the SPEI package.
#'
#' @return Spatially-distributed SPEI values.
#' @export
#'
#' @examples
spatial_spei <- function(P_data,
                         PE_data,
                        scale,
                        ref_start = NULL,
                        ref_end = NULL,
                        distribution = "log-Logistic",
                        fit = "ub-pwm",
                        na.rm = TRUE,
                        package = "SCI",
                        ...){
  
  # Check P_data
  if(class(P_data) != "SpatRaster")
    stop("The object 'P_data' must be a SpatRaster")
  
  # Check scale
  if(!class(scale) %in% c("integer", "numeric"))
    stop("The object 'scale' must be a integer that represents the time scale at which the SPEI will be computed")
  
  # Check ref_start object
  if(class(ref_start) != "character" & !is.null(ref_start))
    stop("If object 'ref_start' is not set to NULL, it must be a character object that indicates the starting point of the reference period used for computing the SPEI. The format should be '%Y-%m'")
  
  # Check ref_end object
  if(class(ref_end) != "character" & !is.null(ref_end))
    stop("If object 'ref_end' is not set to NULL, it must be a character object that indicates the ending point of the reference period used for computing the SPEI. The format should be '%Y-%m'")
  
  # Check ref_start and ref_end
  if(is.null(ref_start) & !is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  if(!is.null(ref_start) & is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  # Extract dates from P_data object
  dates_p <- terra::time(P_data)
  
  # Extract dates from PE_data object
  dates_pe <- terra::time(PE_data)
  
  # Checking temporal dimension
  if(length(dates_p) != length(dates_pe))
    stop("The lengths of 'P_data' and 'PE_data' differ")
  
  # Checking that P and PE have the same period
  if(!all(dates_p %in% dates_pe))
    stop("'P_data' and 'PE_data' have different periods")
    
  # Check extent of products to resample
  if(!terra::compareGeom(P_data, PE_data)){
    # P higher resolution that PE
    if(terra::res(P_data)[1] < terra::res(PE_data)[1]){
      PE_data <- terra::resample(PE_data, P_data, method = "near")
    } else {
      P_data <- terra::resample(P_data, PE_data, method = "near")
    }
  }
  
  # Checking the 'package' object
  if(!package %in% c("SPEI", "SCI"))
    stop("The 'package' object must be either 'SPEI' or 'SCI'")
  
  # Compute P-PE
  diff <- P_data - PE_data
  
  # Apply SPEI
  if(package == "SPEI"){
    
    idx <- terra::app(diff, .spei.spei, scale = scale, dates = dates_p, distribution = distribution, fit = fit,
                      ref_start = ref_start, ref_end = ref_end, na.rm = na.rm, ...)
    
  } else {
    
    idx <- terra::app(diff, .spei.sci, scale = scale, dates = dates_p, distribution = distribution, fit = fit,
                      ref_start = ref_start, ref_end = ref_end, na.rm = na.rm, ...)
    
  }
  
  
  ## set dates and return
  names(idx)        <- paste0(substr(dates_p, 1, 7)) 
  terra::time(idx)  <- dates_p
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  idx[is.infinite(idx)] <- NA
  
  return(idx)
  
}





