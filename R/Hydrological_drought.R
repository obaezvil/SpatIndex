################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: calculate hydrological drought indices
################################################################################
##
## Creation date: 2022-08-23
##
################################################################################
################################################################################

#' Wrapper function from the SCI package to calculate the Standardized Streamflow Index (SSI)
#'
#' @param x Zoo object with one or multiple monthly time series.
#' @param scale Integer value that represents the time scale at which the SSI will be computed.
#' @param ref_start optional value that represents the starting point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.
#' @param ref_end Optional value that represents the ending point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.
#' @param distribution Optional value indicating the name of the distribution function to be used for computing the SSI 
#'  (one of 'log-Logistic', 'GEV', 'Gamma', and 'PearsonIII').
#'  #'
#' @return Standardised Streamflow Index.
#' @export
#'
#' @examples
ssi <- function(x, 
                scale, 
                ref_start = NULL, 
                ref_end = NULL, 
                distribution){
  
  # Check zoo data
  if(!zoo::is.zoo(x))
    stop("x must be a zoo object!")
  
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
  
  # Extract dates from the x object
  dates <- zoo::index(x)
  
  # Apply the SSI
  idx <- .ssi(x, scale = scale, ref_start = ref_start, ref_end = ref_end, distribution = distribution)

  return(idx)
  
}



#' Snow Water Equivalent Index (SWEI)
#'
#' @param SWE_data 'SpatRaster' object that contains spatially-distributed snow water equivalent data that will be used to calculate the SWEI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param scale Integer value that represents the time scale at which the ESSMI will be computed.
#' @param start_month Optional value that represents the starting point of the period that will be returned from the Index.
#'  The idea is not to return periods where SWE is not present. 
#'  The date should be introduced as '\%m-\%d'. For example: "02-01".
#'  The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.
#' @param end_month Optional value that represents the ending point of the period that will be returned from the Index.
#'  The idea is not to return periods where SWE is not present. 
#'  The date should be introduced as '\%m-\%d'. For example: "02-01".
#'  The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.
#' @param missing_ratio Ratio of missing data that is acceptable for the computation. Set to 0.2 by default (20\%).
#' @param ... Additional variables that can be used for the 'density' function.
#'
#' @return Spatially-distributed SWEI values.
#' @export
#'
#' @examples
spatial_swei <- function(SWE_data,
                         scale,
                         start_month = NULL,
                         end_month = NULL,
                         missing_ratio = 0.2,
                          ...){
  
  # Check SWE_data
  if(class(SWE_data) != "SpatRaster")
    stop("The object 'SWE_data' must be a SpatRaster")
  
  # Check scale
  if(!class(scale) %in% c("integer", "numeric"))
    stop("The object 'scale' must be a integer that represents the time scale at which the SWEI will be computed")
  
  # Check ref_start object
  if(class(start_month) != "character" & !is.null(start_month))
    stop("If object 'start_month' is not set to NULL, it must be a character object that indicates the starting point of the reference period used for computing the SPI. The format should be '%Y-%m'")
  
  # Check ref_end object
  if(class(end_month) != "character" & !is.null(end_month))
    stop("If object 'end_month' is not set to NULL, it must be a character object that indicates the ending point of the reference period used for computing the SPI. The format should be '%Y-%m'")
  
  if(class(start_month) != class(end_month))
    stop("The class of 'start_month' and 'end_month' must be the same (i.e., NULL or character)")
  
  # Check start_month and end_month
  if(is.null(start_month) & !is.null(end_month))
    stop("The objects 'start_month' and 'end_month' should be either both NULL or both character!")
  
  if(!is.null(start_month) & is.null(end_month))
    stop("The objects 'start_month' and 'end_month' should be either both NULL or both character!")
  
  # Extract dates from object
  dates <- terra::time(SWE_data)
  
  # Apply SWEI
  idx <- terra::app(SWE_data, .swei, scale = scale, dates = dates,
                    missing_ratio = missing_ratio)
  
  ## set dates and return
  names(idx)        <- paste0(substr(dates, 1, 7)) 
  terra::time(idx)  <- dates
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  idx[is.infinite(idx)] <- NA
  
  # Subset the selected months in the case that 'start_month' and 'end_month' are not NULL
  if(!is.null(start_month)){
    
    period   <- substr(dates, 6, 12)
    u_period <- unique(period)
    pos_ini  <- grep(start_month, u_period)
    pos_fin  <- grep(end_month, u_period)
    pos_per  <- u_period[pos_ini:pos_fin]
    
    final_pos <- which(period %in% pos_per)
    idx       <- idx[[final_pos]]
    
  }
  
  # Masking regions full with zeroes where the computation of the index does
  #  not make sense
  
  return(idx)
  
}
