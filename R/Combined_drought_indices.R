################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: calculate combined drought indices
################################################################################
##
## Creation date: 2022-11-27
##
################################################################################
################################################################################

#' Combined Drought Indicator (CDI; Cammalleri et al. 2021)
#'
#' @param SPI1_data 'SpatRaster' object that contains spatially-distributed SPI-1 data that will be used to calculate the CDI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param SPI3_data 'SpatRaster' object that contains spatially-distributed SPI-3 data that will be used to calculate the CDI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param zSM 'SpatRaster' object that contains spatially-distributed zSM data that will be used to calculate the CDI (to calculate this index see the 'spatial_zscore' function). 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param zfAPAR 'SpatRaster' object that contains spatially-distributed zfAPAR data that will be used to calculate the CDI (to calculate this index see the 'spatial_zscore' function). 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#'  @param threshold_spi1 Threshold to convert the SPI-1 data to a Boolean 'SpatRaster' (1 values > threshold; 0, otherwise). The default value of this threshold is -2 (Cammalleri et al. 2021).
#'  #'  @param threshold_spi3 Threshold to convert the SPI-3 data to a Boolean 'SpatRaster' (1 values > threshold; 0, otherwise). The default value of this threshold is -1 (Cammalleri et al. 2021).
#' @return Spatially-distributed CDI values.
#' @export
#'
#' @examples
spatial_cdi <- function(SPI1_data, 
                        SPI3_data, 
                        zSM, zfAPAR, 
                        threshold_spi1 = -2, 
                        threshold_spi3 = -1){
  
  # Check SPI1_data
  if(class(SPI1_data) != "SpatRaster")
    stop("The object 'SPI1_data' must be a SpatRaster")
  
  # Check SPI3_data
  if(class(SPI3_data) != "SpatRaster")
    stop("The object 'SPI3_data' must be a SpatRaster")
  
  # Check zSM
  if(class(zSM) != "SpatRaster")
    stop("The object 'zSM' must be a SpatRaster")
  
  # Check zfAPAR
  if(class(zfAPAR) != "SpatRaster")
    stop("The object 'zfAPAR' must be a SpatRaster")
  
  # Check that all spatial objects have the same nlyrs.
  n_layers <- unique(c(terra::nlyr(SPI1_data), terra::nlyr(SPI3_data), 
                terra::nlyr(zSM), terra::nlyr(zfAPAR)))
  if(length(n_layers) > 1)
    stop("The objects 'SPI1_data', 'SPI3_data', 'zSM', and 'zfAPAR' does not have the same number of layers!")
  
  # Check if all datasets have the same raster geometry
  if(!terra::compareGeom(SPI1_data, SPI3_data, zSM))
    stop("The geometry from 'SPI1_data', 'SPI3_data', and 'zSM' is not the same!")
  
  if(!terra::compareGeom(SPI1_data, SPI3_data, zfAPAR))
    stop("The geometry from 'SPI1_data', 'SPI3_data', and 'zfAPAR' is not the same!")
 
  # Extract dates from the SPI1_data object
  dates <- terra::time(SPI1_data)
  
  # Calculating the zSPI using the SPI1_data and SPI3_data
  zSPI     <- .zspi(SPI1_data, SPI3_data, threshold_spi1, threshold_spi3)
  
  # Calculating the 'SpatRaster' that contains the zSPI, zSM, and zfAPAR values 
  #   (see the '.cdi_vals' function from the Combined_drought_indices_utils.R')
  cdi_vals <- .cdi_vals(zSPI, zSM, zfAPAR)
  
  # Apply the CDI
  idx <- terra::app(cdi_vals, .get_cdi)

  ## set dates and return
  names(idx)        <- paste0(substr(dates, 1, 7)) 
  terra::time(idx)  <- dates
  
  # Avoid NaNs
  idx[is.nan(idx)]      <- NA
  
  return(idx)
  
}




