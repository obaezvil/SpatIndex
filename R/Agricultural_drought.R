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
## Creation date: 2022-08-11
##
################################################################################
################################################################################

#' Vegetation Condition Index (VCI; Kogan 1995)
#'
#' @param NDVI_data 'SpatRaster' object that contains spatially-distributed NDVI data that will be used to calculate the VCI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @return Spatially-distributed VCI values.
#' @export
#'
#' @examples
spatial_vci <- function(NDVI_data){
  
  # Check P_data
  if(class(NDVI_data) != "SpatRaster")
    stop("The object 'NDVI_data' must be a SpatRaster")
  
  # Extract dates from object
  dates <- terra::time(NDVI_data)
  
  # Apply VCI
  idx <- terra::app(NDVI_data, .vci, na.rm = TRUE)
  
  ## set dates and return
  names(idx)        <- paste0(substr(dates, 1, 7)) 
  terra::time(idx)  <- dates
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  
  return(idx)
  
}

#' Temperature Condition Index (TCI; Kogan 1995)
#'
#' @param BT_data 'SpatRaster' object that contains spatially-distributed BT data that will be used to calculate the TCI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param na.rm Should the NA values be removed? Set to TRUE.
#' @return Spatially-distributed TCI values.
#' @export
#'
#' @examples
spatial_tci <- function(BT_data){
  
  # Check P_data
  if(class(BT_data) != "SpatRaster")
    stop("The object 'BT_data' must be a SpatRaster")
  
  # Extract dates from object
  dates <- terra::time(BT_data)
  
  # Apply TCI
  idx <- terra::app(BT_data, .tci, na.rm = TRUE)
  
  ## set dates and return
  names(idx)        <- paste0(substr(dates, 1, 7)) 
  terra::time(idx)  <- dates
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  
  return(idx)
  
}

#' Vegetation Health Index (VHI; Kogan 1995)
#'
#' @param VCI_data 'SpatRaster' object that contains spatially-distributed VCI data that will be used to calculate the VHI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param TCI_data 'SpatRaster' object that contains spatially-distributed TCI data that will be used to calculate the VHI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#'  @param alpha empirical weight that is applied to the VCI (set to 0.5). Please note that the weight of TCI is 1 - alpha.
#' @return Spatially-distributed VHI values.
#' @export
#'
#' @examples
spatial_vhi <- function(VCI_data, TCI_data, alpha = 0.5, resampling2low = TRUE){
  
  # Check VCI_data
  if(class(VCI_data) != "SpatRaster")
    stop("The object 'VCI_data' must be a SpatRaster")
  
  # Check TCI_data
  if(class(TCI_data) != "SpatRaster")
    stop("The object 'TCI_data' must be a SpatRaster")
  
  # Extract dates from VCI_data object
  dates_vci <- terra::time(VCI_data)
  
  # Extract dates from TCI_data object
  dates_tci <- terra::time(TCI_data)
  
  # Checking temporal dimension
  if(length(dates_vci) != length(dates_tci))
    stop("The lengths of 'VCI_data' and 'TCI_data' differ")
  
  # Checking that VCI and TCI have the same period
  if(!all(dates_vci %in% dates_tci))
    stop("'dates_tci' and 'TCI_data' have different periods")
  
  # Evaluating spatial resolution
  if(any(terra::res(VCI_data) != terra::res(TCI_data))){
    
    if(resampling2low){
      
      # condition if VCI_data has higher resolution
      if(terra::res(VCI_data)[1] < terra::res(TCI_data)[1] | 
         terra::res(VCI_data)[2] < terra::res(TCI_data)[2]){
        
        VCI_data <- terra::resample(VCI_data, TCI_data, method = "near")
        
        # condition if TCI_data has higher resolution
      } else {
        
        TCI_data <- terra::resample(TCI_data, VCI_data, method = "near")
        
      }
      
      # condition if resampling2low = FALSE 
    } else {
      
      # condition if VCI_data has higher resolution
      if(terra::res(VCI_data)[1] < terra::res(TCI_data)[1] | 
         terra::res(VCI_data)[2] < terra::res(TCI_data)[2]){
        
        TCI_data <- terra::resample(TCI_data, VCI_data, method = "bilinear")
        
        # condition if TCI_data has higher resolution
      } else {
        
        VCI_data <- terra::resample(VCI_data, TCI_data, method = "bilinear")
        
      }
      
    }
    
  }
  
  # Apply VHI
  idx <- ( VCI_data * alpha ) + ( (1 - alpha ) * TCI_data)
  
  ## set dates and return
  names(idx)        <- paste0(substr(dates, 1, 7)) 
  terra::time(idx)  <- dates
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  
  return(idx)
  
}


