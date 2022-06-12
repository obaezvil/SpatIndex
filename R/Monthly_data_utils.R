################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: import data from file paths in the required fromat to apply 
##            monthly indices
################################################################################
##
## Creation date: 2022-06-09
##
################################################################################
################################################################################


#' Import a monthly SpatRaster from monthly NCDF or GeoTIFF files.
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Monthly raster files that include the date in the file names (e.g., 'ProductX_1989-02.tif'). The dates must have the format '\%Y-\%m'.
#'  \item'NCDF files': Annual raster files that contain the monthly layers respective for that year. The file name must include the date (e.g., 'ProductX_1989.nc'). The dates must have the format '\%Y'.
#'}
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 16 for 'ProductX_1989-02.tif', and end_date = 13 for 'ProductX_1989.nc'.
#' @param years Integer vector indicating the years that should be imported. For example, years <- 2000:2010 for importing monthly data from 2000 to 2010.
#'
#' @return SpatRaster of the monthly data that will be used in further analysis.
#' @export
#'
#' @examples
import_monthly_data <- function(rst_path,
                                start_date,
                                end_date,
                                years = NULL
){
  
  # Checking the 'years' element
  if(!class(years) == "integer" & !is.null(years))
    stop("If object 'years' is not set to NULL, it must be an integer vector that contains all the years that must be considered!")
  
  # Getting product and dates
  res     <- .get_names(rst_path, start_date, end_date)
  product <- res$product
  dates   <- res$dates
  ext     <- res$format
  
  # Checking dates
  n <- nchar(dates)
  n <- unique(n)
  if(length(n) > 1)
    stop("The dates of the files have an ambiguous format. If you are importing data from tif files,
           the dates should have the format '%Y-%m', while for nc files, the format should be '%Y'.")
  
  # Setting dates with the full format
  if(ext == "tif"){
    
    dates <- paste0(dates, "-01")
    
  } else if(ext %in% c("nc", "nc4")){
    
    s     <- paste0(dates[1], "-01-01")
    e     <- paste0(dates[length(dates)], "-12-01")
    dates <- hydroTSM::mip(as.Date(s), as.Date(e))
    
  }
  
  ## Importing gridded data
  product <- terra::rast(product)
  
  # Subsetting the dates according to selection
  if(!is.null(years)){
    
    y       <- substr(dates, 1, 4)
    pos     <-which(as.numeric(y) %in% years)
    dates   <- dates[pos]
    product <- product[[pos]]
    
  }
  
  # Setting dates to SpatRaster
  terra::time(product) <-as.Date(dates)

  return(product)
  
}


