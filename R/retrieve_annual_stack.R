
# rst_path only tif or nc, tif daily and nc monthly with daily layers

#' Retrieve SpatRaster for a selected year
#'
#' @param rst.path File path to either:
#' #' \itemize{
#'   \item{'GTiff files': }{Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format '%Y-%m-%d'.}
#'   \item{'NCDF files': }{Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format '%Y-%m'.}
#'}
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param selected_year The year that will be used to construct the 'SpatRaster' that will be returned by the function.
#'
#' @return 'SpatRaster' object for the designed year according to the initial mont that has been selected.
#'
#' @examples
.get_yearly <- function(rst_path,
                        init_month = NULL,
                        start_date,
                        end_date,
                        selected_year
){
  
  # Checking the init_month parameter
  if(!is.null(init_month)){
    
    if(!is.numeric(init_month))
      stop("The parameter 'init_month' should be an integer going from 1 to 12 according to the month that should be used as the start of the year.")
    
    if(init_month < 0 | init_month > 12)
      stop("The parameter 'init_month' should be an integer going from 1 to 12 according to the month that should be used as the start of the year.")
    
  }

  # Getting product and dates
  res     <- .get_names(rst_path, start_date, end_date)
  product <- res$product
  dates   <- res$dates
  ext     <- res$format
  
  # Subset values for tif files (single files)
  if(ext == "tif"){
    
    if(is.null(init_month)){
      
      # Subsetting product and importing it
      pos      <- grep(selected_year, dates)
      selected <- product[pos]
      selected <- terra::rast(selected)
      terra::time(selected) <- as.Date(dates[pos])
      
    } else {
      
      # Subsetting dates according to 'init_month'
      months_ex <- sprintf("%02d", 1:12)
      months    <- substr(dates, 1, 7)
      ini       <- paste0(selected_year, "-", months_ex[init_month])
      fin       <- paste0(selected_year + 1, "-", months_ex[init_month - 1])
      ini.pos   <- min(grep(ini, dates))
      fin.pos   <- max(grep(fin, dates))
      
      # Subsetting product
      selected <- product[ini.pos:fin.pos]
      selected <- terra::rast(selected)
      terra::time(selected) <- as.Date(dates[ini.pos:fin.pos])
    }

  } else if(ext %in% c("nc", "nc4")){
    
    if(is.null(init_month)){
      
      # Subsetting product and importing it
      pos      <- grep(selected_year, dates)
      selected <- product[pos]
      selected <- terra::rast(selected)
      ini.dts  <- paste0(dates[pos[1]], "-01")
      fin.dts  <- paste0(selected_year, "-12-31")
      dates    <- hydroTSM::dip(ini.dts, fin.dts)
      terra::time(selected) <- as.Date(dates)
      
    } else {
      
      months_ex <- sprintf("%02d", 1:12)
      months    <- substr(dates, 1, 7)
      ini       <- paste0(selected_year, "-", months_ex[init_month])
      fin       <- paste0(selected_year + 1, "-", months_ex[init_month - 1])
      ini.pos   <- grep(ini, dates)
      fin.pos   <- grep(fin, dates)
      
      # Subsetting product
      selected <- product[ini.pos:fin.pos]
      selected <- terra::rast(selected)
      ini.dts  <- paste0(ini, "-01")
      fin.dts  <- paste0(selected_year + 2, "-12-31")
      dates    <- hydroTSM::dip(ini.dts, fin.dts)[1:terra::nlyr(selected)]
      terra::time(selected) <- as.Date(dates)
      
    } # end if is.null(init_month)
    
  } # end if format files
  
  return(selected)
  
} # end .get_yearly



#' Get layer names, extension, and dates
#'
#' @param rst.path File path to either:
#' #' \itemize{
#'   \item{'GTiff files': }{Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format '%Y-%m-%d'.}
#'   \item{'NCDF files': }{Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format '%Y-%m'.}
#'}
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#'
#' @return this function retrieves the file names (including the path), their corresponding dates and file extension.
#'
#' @examples
.get_names <- function(rst_path,
                       start_date,
                       end_date
){
  
  # Controling the rst_path parameter
  if(!is.character(rst_path))
    stop("The parameter 'rst_path' must be a character indicating the path of the stored files!")
  
  # reading the data
  product <- list.files(rst_path, full.names = TRUE)
  
  # checking availability of files
  if(length(product) < 1)
    stop("There are no files inside the 'rst_path' directory!")
  
  # Extracting extension of the files
  ext <- tools::file_ext(product)
  
  # Checking file format
  if(!any(ext %in% c("tif", "nc", "nc4")))
    stop("The files in folder: ", rst_path, " must be TIFF or NCDF!")
  
  # Assessing whether we have one single extension
  if(length(unique(ext)) > 1 & "tif" %in% unique(ext))
    product <- list.files(rst_path, full.names = TRUE, pattern = ".tif$")
  
  if(length(unique(ext)) > 1 & "nc" %in% unique(ext))
    product <- list.files(rst_path, full.names = TRUE, pattern = ".nc$")
  
  if(length(unique(ext)) > 1 & "nc4" %in% unique(ext))
    product <- list.files(rst_path, full.names = TRUE, pattern = ".nc4$")
  
  # Extracting final extension of the files
  ext <- unique(tools::file_ext(product))
  
  # Substracting the dates from filenames
  dates <- substr(basename(product), start_date, end_date)
  
  # Creating a list to store results
  res <- list(product = product,
              dates = dates,
              format = ext)
  
  
  return(res)

}

