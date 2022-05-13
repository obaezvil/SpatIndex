################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: to calculate Climdex indices per year
################################################################################
##
## Creation date: 2022-04-20
## Updated: 2022-05-13
##
################################################################################
################################################################################

#' Maximum 1-day precipitation
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format '\%Y-\%m-\%d'.
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format '\%Y-\%m'.
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'.
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the maximum 1-day precipitation for each year of analysis.
#' @export
#'
#' @examples
Rx1day <- function(rst_path, 
                   vct = NULL, 
                   start_date,
                   end_date,
                   init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.rx1day'
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.rx1day",
                      index_args = list())
  
  return(res)
  
}

#' Maximum consecutive 5-day precipitation
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the maximum 5-day precipitation for each year of analysis.
#' @export
#'
#' @examples
Rx5day <- function(rst_path, 
                   vct = NULL, 
                   start_date,
                   end_date,
                   init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.rx5day'
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.rx5day",
                      index_args = list())
  
  return(res)
  
}


#' Simple precipitation intensity index
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the simple precipitation intensity index for each year of analysis.
#' @export
#'
#' @examples
sdii <- function(rst_path, 
                   vct = NULL, 
                   start_date,
                   end_date,
                   init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.sdii'
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.sdii",
                      index_args = list())
  
  return(res)
  
}


#' Maximum length of dry spell: maximum number of consecutive days with RR < 1mm
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the maximum length of dry spell (maximum number of consecutive days with RR < 1mm) for each year of analysis.
#' @export
#'
#' @examples
CDD <- function(rst_path, 
                vct = NULL, 
                start_date,
                end_date,
                init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.cdd'
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.cdd",
                      index_args = list())
  
  return(res)
  
}


#' Maximum length of wet spell: maximum number of consecutive days with RR ≥ 1mm
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the Maximum length of wet spell (maximum number of consecutive days with RR ≥ 1mm) for each year of analysis.
#' @export
#'
#' @examples
CWD <- function(rst_path, 
                vct = NULL, 
                start_date,
                end_date,
                init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.cwd'
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.cwd",
                      index_args = list())
  
  return(res)
  
}

#' Annual count of days when PRCP ≥ 10mm
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the annual count of days when PRCP ≥ 10mm for each year of analysis.
#' @export
#'
#' @examples
R10mm <- function(rst_path, 
                  vct = NULL, 
                  start_date,
                  end_date,
                  init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.rnnmm' with a threshold of 10mm
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.rnnmm",
                      index_args = list(thres = 10))
  
  return(res)
  
}

#' Annual count of days when PRCP ≥ 20mm
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the annual count of days when PRCP ≥ 20mm for each year of analysis.
#' @export
#'
#' @examples
R20mm <- function(rst_path, 
                  vct = NULL, 
                  start_date,
                  end_date,
                  init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.rnnmm' with a threshold of 20mm
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.rnnmm",
                      index_args = list(thres = 20))
  
  return(res)
  
}


#' Annual count of days when PRCP ≥ nn mm
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#' @param thres Threshold value that corresponds to nn.
#' 
#' @return raster layers with the spatial distribution of the annual count of days when PRCP ≥ nn mm (where nn is a user selected threshold) for each year of analysis.
#' @export
#'
#' @examples
Rnnmm <- function(rst_path, 
                  vct = NULL, 
                  start_date,
                  end_date,
                  init_month = NULL,
                  thres){
  
  # Applying the '.apply_index' function for '.idx.rnnmm' with a threshold of nn mm
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.rnnmm",
                      index_args = list(thres = thres))
  
  return(res)
  
}


#' Annual total PRCP when RR > 95th percentile
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the annual total PRCP when RR > 95th percentile for each year of analysis.
#' @export
#'
#' @examples
R95p <- function(rst_path, 
                 vct = NULL, 
                 start_date,
                 end_date,
                 init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.rnnp' with a threshold of 0.95
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.rnnp",
                      index_args = list(thres = 0.95))
  
  return(res)
  
}


#' Annual total PRCP when RR > 99th percentile
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the annual total PRCP when RR > 99th percentile for each year of analysis.
#' @export
#'
#' @examples
R99p <- function(rst_path, 
                 vct = NULL, 
                 start_date,
                 end_date,
                 init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.rnnp' with a threshold of 0.99
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.rnnp",
                      index_args = list(thres = 0.99))
  
  return(res)
  
}


#' Contribution to total precipitation from very wet days (95p)
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the contribution to total precipitation from very wet days (95p) for each year of analysis.
#' @export
#'
#' @examples
R95pTOT <- function(rst_path, 
                    vct = NULL, 
                    start_date,
                    end_date,
                    init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.rnnptot' with a threshold of 0.95
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.rnnptot",
                      index_args = list(thres = 0.95))
  
  return(res)
  
}


#' Contribution to total precipitation from very wet days (99p)
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the contribution to total precipitation from very wet days (99p) for each year of analysis.
#' @export
#'
#' @examples
R99pTOT <- function(rst_path, 
                    vct = NULL, 
                    start_date,
                    end_date,
                    init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.rnnptot' with a threshold of 0.99
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.rnnptot",
                      index_args = list(thres = 0.99))
  
  return(res)
  
}


#' Annual total precipitation on wet days
#'
#' @param rst_path File path to either:
#' \itemize{
#'  \item'GTiff files': Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format  '\%Y-\%m-\%d'. 
#'  \item'NCDF files': Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format  '\%Y-\%m'. 
#'}
#'
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL'. 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#'
#' @return raster layers with the spatial distribution of the contribution to total precipitation from very wet days (99p) for each year of analysis.
#' @export
#'
#' @examples
PRCPTOT <- function(rst_path, 
                    vct = NULL, 
                    start_date,
                    end_date,
                    init_month = NULL){
  
  # Applying the '.apply_index' function for '.idx.prcptot'
  res <- .apply_index(rst_path = rst_path, 
                      vct = vct, 
                      start_date = start_date,
                      end_date = end_date,
                      init_month = init_month,
                      index_fun = ".idx.prcptot",
                      index_args = list())
  
  return(res)
  
}
