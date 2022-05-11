################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: base functions to calculate Climdex indices per year
################################################################################
##
## Date: 2020-04-20
##
################################################################################
################################################################################

#' Maximum 1-day precipitation
#'
#' @param rst.path File path to either:
#' #' \itemize{
#'   \item{'GTiff files': }{Daily raster files that include the date in the file names (e.g., 'ProductX_1989-02-28.tif'). The dates must have the format '%Y-%m-%d'.}
#'   \item{'NCDF files': }{Monthly raster files that contain the daily layers respective for that month. The file name must include the date (e.g., 'ProductX_1989-02.nc'). The dates must have the format '%Y-%m'.}
#'}
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param start_date Position where the dates of the raster files start. For example, start_date = 10 for 'ProductX_1989-02-28.tif' as the date start in the 10th character.
#' @param end_date Position where the dates of the raster files end. For example, end_date = 19 for 'ProductX_1989-02-28.tif', and end_date = 16 for 'ProductX_1989-02.nc'.
#' @param init_month Numeric value that represents the initial month that will be used in the calculation of the indices for each year (1 = Jan, 2 = Feb, ...., 12 = December). 
#' This parameter is set to NULL meaning that if not specified, the year will start in January.
#' @param index_fun The function that will be applied to the raster layers. 
#' @param index_args List of additional arguments to be passed to the selected function ('index_fun'). 
#'
#' @return
#'
#' @examples
.apply_index <- function(rst_path, 
                          vct = NULL, 
                          start_date,
                          end_date,
                          init_month = NULL,
                          index_fun,
                          index_args = list()){
  
  # Getting product and dates
  res     <- .get_names(rst_path, start_date, end_date)
  product <- res$product
  dates   <- res$dates
  format  <- res$format
  
  # Matching the function
  index_fun <- match.fun(index_fun)
  
  # Setting additional parameters
  index_args_default <- formals(index_fun)
  index_args         <- modifyList(index_args_default, index_args)
  
  # Subsetting the years
  years   <- substr(dates, 1, 4)
  u_years <- as.numeric(unique(years))
    
  # Creating an empty list to store the results
  results <- list()
  
  # Checkinf if init_month == NULL
  if(is.null(init_month)){
    
    for(i in 1:length(u_years)){
      
      rst <- .get_yearly(rst_path = rst_path, init_month = init_month,
                          start_date = start_date, end_date = end_date,
                          selected_year = u_years[i])
      
      # Applying the function
      index_args_it  <- modifyList(index_args, list(rst = rst)) 
      r_ind        <- do.call(index_fun, as.list(index_args_it))
      results[[i]] <- r_ind
      
    } # end for
    
    
  } else {
    
    for(i in 1:(length(u_years) - 1)){
      
      rst <- .get_yearly(rst_path = rst_path, init_month = init_month,
                         start_date = start_date, end_date = end_date,
                         selected_year = u_years[i])
      
      # Applying the function
      index_args_it  <- modifyList(index_args, list(rst = rst)) 
      r_ind        <- do.call(index_fun, as.list(index_args_it))
      results[[i]] <- r_ind
      
    } # end for
    
    
  } # end if is.null(init_month)
  
  results        <- terra::rast(results)
  names(results) <- u_years
  
  # cropping to extent if given
  if(!is.null(vct)){
    rst <- terra::crop(rst, vct, snap = "out")
  }
  
  return(results)
  
}


###################################
#' RX1day
#'
#' @param rst 'SpatRaster' object with daily precipitation values.
#'
#' @return Maximum 1-day precipitation
#'
#' @examples
.idx.rx1day <- function(rst){
  r <-terra::app(rst, fun = max, message = FALSE)
}

###################################
#' RX5day
#'
#' @param rst 'SpatRaster' object with daily precipitation values.
#'
#' @return Maximum 5-day precipitation
#'
#' @examples
.idx.rx5day <- function(rst){
  
  ra <- function(x){
    v <- RcppRoll::roll_sum(x, 5)
    v <- max(v)
    return(v)
  }
  
  r <-terra::app(rst, fun = ra)
}


###################################
#' SDII
#'
#' @param rst 'SpatRaster' object with daily precipitation values.
#'
#' @return Simple precipitation intensity index
#'
#' @examples
.idx.sdii <- function(rst){
  
  rcl <- matrix(c(0, 0.99999, 0, 1, 10000, 1), byrow=TRUE, ncol = 3)
  v   <- terra::classify(rst, rcl)
  v   <- sum(v)
  s   <- sum(rst)
  
  rs  <- s / v
  
  return(rs)
  
}


###################################
#' CDD
#'
#' @param rst 'SpatRaster' object with daily precipitation values.
#'
#' @return Maximum length of dry spell: maximum number of consecutive days with RR < 1mm
#'
#' @examples
.idx.cdd <- function(rst){
  
  cdd <- function(x){
    
    x[x < 1]  <- -1
    x[x >= 1] <- 0
    x <-abs(x)
    
    rl <- rle(x)
    s  <- unlist(mapply(function(a,b) b * seq_len(a), rl$lengths, rl$values))
    return(max(s))
    
  }
  
  v <- terra::app(rst, cdd)
  
}

###################################
#' CWD
#'
#' @param rst 'SpatRaster' object with daily precipitation values.
#'
#' @return Maximum length of wet spell: maximum number of consecutive days with RR >= 1mm
#'
#' @examples
.idx.cwd <- function(rst){
  
  cwd <- function(x){
    
    x[x < 1]  <- 0
    x[x >= 1] <- 1
    x <-abs(x)
    
    rl <- rle(x)
    s  <- unlist(mapply(function(a,b) b * seq_len(a), rl$lengths, rl$values))
    return(max(s))
    
  }
  
  v <- terra::app(rst, cwd)
  
}

###################################
#' Rnnmm (threshold should be applied: for R10mm thres = 10; for R20mm thres = 20;
#' for Rnnmm thres = nn)
#'
#' @param rst 'SpatRaster' object with daily precipitation values.
#'
#' @return Annual count of days when PRCP >= a designed threshold
#'
#' @examples
.idx.rnnmm <- function(rst, thres){
  
  rcl <- matrix(c(0, thres - 0.0000001, 0, thres, 10000, 1), byrow=TRUE, ncol = 3)
  v   <- terra::classify(rst, rcl)
  v   <- sum(v)
  
  return(v)
}


###################################
#' Rnnp (threshold should be applied: for R95p thres = 0.95; for R99p thres = 0.99)
#'
#' @param rst 'SpatRaster' object with daily precipitation values.
#'
#' @return Annual total PRCP when RR > a designed percentile
#'
#' @examples
.idx.rnnp <- function(rst, thres){
  
  threspp <-function(x, thres){
    
    x <- x[which(x >= 1)]
    s <- quantile(x, thres)
    s <- as.numeric(s)
    s <- sum(x[which(x >= s)])
    return(s)
  }
  
  v <- terra::app(rst, threspp, thres = thres)
  
  return(v)
}

###################################
#' Rnnp (threshold should be applied: for R95pTOT thres = 0.95; for R99pTOT thres = 0.99)
#'
#' @param rst 'SpatRaster' object with daily precipitation values.
#'
#' @return Contribution to total precipitation from very wet days (0.95) and extremely wet days (0.99)
#'
#' @examples
.idx.rnnptot <- function(rst, thres){
  
  threspp <-function(x, thres){
    
    x <- x[which(x >= 1)]
    s <- quantile(x, thres)
    s <- as.numeric(s)
    s <- sum(x[which(x >= s)])
    return(s)
  }
  
  v <- terra::app(rst, threspp, thres = thres)
  s <- sum(rst)
  
  m <- (100 * v) / s
  
  return(m)
}


###################################
#' PRCPTOT 
#'
#' @param rst 'SpatRaster' object with daily precipitation values.
#'
#' @return Annual total precipitation on wet days
#'
#' @examples
.idx.prcptot <- function(rst){
  
  rst <- sum(rst)
  
}





