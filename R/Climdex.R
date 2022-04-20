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
## Date: 2020-04-20
##
################################################################################
################################################################################

#' Maximum 1-day precipitation
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
Rx1day <- function(rst.path, 
                   vct = NULL, 
                   temporal.scale = c("total", "annual"),
                   start.date,
                   end.date,
                   date.fmt = "%Y-%m-%d",
                   pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.rx1day <- function(rst){
    r <-terra::app(rst, fun = max)
  }
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.rx1day(rst = rst)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.rx1day(rst = rst)
    
  }
  
  return(results)
  
}


#' Maximum consecutive 5-day precipitation
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
Rx5day <- function(rst.path, 
                   vct = NULL, 
                   temporal.scale = c("total", "annual"),
                   start.date,
                   end.date,
                   date.fmt = "%Y-%m-%d",
                   pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.rx5day <- function(rst){
    
    ra <- function(x){
      v <- RcppRoll::roll_sum(x, 5)
      v <- max(v)
      return(v)
    }
    
    r <-terra::app(rst, fun = ra)
  }
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.rx5day(rst = rst)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.rx5day(rst = rst)
    
  }
  
  return(results)
  
}

#' Maximum consecutive 5-day precipitation
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
sdii <- function(rst.path, 
                 vct = NULL, 
                 temporal.scale = c("total", "annual"),
                 start.date,
                 end.date,
                 date.fmt = "%Y-%m-%d",
                 pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.sdii <- function(rst){
    
    rcl <- matrix(c(0, 0.9999, 0, 1, 10000, 1), byrow=TRUE, ncol = 3)
    v   <- terra::classify(rst, rcl)
    v   <- sum(v)
    s   <- sum(rst)
    
    rs  <- s / v
    
    return(rs)
    
  }
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.sdii(rst = rst)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.sdii(rst = rst)
    
  }
  
  return(results)
  
}

#' Maximum length of dry spell: maximum number of consecutive days with RR < 1mm
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
CDD <- function(rst.path, 
                vct = NULL, 
                temporal.scale = c("total", "annual"),
                start.date,
                end.date,
                date.fmt = "%Y-%m-%d",
                pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.cdd <- function(rst){
    
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
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.cdd(rst = rst)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.cdd(rst = rst)
    
  }
  
  return(results)
  
}

#' Maximum length of dry spell: maximum number of consecutive days with RR < 1mm
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
CWD <- function(rst.path, 
                vct = NULL, 
                temporal.scale = c("total", "annual"),
                start.date,
                end.date,
                date.fmt = "%Y-%m-%d",
                pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.cwd <- function(rst){
    
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
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.cwd(rst = rst)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.cwd(rst = rst)
    
  }
  
  return(results)
  
}

#' Annual count of days when PRCP ≥ 10mm
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
R10mm <- function(rst.path, 
                  vct = NULL, 
                  temporal.scale = c("total", "annual"),
                  start.date,
                  end.date,
                  date.fmt = "%Y-%m-%d",
                  pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.rnnmm <- function(rst, thres){
    
    rcl <- matrix(c(0, thres - 0.0000001, 0, thres, 10000, 1), byrow=TRUE, ncol = 3)
    v   <- terra::classify(rst, rcl)
    v   <- sum(v)
    
    return(v)
  }
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.rnnmm(rst = rst, thres = 10)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.rnnmm(rst = rst, thres = 10)
    
  }
  
  return(results)
  
}

#' Annual count of days when PRCP ≥ 20mm
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
R20mm <- function(rst.path, 
                  vct = NULL, 
                  temporal.scale = c("total", "annual"),
                  start.date,
                  end.date,
                  date.fmt = "%Y-%m-%d",
                  pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.rnnmm <- function(rst, thres){
    
    rcl <- matrix(c(0, thres - 0.0000001, 0, thres, 10000, 1), byrow=TRUE, ncol = 3)
    v   <- terra::classify(rst, rcl)
    v   <- sum(v)
    
    return(v)
  }
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.rnnmm(rst = rst, thres = 20)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.rnnmm(rst = rst, thres = 20)
    
  }
  
  return(results)
  
}

#' Annual count of days when PRCP ≥ nn mm, where nn is a user-defined threshold
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#'  @param thres user-defined threshold
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
Rnnmm <- function(rst.path, 
                  vct = NULL, 
                  thres,
                  temporal.scale = c("total", "annual"),
                  start.date,
                  end.date,
                  date.fmt = "%Y-%m-%d",
                  pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.rnnmm <- function(rst, thres){
    
    rcl <- matrix(c(0, thres - 0.0000001, 0, thres, 10000, 1), byrow=TRUE, ncol = 3)
    v   <- terra::classify(rst, rcl)
    v   <- sum(v)
    
    return(v)
  }
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.rnnmm(rst = rst, thres = thres)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.rnnmm(rst = rst, thres = thres)
    
  }
  
  return(results)
  
}


#' Annual total PRCP when RR > 95th percentile
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
R95p <- function(rst.path, 
                  vct = NULL, 
                  temporal.scale = c("total", "annual"),
                  start.date,
                  end.date,
                  date.fmt = "%Y-%m-%d",
                  pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.rnnp <- function(rst, thres){
    
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
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.rnnp(rst = rst, thres = 0.95)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.rnnp(rst = rst, thres = 0.95)
    
  }
  
  return(results)
  
}

#' Annual total PRCP when RR > 99th percentile
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
R99p <- function(rst.path, 
                 vct = NULL, 
                 temporal.scale = c("total", "annual"),
                 start.date,
                 end.date,
                 date.fmt = "%Y-%m-%d",
                 pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.rnnp <- function(rst, thres){
    
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
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.rnnp(rst = rst, thres = 0.99)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.rnnp(rst = rst, thres = 0.99)
    
  }
  
  return(results)
  
}


#' Contribution to total precipitation from very wet days (95p)
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#'  @param thres user-defined threshold
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
R95pTOT <- function(rst.path, 
                    vct = NULL, 
                    temporal.scale = c("total", "annual"),
                    start.date,
                    end.date,
                    date.fmt = "%Y-%m-%d",
                    pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.rnnptot <- function(rst, thres){
    
    threspp <-function(x, thres){
      
      x <- x[which(x >= 1)]
      s <- quantile(x, thres)
      s <- as.numeric(s)
      s <- sum(x[which(x >= s)])
      return(s)
    }
    
    v <- terra::app(rst, threspp, thres = thres)
    s <- sum(rast)
    
    m <- (100 * v) / s
    
    return(m)
  }
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.rnnptot(rst = rst, thres = 0.95)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.rnnptot(rst = rst, thres = 0.95)
    
  }
  
  return(results)
  
}


#' Contribution to total precipitation from extremely wet days (95p)
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#'  @param thres user-defined threshold
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
R99pTOT <- function(rst.path, 
                    vct = NULL, 
                    temporal.scale = c("total", "annual"),
                    start.date,
                    end.date,
                    date.fmt = "%Y-%m-%d",
                    pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # defining function to be applied
  idx.rnnptot <- function(rst, thres){
    
    threspp <-function(x, thres){
      
      x <- x[which(x >= 1)]
      s <- quantile(x, thres)
      s <- as.numeric(s)
      s <- sum(x[which(x >= s)])
      return(s)
    }
    
    v <- terra::app(rst, threspp, thres = thres)
    s <- sum(rast)
    
    m <- (100 * v) / s
    
    return(m)
  }
  
  # analysing temporal scale to be used
  if(temporal.scale == "annual"){
    
    # Subsetting the years
    years   <- substr(dts, 1, 4)
    u.years <- unique(years)
    
    # creating an empty list to store the results
    results <- list()
    
    for(i in 1:length(u.years)){
      
      # Subsetting files
      pos <- which(years %in% u.years[i])
      rst <- r[pos]
      rst <- terra::rast(rst)
      
      # cropping to extent if given
      if(!is.null(vct)){
        rst <- terra::crop(rst, vct, snap = "out")
      }
      
      r.ind        <- idx.rnnptot(rst = rst, thres = 0.99)
      names(r.ind) <- u.years[i]
      results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
    
  } else {
    
    rst <- terra::rast(r)
    
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
    
    results <- idx.rnnptot(rst = rst, thres = 0.99)
    
  }
  
  return(results)
  
}


#' Contribution to total precipitation from extremely wet days (95p)
#'
#' @param rst.path File path to daily raster files for the period of analysis. 
#'  These files should include the date in any format.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#'  @param thres user-defined threshold
#' @param temporal.scale either 'total' to use all the period of record, or
#'  'annual' to compute the index annually (From Jan to Dec)
#' @param start.date Position where the dates of the raster files start
#' @param end.date Position where the dates of the raster files end
#' @param date.fmt Format of the dates included in the file names (default = 
#'  "%Y-%m-%d")
#' @param pattern Set to NULL. Is there a specific pattern to list the 
#'  raster files?
#'
#' @return
#' @export
#'
#' @examples
PRCPTOT <- function(rst.path, 
                    vct = NULL, 
                    temporal.scale = c("total", "annual"),
                    start.date,
                    end.date,
                    date.fmt = "%Y-%m-%d",
                    pattern = NULL){
  
  # checking temporal.scale input
  if(!temporal.scale %in% c("total", "annual"))
    stop("Please provide a valid 'temporal.scale' (i.e., 'total' or 'annual'.")
  
  # listing raster files 
  if(is.null(pattern)){
    r <- list.files(rst.path, full.names = TRUE)
  } else {
    r <- list.files(rst.path, full.names = TRUE, pattern = pattern)
  }
  
  # checking availability of files
  if(length(r) < 1)
    stop("There are no files inside the 'rst.path' directory!")
  
  # getting the dates from the files
  dts <- substr(basename(r), start.date, end.date)
  dts <- as.Date(dts, format = date.fmt)
  
  # analysing temporal scale to be used (this applies for any scale)
  # Subsetting the years
  years   <- substr(dts, 1, 4)
  u.years <- unique(years)
    
  # creating an empty list to store the results
  results <- list()
    
  for(i in 1:length(u.years)){
      
    # Subsetting files
    pos <- which(years %in% u.years[i])
    rst <- r[pos]
    rst <- terra::rast(rst)
      
    # cropping to extent if given
    if(!is.null(vct)){
      rst <- terra::crop(rst, vct, snap = "out")
    }
      
    r.ind        <- sum(rst)
    names(r.ind) <- u.years[i]
    results[[i]] <- r.ind
    } # end for
    
    results <- terra::rast(results)
  
}



