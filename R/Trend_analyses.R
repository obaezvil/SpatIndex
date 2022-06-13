################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: function to perform trend analysis based on the Man-Kendall test
################################################################################
##
## Creation date: 2022-05-13
##
################################################################################
################################################################################

#' Non-parametric Man-Kendall trend test
#'
#' @param rst 'SpatRaster' object.
#' @param vct Vector file of the study area (Optional). It will be used to crop 
#'  the spatial extent of the raster files if required. This parameter is set 
#'  to 'NULL' 
#' @param conf_level Level of significance. Set to 0.95 by default.
#' @param p_thres threshold of p-value. Set to 0.05 by default.
#'
#' @return returns gridded layers(i.e., grid-cells where p-value < p_thres, 
#'  p-value, Kendall_Score-S, Variance_of_Kendall_Score-varS, Tau values, 
#'  "Sen's_Slope", "Total_change") and a vector layer (i.e., grid-cells where p-value < p_thres).
#' @export
#'
#' @examples
mk_spatial <- function(rst, 
                       vct = NULL,
                       conf_level = 0.95, 
                       p_thres = 0.05){
  
  # Controlling the class of the object
  !if(class(rst) %in% c("SpatRaster", "RasterBrick", "RasterStack"))
    stop("The object 'rst' must have class 'SpatRaster', 'RasterBrick', or 'RasterStack'!")
  
  # Checking the 'vct' object
  if(!is.null(vct) & class(vct) != "SpatVector")
    stop("The object 'vct' is not a 'SpatVector'!")
  
  # Converting raster object to SpatRaster
  if(class(rst) %in% c("RasterBrick", "RasterStack"))
     rst <- terra::rast(rst)
  
  # cropping to extent if given
  if(!is.null(vct)){
    rst <- terra::crop(rst, vct, snap = "out")
    rst <- terra::mask(rst, vct)
  }
  
  # Applying Man-Kendall trend test
  mk <- terra::app(rst, .MKpval)
  
  # Setting names of reult layers
  names(mk) <- c("p-value", "Kendall_Score-S", "Variance_of_Kendall_Score-varS", "Tau")
  
  # Applying the Sen's Slope
  ss <- terra::app(rst, .Senslope, conf.level = conf_level)
  
  # Setting names of reult layers
  names(ss) <- c("Sen's_Slope", "Total_change")
  
  # Creating a SpatVect with the spatRast geometry
  points <- terra::as.points(mk[[1]])
  
  # Excluding points outside of the 'p_thres' argument
  points <- points[points$`p-value` < p_thres,]
  names(points) <- "Significant_points"
  
  # Creating boolean layer of significant areas
  significant_cells        <- mk$`p-value` < p_thres
  names(significant_cells) <- "Significant_cells"
  
  # Aggregating results
  results        <- significant_cells
  terra::add(results) <- c(mk, ss)
  
  if(length(points) < 1){
    warning("The 'Significant_points' layer was not generated because any p-value < ", p_thres, "!")
    results <- list(gridded = results,
                    vector = NULL)
  } else {
    results <- list(gridded = results,
                    vector = points)
  }

  return(results)
}

#' Utils non-parametric man-Kendall trend test
#'
#' @param x numerical vector.
#'
#' @return returns a numerical vector representing the p-value, Kendall_Score-S,
#'   Variance_of_Kendall_Score-varS, and Tau.
#'
#' @examples
.MKpval <- function(x){
  
  # Condition to deal with NA values
  if(!any(is.na(x))){
    
    # Applying Man-Kendall trend test
    mk <- trend::mk.test(x)
    
    # Extracting results
    pval <- mk$p.value
    S    <- mk$estimates[1]
    varS <- mk$estimates[2]
    tau  <- mk$estimates[3]
    
  } else {
    pval <- S <- varS <- tau <- NA
  }
  
  return(as.numeric(c(pval, S, varS, tau)))
  
}

#' Utils function that calculates the Sen's slope
#'
#' @param x numerical vector.
#' @param conf_level Level of significance. Set to 0.95 by default.
#'
#' @return
#' @export
#'
#' @examples
.Senslope <- function(x, conf_level = 0.95){
  
  # Condition to deal with NA values
  if(!any(is.na(x))){
    
    # Calculating the Sen's Slope
    ss     <- trend::sens.slope(x, conf.level = conf_level)
    ss     <- as.numeric(ss$estimates)
    change <- ss * length(x)
    
  } else {
    ss <- change <- NA
  }
  
  return(as.numeric(c(ss, change)))
  
}
