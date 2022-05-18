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

mk_spatial <- function(rst, 
                       vct = NULL,
                       conf.level = 0.95, 
                       pval.thres = 0.05){
  
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
  ss <- terra::app(rst, .Senslope, conf.level = conf.level)
  
  # Setting names of reult layers
  names(ss) <- c("Sen's_Slope", "Total_change")
  
  # Creating a SpatVect with the spatRast geometry
  points <- terra::as.points(results[[1]])
  
  # Excluding points outside of the 'pval.thres' argument
  points <- points[points$`p-value` < pval.thres,]
  names(points) <- "Significant_points"
  
  # Creating boolean layer of significant areas
  significant_cells        <- mk$`p-value` < pval.thres
  names(significant_cells) <- "Significant_cells"
  
  # Aggregating results
  results        <- significant_cells
  terra::add(mk) <- c(mk, ss)
  
  if(length(points) < 1){
    warning("The 'Significant_points' layer was not generated because any p-value < ", pval.thres, "!")
    results <- list(gridded = results,
                    vector = NULL)
  } else {
    results <- list(gridded = results,
                    vector = points)
  }

  return(results)
}

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

.Senslope <- function(x, conf.level = 0.95){
  
  # Condition to deal with NA values
  if(!any(is.na(x))){
    
    # Calculating the Sen's Slope
    ss     <- trend::sens.slope(x, conf.level = conf.level)
    ss     <- as.numeric(ss$estimates)
    change <- ss * length(x)
    
  } else {
    ss <- change <- NA
  }
  
  return(as.numeric(c(ss, change)))
  
}
