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
## Creation date: 2022-06-12
##
################################################################################
################################################################################

spatial_spi <- function(P_data,
                        scale,
                        ref_start = NULL,
                        ref_end = NULL,
                        distribution = "Gamma",
                        fit = "ub-pwm",
                        ...){
  
  # Check P_data
  if(class(P_data) != "SpatRaster")
    stop("The object 'P_data' must be a SpatRaster")
  
  # Check scale
  if(class(scale) != "integer")
    stop("The object 'scale' must be a integer rhat represents the time scale at which the SPI will be computed")
  
  # Check ref_start object
  if(class(ref_start) != "character" & !is.null(ref_start))
    stop("If object 'ref_start' is not set to NULL, it must be a character object that indicates the starting point of the reference period used for computing the SPI. The format should be '%Y-%m'")
 
  # Check ref_end object
  if(class(ref_end) != "character" & !is.null(ref_end))
    stop("If object 'ref_end' is not set to NULL, it must be a character object that indicates the ending point of the reference period used for computing the SPI. The format should be '%Y-%m'")
  
}





