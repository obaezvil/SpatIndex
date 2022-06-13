################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: util functions for meteorological drought indices
################################################################################
##
## Creation date: 2022-06-12
##
################################################################################
################################################################################


#' Utils function to calculate the SPI
#'
#' @param x Numerical vector.
#' @param dates Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_spi' function.
#' @param scale Integer value that represents the time scale at which the SPI will be computed.
#' @param ref_start optional value that represents the starting point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.
#' @param ref_end Optional value that represents the ending point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.
#' @param distribution Optional value indicating the name of the distribution function to be used for computing the SPI 
#'  (one of 'log-Logistic', 'Gamma' and 'PearsonIII'). Defaults to 'Gamma' for SPI.
#' @param fit Optional value indicating the name of the method used for computing the distribution function parameters 
#'  (one of 'ub-pwm', 'pp-pwm' and 'max-lik'). Defaults to 'ub-pwm'.
#' @param ... Additional variables that can be used for the 'spi' function of the SPEI package.

#'
#' @return
#'
#' @examples
.spi <- function(x,
                 dates,
                 scale,
                 ref_start = NULL,
                 ref_end = NULL,
                 distribution = "Gamma",
                 fit = "ub-pwm",
                 ...){
  
  # Format ref_start if it is not NULL
  if(!is.null(ref_start)){
    
    yr        <- as.numeric(substr(ref_start, 1, 4))
    mo        <- as.numeric(substr(ref_start, 6, 7))
    ref_start <- c(yr, mo)
    
  }
  
  # Format ref_end if it is not NULL
  if(!is.null(ref_start)){
    
    yr      <- as.numeric(substr(ref_end, 1, 4))
    mo      <- as.numeric(substr(ref_end, 6, 7))
    ref_end <- c(yr, mo)
    
  }
  
  # Format first date
  yr      <- as.numeric(substr(dates[1], 1, 4))
  mo      <- as.numeric(substr(dates[1], 6, 7))
  dt      <- c(yr, mo)
  
  # Format y to a ts object
  x <- ts(x, start=dt, frequency = 12)
  
  # Apply the SPI
  spi <- SPEI::spi(x, scale = scale, distribution = distribution, fit = fit,
                   ref.start = ref_start, ref.end = ref_end, ...)
  
  spi <- as.numeric(spi$fitted)
  
  return(spi)
  
}




#' Utils function to calculate the SPI
#'
#' @param x Numerical vector.
#' @param dates Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_spei' function.
#' @param scale Integer value that represents the time scale at which the SPI will be computed.
#' @param ref_start optional value that represents the starting point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.
#' @param ref_end Optional value that represents the ending point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.
#' @param distribution Optional value indicating the name of the distribution function to be used for computing the SPI 
#'  (one of 'log-Logistic', 'Gamma' and 'PearsonIII'). Defaults to 'log-Logistic' for SPEI.
#' @param fit Optional value indicating the name of the method used for computing the distribution function parameters 
#'  (one of 'ub-pwm', 'pp-pwm' and 'max-lik'). Defaults to 'ub-pwm'.
#' @param ... Additional variables that can be used for the 'spi' function of the SPEI package.

#'
#' @return
#'
#' @examples
.spei <- function(x,
                 dates,
                 scale,
                 ref_start = NULL,
                 ref_end = NULL,
                 distribution = "log-Logistic",
                 fit = "ub-pwm",
                 ...){
  
  # Format ref_start if it is not NULL
  if(!is.null(ref_start)){
    
    yr        <- as.numeric(substr(ref_start, 1, 4))
    mo        <- as.numeric(substr(ref_start, 6, 7))
    ref_start <- c(yr, mo)
    
  }
  
  # Format ref_end if it is not NULL
  if(!is.null(ref_start)){
    
    yr      <- as.numeric(substr(ref_end, 1, 4))
    mo      <- as.numeric(substr(ref_end, 6, 7))
    ref_end <- c(yr, mo)
    
  }
  
  # Format first date
  yr      <- as.numeric(substr(dates[1], 1, 4))
  mo      <- as.numeric(substr(dates[1], 6, 7))
  dt      <- c(yr, mo)
  
  # Format y to a ts object
  x <- ts(x, start=dt, frequency = 12)
  
  # Apply the SPI
  spei <- SPEI::spei(x, scale = scale, distribution = distribution, fit = fit,
                   ref.start = ref_start, ref.end = ref_end, ...)
  
  spei <- as.numeric(spei$fitted)
  
  return(spei)
  
}


