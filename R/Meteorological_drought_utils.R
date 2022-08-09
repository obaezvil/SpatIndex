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


#' Utils function to calculate the SPI (SPEI package)
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
.spi.spei <- function(x,
                      dates,
                      scale,
                      ref_start,
                      ref_end,
                      distribution,
                      fit,
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
  
  # Trasforming NaNs to NAs
  spi[is.nan(spi)] <- NA
  
  return(spi)
  
}


#' Utils function to calculate the SPI (SCI package)
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
.spi.sci <- function(x,
                      dates,
                      scale,
                      ref_start,
                      ref_end,
                      distribution,
                      fit,
                      ...){
  
  # Setting by default 'pos_ini' and 'pos_fin' as the first and last layer
  pos_ini <- 1
  pos_fin <- length(x)
  
  # Format ref_start if it is not NULL
  if(!is.null(ref_start)){
    
    pos_ini <- grep(ref_start, dates)
    
  }
  
  # Format ref_end if it is not NULL
  if(!is.null(ref_start)){
    
    pos_fin <- grep(ref_end, dates)
    
  }
  
  # Checking distribution
  if(!distribution %in% c("Gamma", "log-Logistic", "PearsonIII"))
    stop("The accepted distributions for the SCI package are 'Gamma', 'log-Logistic', and 'PearsonIII'")
    
  # Adjusting distributions for SCI package
  if(distribution == "Gamma")
    distribution <- tolower(distribution)
  
  if(distribution == "log-Logistic")
    distribution <- 'genlog'
  
  if(distribution == "PearsonIII")
    distribution <- 'pe3'
  
  # Obtaining parameters for the reference period
  spi_ref_params <- SCI::fitSCI(x[pos_ini:pos_fin], first.mon = 1, 
                                distr = distribution, 
                                time.scale = scale, p0 = TRUE)
  
  # apply the the parameters of the reference period to all data
  spi <- SCI::transformSCI(x, first.mon = 1, obj = spi_ref_params)
  
  # Trasforming NaNs to NAs
  spi[is.nan(spi)] <- NA
  
  return(spi)
  
}




#' Utils function to calculate the SPEI (SPEI package)
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
.spei.spei <- function(x,
                       dates,
                       scale,
                       ref_start,
                       ref_end,
                       distribution,
                       fit,
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
  
  # Trasforming NaNs to NAs
  spei[is.nan(spei)] <- NA
  
  return(spei)
  
}

#' Utils function to calculate the SPEI (SCI package)
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
.spei.sci <- function(x,
                       dates,
                       scale,
                       ref_start,
                       ref_end,
                       distribution,
                       fit,
                       ...){
  
  # Setting by default 'pos_ini' and 'pos_fin' as the first and last layer
  pos_ini <- 1
  pos_fin <- length(x)
  
  # Format ref_start if it is not NULL
  if(!is.null(ref_start)){
    
    pos_ini <- grep(ref_start, dates)
    
  }
  
  # Format ref_end if it is not NULL
  if(!is.null(ref_start)){
    
    pos_fin <- grep(ref_end, dates)
    
  }
  
  # Checking distribution
  if(!distribution %in% c("Gamma", "log-Logistic", "PearsonIII"))
    stop("The accepted distributions for the SCI package are 'Gamma', 'log-Logistic', and 'PearsonIII'")
  
  # Adjusting distributions for SCI package
  if(distribution == "Gamma")
    distribution <- tolower(distribution)
  
  if(distribution == "log-Logistic")
    distribution <- 'genlog'
  
  if(distribution == "PearsonIII")
    distribution <- 'pe3'
  
  # Obtaining parameters for the reference period
  spei_ref_params <- SCI::fitSCI(x[pos_ini:pos_fin], first.mon = 1, 
                                distr = distribution, 
                                time.scale = scale, p0 = TRUE)
  
  # apply the the parameters of the reference period to all data
  spei <- SCI::transformSCI(x, first.mon = 1, obj = spei_ref_params)
  
  # Trasforming NaNs to NAs
  spei[is.nan(spei)] <- NA
  
  return(spei)
  
}


