################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: calculate combined drought indices
################################################################################
##
## Creation date: 2022-11-27
##
################################################################################
################################################################################

#' Utils function to calculate the boolean layers to calculate the zSPI
#'
#' @param SPI_data 'SpatRaster' object that contains the SPI data that will be used to calculate the zSPI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param threshold Threshold to construct the Boolean layers.
#'
#'
#' @examples
.spi_bool <- function(SPI_data, threshold){
  
  # Creating the matrix for reclassifying the data
  rcl <- matrix(c(threshold,   Inf, 0,
                  -Inf, threshold, 1), ncol = 3, byrow = TRUE)
  
  # Classifying the SPI data according to the selected threshold
  spi_bol <- terra::classify(SPI_data, rcl, right = FALSE)
  
  return(spi_bol)
  
}

#' Utils function to calculate the boolean SPI (zSPI) to be used in the CDI
#'
#' @param SPI1_data 'SpatRaster' object that contains the SPI-1 data that will be used to calculate the zSPI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param SPI3_data 'SpatRaster' object that contains the SPI-3 data that will be used to calculate the zSPI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param threshold_spi1 Threshold to construct the Boolean layer for SPI-1. Set to -2 by default.
#' @param threshold_spi3 Threshold to construct the Boolean layer for SPI-3. Set to -1 by default.
#'
#' @return Returns the zSPI
#'
#' @examples
.zspi <- function(SPI1_data, SPI3_data, threshold_spi1, threshold_spi3){   #### ADDD threshold options in the main code XXXXXXXX
  
  # Calculating zSPI1
  spi1_bool <- .spi_bool(SPI1_data, threshold_spi1)
  
  # Calculating zSPI3
  spi3_bool <- .spi_bool(SPI3_data, threshold_spi3)
  
  # Creating unitary 'SpatRaster' as should fulfill the following conditions:
  #  1 when SPI-1 < -2 or SPI-3 < -1; otherwise 0
  zspi           <- spi1_bool + spi3_bool
  zspi[zspi > 0] <- 1
  
  return(zspi)
  
}



#' Utils function to calculate the values to be used in the CDI
#'
#' @param zSPI 'SpatRaster' object that contains the zSPI data that will be used to calculate the CDI 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param zSM 'SpatRaster' object that contains the zSM data that will be used to calculate the CDI 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param zfAPAR 'SpatRaster' object that contains the zSfAPAR data that will be used to calculate the CDI 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#'
#' @return Return the 3-digit values to compute the CDI.
#'
#' @examples
.cdi_vals <- function(zSPI, zSM, zfAPAR){
  
  # Create reclassification matrices according to the different thresholds of the variables (Cammalleri et al., 2021):
  #   For zSPI the values could be 1 or 0 -> values given: 1 and 2, respectively
  #   For zSM the ranges could be [0, Inf); [-1, 0); and (-Inf, -1) -> values given: 3, 2, and 1, respectively
  #   For zfAPAR the ranges could be [0, Inf); [-1, 0); and (-Inf, -1) -> values given: 3, 2, and 1, respectively
  rcl <- matrix(c(   0, Inf,   3,
                    -1,   0,   2,
                  -Inf,  -1,   1), ncol = 3, byrow = TRUE)
  
  rcl_spi <- matrix(c(1,  1,
                      0,  2), ncol = 2, byrow = TRUE)
  
  # Classifying the variables:
  #   If zSPI == 1; 1, otherwise 2 (zSPI ==0)
  #   If zSM or zfAPAR == [0, Inf) -> 3; [-1, 0) -> 2; and (-Inf, -1) -> 1
  zSPI   <- terra::classify(zSPI, rcl_spi)
  zSM    <- terra::classify(zSM, rcl, right = FALSE)
  zfAPAR <- terra::classify(zfAPAR, rcl, right = FALSE)
  
  # A layer with all results will be generated. 
  #   The value for the zSPI will be stored in the hundredths (either 100 or 200)
  #   The value for the zSM will be stored in the tenths (either 10, 20, or 30)
  #   The value for the zfAPAR will be stored in the unity (either 100, 200, or 300)
  #   Therefore, the grid-cells will be composed by one 3-digit value: hundredths (zSPI); tenths (zSM); and unit (zfAPAR)
  #                1st     2nd      3rd
  #           |--------|--------|---------|
  #           | zSPI   |  zSM   |  zfAPAR |
  #           | [1, 2] | [1, 3] | [1, 4]  |   
  #           |--------|--------|---------|
  zSPI <- zSPI * 100
  zSM  <- zSM  * 10
  
  cdi_vals <- zSPI + zSM + zfAPAR
  
  return(cdi_vals)
  
}



#' Calculate the CDI index
#'
#' @param x vector object with the three-digit CDI values.
#'
#' @return Return a vector with CDI values.
#'
#' @examples
.get_cdi <- function(x){
  
    # Imputing missing gaps where NAs are present, the maximum consecutive gap acceptable is of 4 elements   
    n <- which(is.na(x))
    n <- length(n)
    if(n != length(x)){
      x_imp <- as.numeric(imputeTS::na_locf(x, maxgap = 6))
    } else {
      x_imp <- x
    }
      
    
    # Setting to NA the result if there is still a NA after the TS imputation
    nas <- which(is.na(x_imp))
    nas <- length(nas)
    if(nas == 0){
    
    # Starting the CDI computation considering cdi(t-1) as zero
    cdi_tm1 <- 0
    cdi     <- NA
    
    for(i in 1:length(x_imp)){
      
      x_ti     <- x_imp[i]
      cdi_prev <- cdi_tm1[i]
      
      # The following values and ranges will take values from 1 to 3 to make the computation easier
      # zSPI: 2 = 0; 1 = 1
      # zSM: [0, Inf) = 3; [-1, 0) = 2; and (-Inf, -1) = 1
      # zfAPAR: [0, Inf) = 3; [-1, 0) = 2; and (-Inf, -1) = 1
      #
      #                                       zSPI
      #              |-------------------------|-------------------------|
      #              2                                                   1  
      #             zSM                                                 zSM
      #     |--------|--------|                                 |--------|--------|
      #     3        2        1                                 3        2        1
      #   zfAPAR   zfAPAR   zfAPAR                            zfAPAR   zfAPAR   zfAPAR
      #  |--|--|  |--|--|  |--|--|                           |--|--|  |--|--|  |--|--| 
      #  3  2  1  3  2  1  3  2  1                           3  2  1  3  2  1  3  2  1
      #  
      #  ----------------------------------------------------------------------------------------------------------------------------
      #  Number             233   232   231   223   222   221   213   212   211   133   132   131   123   122   121   113   112   111
      #  ----------------------------------------------------------------------------------------------------------------------------
      #  CDIt-1 = 0,4         0     0     0     0     0     0     0     0     0     1     1     3     1     1     3     2     2     3
      #  CDIt-1 = 1           4     4     3     4     4     3     2     2     3     1     1     3     1     1     3     2     2     3
      #  CDIt-1 = 2,5         5     4     3     5     4     3     2     2     3     5     1     3     5     1     3     2     2     3
      #  CDIt-1 = 3,6         6     4     3     6     4     3     6     2     3     6     1     3     6     1     3     6     2     3
      #  ----------------------------------------------------------------------------------------------------------------------------
      
      # Condition if 'cdi_prev' is zero or four (if TRUE, there are four values possible: 0, 1, 2, and 3)
      if(cdi_prev == 0 | cdi_prev == 4){
        
        if(x_ti %in% c(233, 232, 231, 223, 222, 221, 213, 212, 211)){
          cdi[i] <- 0
        } else if(x_ti %in% c(133, 132, 123, 122)){
          cdi[i] <- 1
        } else if(x_ti %in% c(113, 112, 111)){
          cdi[i] <- 2
        } else if(x_ti %in% c(131, 121)){
          cdi[i] <- 3
        } else {
          cdi[i] <- NA
        }
        
      } # end condition if 'cdi_prev' is zero or four
      
      # Condition if 'cdi_prev' is one (if TRUE, there are four values possible: 1, 2, 3, and 4)
      if(cdi_prev == 1){
        
        if(x_ti %in% c(133, 132, 123, 122)){
          cdi[i] <- 1
        } else if(x_ti %in% c(213, 212, 113, 112)){
          cdi[i] <- 2
        } else if(x_ti %in% c(231, 221, 211, 131, 121, 111)){
          cdi[i] <- 3
        } else if(x_ti %in% c(233, 232, 223, 222)){
          cdi[i] <- 4
        } else {
          cdi[i] <- NA
        }
        
      } # end condition if 'cdi_prev' is zero
      
      # Condition if 'cdi_prev' is two or five (if TRUE, there are four values possible: 1, 2, 3, 4 and 5)
      if(cdi_prev == 2 | cdi_prev == 5){
        
        if(x_ti %in% c(132, 122)){
          cdi[i] <- 1
        } else if(x_ti %in% c(213, 212, 113, 112)){
          cdi[i] <- 2
        } else if(x_ti %in% c(231, 221, 211, 131, 121, 111)){
          cdi[i] <- 3
        } else if(x_ti %in% c(232, 222)){
          cdi[i] <- 4
        } else if(x_ti %in% c(233, 223, 133, 123)){
          cdi[i] <- 5
        } else {
          cdi[i] <- NA
        }
        
      } # end condition if 'cdi_prev' is two or five
      
      # Condition if 'cdi_prev' is three or six (if TRUE, there are four values possible: 1, 2, 3, 4 and 6)
      if(cdi_prev == 3 | cdi_prev == 6){
        
        if(x_ti %in% c(132, 122)){
          cdi[i] <- 1
        } else if(x_ti %in% c(212, 112)){
          cdi[i] <- 2
        } else if(x_ti %in% c(231, 221, 211, 131, 121, 111)){
          cdi[i] <- 3
        } else if(x_ti %in% c(232, 222)){
          cdi[i] <- 4
        } else if(x_ti %in% c(233, 223, 213, 133, 123, 113)){
          cdi[i] <- 6
        } else {
          cdi[i] <- NA
        }
        
      } # end condition if 'cdi_prev' is two or five
      
      # Setting the value to the 'cdi_tm1' for the next iteration
      cdi_tm1[i+1] <- cdi[i]
      
    } # enf for 'i'
    
  } else { # condition that states that the TS is full of NAs
    
    cdi <- rep(NA, length(x_imp))
    
  }

  return(cdi)
  
}







