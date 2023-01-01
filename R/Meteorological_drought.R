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

#' Wrapper function from the SCI and SPEI packages to calculate the spatially-distributed Standardized Precipitation Index (SPI)
#'
#' @param P_data 'SpatRaster' object that contains spatially-distributed monthly precipitation data that will be used to calculate the SPI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
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
#' @param na.rm Should the NA values be removed? Set to TRUE.
#' @param package Either 'SCI' or 'SPEI'. Should the SCI or SPEI package be used in the implementation?
#' @param ... Additional variables that can be used for the 'spi' function of the SPEI package.
#'
#' @return Spatially-distributed SPI values.
#' @export
#'
#' @examples
spatial_spi <- function(P_data,
                        scale,
                        ref_start = NULL,
                        ref_end = NULL,
                        distribution = "Gamma",
                        fit = "ub-pwm",
                        na.rm = TRUE,
                        package = "SCI",
                        ...){
  
  # Check P_data
  if(class(P_data) != "SpatRaster")
    stop("The object 'P_data' must be a SpatRaster")
  
  # Check scale
  if(!class(scale) %in% c("integer", "numeric"))
    stop("The object 'scale' must be a integer that represents the time scale at which the SPI will be computed")
  
  # Check ref_start object
  if(class(ref_start) != "character" & !is.null(ref_start))
    stop("If object 'ref_start' is not set to NULL, it must be a character object that indicates the starting point of the reference period used for computing the SPI. The format should be '%Y-%m'")
  
  # Check ref_end object
  if(class(ref_end) != "character" & !is.null(ref_end))
    stop("If object 'ref_end' is not set to NULL, it must be a character object that indicates the ending point of the reference period used for computing the SPI. The format should be '%Y-%m'")
  
  # Check ref_start and ref_end
  if(is.null(ref_start) & !is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  if(!is.null(ref_start) & is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  # Checking the 'package' object
  if(!package %in% c("SPEI", "SCI"))
    stop("The 'package' object must be either 'SPEI' or 'SCI'")
  
  # Extract dates from object
  dates <- terra::time(P_data)
  
  # Apply SPI
  if(package == "SPEI"){
    
    idx <- terra::app(P_data, .spi.spei, scale = scale, dates = dates, distribution = distribution, fit = fit,
                      ref_start = ref_start, ref_end = ref_end, na.rm = na.rm, ...)
    
  } else {
    
    idx <- terra::app(P_data, .spi.sci, scale = scale, dates = dates, distribution = distribution, fit = fit,
                      ref_start = ref_start, ref_end = ref_end, na.rm = na.rm, ...)
    
  }
  
  
  ## set dates and return
  names(idx)        <- paste0(substr(dates, 1, 7)) 
  terra::time(idx)  <- dates
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  idx[is.infinite(idx)] <- NA
  
  return(idx)
  
}

#' Wrapper function from the SCI and SPEI packages to calculate the spatially-distributed Standardized Precipitation-Evapotranspiration Index (SPEI)
#'
#' @param P_data 'SpatRaster' object that contains spatially-distributed monthly precipitation data that will be used to calculate the SPEI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param PE_data 'SpatRaster' object that contains spatially-distributed monthly potential evaporation (also referred to as potential evapotranspiation) 
#'  data that will be used to calculate the SPEI. This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. 
#'  They can be set with the function time of the terra package.
#' @param scale Integer value that represents the time scale at which the SPI will be computed.
#' @param ref_start Optional value that represents the starting point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.
#' @param ref_end Optional value that represents the ending point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.
#' @param distribution Optional value indicating the name of the distribution function to be used for computing the SPEI 
#'  (one of 'log-Logistic', 'Gamma' and 'PearsonIII'). Defaults to 'log-Logistic' for SPEI.
#' @param fit Optional value indicating the name of the method used for computing the distribution function parameters 
#'  (one of 'ub-pwm', 'pp-pwm' and 'max-lik'). Defaults to 'ub-pwm'.
#' @param na.rm Should the NA values be removed? Set to TRUE.
#' @param package Either 'SCI' or 'SPEI'. Should the SCI or SPEI package be used in the implementation?
#' @param ... Additional variables that can be used for the 'spi' function of the SPEI package.
#'
#' @return Spatially-distributed SPEI values.
#' @export
#'
#' @examples
spatial_spei <- function(P_data,
                         PE_data,
                        scale,
                        ref_start = NULL,
                        ref_end = NULL,
                        distribution = "log-Logistic",
                        fit = "ub-pwm",
                        na.rm = TRUE,
                        package = "SCI",
                        ...){
  
  # Check P_data
  if(class(P_data) != "SpatRaster")
    stop("The object 'P_data' must be a SpatRaster")
  
  # Check PE_data
  if(class(PE_data) != "SpatRaster")
    stop("The object 'PE_data' must be a SpatRaster")
  
  # Check scale
  if(!class(scale) %in% c("integer", "numeric"))
    stop("The object 'scale' must be a integer that represents the time scale at which the SPEI will be computed")
  
  # Check ref_start object
  if(class(ref_start) != "character" & !is.null(ref_start))
    stop("If object 'ref_start' is not set to NULL, it must be a character object that indicates the starting point of the reference period used for computing the SPEI. The format should be '%Y-%m'")
  
  # Check ref_end object
  if(class(ref_end) != "character" & !is.null(ref_end))
    stop("If object 'ref_end' is not set to NULL, it must be a character object that indicates the ending point of the reference period used for computing the SPEI. The format should be '%Y-%m'")
  
  # Check ref_start and ref_end
  if(is.null(ref_start) & !is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  if(!is.null(ref_start) & is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  # Extract dates from P_data object
  dates_p <- terra::time(P_data)
  
  # Extract dates from PE_data object
  dates_pe <- terra::time(PE_data)
  
  # Checking temporal dimension
  if(length(dates_p) != length(dates_pe))
    stop("The lengths of 'P_data' and 'PE_data' differ")
  
  # Checking that P and PE have the same period
  if(!all(dates_p %in% dates_pe))
    stop("'P_data' and 'PE_data' have different periods")
    
  # Check extent of products to resample
  if(!terra::compareGeom(P_data, PE_data)){
    # P higher resolution that PE
    if(terra::res(P_data)[1] < terra::res(PE_data)[1]){
      PE_data <- terra::resample(PE_data, P_data, method = "near")
    } else {
      P_data <- terra::resample(P_data, PE_data, method = "near")
    }
  }
  
  # Checking the 'package' object
  if(!package %in% c("SPEI", "SCI"))
    stop("The 'package' object must be either 'SPEI' or 'SCI'")
  
  # Compute P-PE
  diff <- P_data - PE_data
  
  # Apply SPEI
  if(package == "SPEI"){
    
    idx <- terra::app(diff, .spei.spei, scale = scale, dates = dates_p, distribution = distribution, fit = fit,
                      ref_start = ref_start, ref_end = ref_end, na.rm = na.rm, ...)
    
  } else {
    
    idx <- terra::app(diff, .spei.sci, scale = scale, dates = dates_p, distribution = distribution, fit = fit,
                      ref_start = ref_start, ref_end = ref_end, na.rm = na.rm, ...)
    
  }
  
  
  ## set dates and return
  names(idx)        <- paste0(substr(dates_p, 1, 7)) 
  terra::time(idx)  <- dates_p
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  idx[is.infinite(idx)] <- NA
  
  return(idx)
  
}


#' Function to calculate the Deciles method
#'
#' @param P_data 'SpatRaster' object that contains spatially-distributed monthly precipitation data that will be used to calculate the SPI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @return Spatially-distributed Decile values.
#' @export
#'
#' @examples
spatial_deciles <- function(P_data){
  
  # Check P_data
  if(class(P_data) != "SpatRaster")
    stop("The object 'P_data' must be a SpatRaster")
  
  # Extract dates from object
  dates <- terra::time(P_data)
  
  # Apply Deciles
  idx <- terra::app(P_data, .deciles)
    
  ## Set dates and return
  names(idx)        <- paste0(substr(dates, 1, 7)) 
  terra::time(idx)  <- dates
  
  # Avoid NaN values
  idx[is.nan(idx)]      <- NA
  
  return(idx)
  
}

#' Function to calculate the Percent of Normal Index
#'
#' @param P_data 'SpatRaster' object that contains spatially-distributed monthly precipitation data that will be used to calculate the 
#'  Percent of Normal Index. This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can 
#'  be set with the function time of the terra package.
#'
#' @return Spatially-distributed Percent of Normal Index.
#' @export
#'
#' @examples
spatial_pni <- function(P_data){
  
  # Check P_data
  if(class(P_data) != "SpatRaster")
    stop("The object 'P_data' must be a SpatRaster")
  
  # Extract dates from P_data object
  dates <- terra::time(P_data)
  
  # Apply PNI
  idx <- terra::app(P_data, .pni, dates = dates)
  
  
  ## set dates and return
  names(idx)        <- paste0(substr(dates, 1, 7)) 
  terra::time(idx)  <- dates
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  idx[is.infinite(idx)] <- NA
  
  return(idx)
  
}



#' Calculation of the parameters of a selected probability distribution for a selected day (SPI).
#'
#' @param Prod_data 'SpatRaster' object that contains spatially-distributed monthly data for a specific date (e.g., "%Y-02-28).
#'   This object contains the daily accumulated values of a specific day according to the scale provided. This object can be computed with
#'   the 'spi_agregate_daily' function.
#'   This 'SpatRaster' must only contain the days that corresponds to the specific selection.
#' @param trgt The day for which the function will be computed. The default is NULL, indicating that the last day of the data
#'  will be used to return the SPI values.
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
#' @param package Either 'SCI' or 'SPEI'. Should the SCI or SPEI package be used in the implementation?
#'
#' @return This function returns the fit parameters of the selected distribution according to a selected day. The input object 'Prod_data'
#'    has to be created before, which can be achieved with the 'spi.agregate_daily' function. 
#' @export
#'
#' @examples
params_spi <- function(Prod_data,  
                       trgt = NULL, 
                       ref_start = NULL,
                       ref_end = NULL, 
                       distribution = "Gamma", 
                       fit = "ub-pwm",
                       package = "SCI"){
  
  # Check Prod_data
  if(class(Prod_data) != "SpatRaster")
    stop("The object 'Prod_data' must be a SpatRaster")
  
  # Extract dates from Prod_data object
  dates <- terra::time(Prod_data)
  
  if(is.null(trgt))
    trgt <- dates[length(dates)]
  
  # Evaluate the 'trgt' element
  if(!is.null(trgt) & !as.Date(trgt) %in% dates)
    stop("The 'trgt' date should be included in 'Prod_data'.")
  
  # Check ref_start object
  if(class(ref_start) != "character" & !is.null(ref_start))
    stop("If object 'ref_start' is not set to NULL, it must be a character object that indicates the starting point of the reference period used for computing the SPEI. The format should be '%Y-%m'")
  
  # Check ref_end object
  if(class(ref_end) != "character" & !is.null(ref_end))
    stop("If object 'ref_end' is not set to NULL, it must be a character object that indicates the ending point of the reference period used for computing the SPEI. The format should be '%Y-%m'")
  
  # Check ref_start and ref_end
  if(is.null(ref_start) & !is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  if(!is.null(ref_start) & is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  # Checking the 'package' object
  if(!package %in% c("SPEI", "SCI"))
    stop("The 'package' object must be either 'SPEI' or 'SCI'")
  
  # Apply daily SPI
  if(package == "SPEI"){
    
    idx <- terra::app(Prod_data, .get_params.spei, trgt = trgt, dates = dates, 
                      ref_start = ref_start, ref_end =ref_end, distribution = distribution, 
                      fit = fit)
    
  } else {
    
    idx <- terra::app(Prod_data, .get_params.sci, trgt = trgt, dates = dates, 
                      ref_start = ref_start, ref_end =ref_end, distribution = distribution, 
                      fit = fit)
    
  }
  
  ## set dates and return
  terra::time(idx)  <- rep(as.Date(trgt), terra::nlyr(idx))
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  idx[is.infinite(idx)] <- NA
  
  # Adding names to the parameters
  names <- switch(distribution,
                  "Gamma" = c('shape','rate'),
                  "PearsonIII" = c('mu','sigma','gamma'),
                  "log-Logistic" = c('xi','alpha','kappa')
  )
  
  if(package == "SCI" & distribution == "log-Logistic")
    names = c("location",  "scale")
  
  if(package == "SCI" & distribution == "PearsonIII")
    names = c("shape",  "scale", "location")
  
  
  names(idx) <- names
  
  
  return(idx)
  
}

#' Calculation of the parameters of a selected probability distribution  for a selected day (SPEI).
#'
#' @param Prod_data 'SpatRaster' object that contains spatially-distributed monthly data for a specific date (e.g., "%Y-02-28).
#'   This object contains the daily accumulated values of a specific day according to the scale provided. This object can be computed with
#'   the 'spi_agregate_daily' function.
#'   This 'SpatRaster' must only contain the days that corresponds to the specific selection.
#' @param trgt The day for which the function will be computed. The default is NULL, indicating that the last day of the data
#'  will be used to return the SPI values.
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
#' @param package Either 'SCI' or 'SPEI'. Should the SCI or SPEI package be used in the implementation?
#'
#' @return This function returns the fit parameters of the selected distribution according to a selected day. The input object 'Prod_data'
#'    has to be created before, which can be achieved with the 'spi.agregate_daily' function. 
#' @export
#'
#' @examples
params_spei <- function(Prod_data,  
                       trgt = NULL, 
                       ref_start = NULL,
                       ref_end = NULL, 
                       distribution = "log-Logistic", 
                       fit = "ub-pwm",
                       package = "SCI"){
  
  
  idx <- params_spi(Prod_data, trgt = trgt, ref_start = ref_start, 
                   ref_end = ref_end, distribution = distribution, fit = fit, package = package)
  
  
  return(idx)
  
}

#' Calculation of the parameters of a selected probability distribution for all days (SPI and SPEI).
#'
#' @param P_data ''SpatRaster' object that contains spatially-distributed daily precipitation data that will be used to calculate the 
#'  accumulations according a selected scale. This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can 
#'  be set with the function time of the terra package.
#' @param scale Integer value that represents the time scale at which the SPI will be computed.
#' @param ref_start optional value that represents the starting point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.
#' @param ref_end Optional value that represents the ending point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.
#' @param distribution Optional value indicating the name of the distribution function to be used for calculating the parameters
#'  (one of 'log-Logistic', 'Gamma' and 'PearsonIII').
#' @param fit Optional value indicating the name of the method used for computing the distribution function parameters 
#'  (one of 'ub-pwm', 'pp-pwm' and 'max-lik'). Defaults to 'ub-pwm'.
#' @param package Either 'SCI' or 'SPEI'. Should the SCI or SPEI package be used in the implementation?
#'
#' @return This function returns the fit parameters of the selected distribution for all days. 
#' Please note that this function actively excludes the 29th of February for consistency reasons. The dates of the return 'SpatRaster'
#' represent the days for the calculated parameters. The last year in the computation is used in the dates (in case that ref_end exist, that year is used).
#' @export
#'
#' @examples
params_spi <- function(P_data,  
                       scale,
                       ref_start = NULL,
                       ref_end = NULL, 
                       distribution, 
                       fit = "ub-pwm",
                       package = "SCI"){
  
  # Setting values to be used in the iteration process
  dates   <- terra::time(P_data)
  periods <- substr(dates, 6, 10)
  periods <- unique(periods)
  
  # Excluding the 29th of February
  pos <- which(periods == "02-29")
  if(length(pos) == 1)
    periods <- periods[-pos]
  
  # Inisialising an empty list to store the results
  res <- list()
  
  cat("| Initialisation of iteration process. This may require a lot of time! | \n")
  
  # start of the iteration process
  for(i in 1:length(periods)){
    
    # Setting the 'trgt' day
    year <- max(as.numeric(substr(dates, 1, 4)))
    if(!is.null(ref_start))
      year <- as.numeric(substr(ref_start), 1, 4)
    
    target <- paste0(year, "-", periods[i])
    
    cat("Processing", periods[i], ":")
    
    # Computing the accumulations
    Prod_data <- aggregate_days4spi(P_data, scale = scale, trgt = target)
    
    # Calculating the parameters for the respective day
    res[[i]] <- params_spi(Prod_data, trgt = target, ref_start = ref_start, 
                        ref_end = ref_end, distribution = distribution, fit = fit, package = package)
    
    cat("Done! \n")
  }
  
  # Stacking the list
  res <-terra::rast(res)
  
  return(res)
  
}





#' Standardised precipitation index calculated for a specific day.
#'
#' @param Prod_data 'SpatRaster' object that contains spatially-distributed monthly data for a specific date (e.g., "%Y-02-28).
#'   This object contains the daily accumulated values of a specific day according to the scale provided. This object can be computed with
#'   the 'spi_agregate_daily' function.
#'   This 'SpatRaster' must only contain the days that corresponds to the specific selection.
#' @param trgt The day for which the function will be computed. The default is NULL, indicating that the last day of the data
#'  will be used to return the SPI values.
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
#'  @param params Should the parameters of the selected distributions be returned? Set to FALSE as the default.
#' @param package Either 'SCI' or 'SPEI'. Should the SCI or SPEI package be used in the implementation?
#'
#' @return This function returns one layer of SPI-n according to a selected day. The input object 'Prod_data'
#'    has to be created before, which can be achieved with the 'spi.agregate_daily' function. 
#' @export
#'
#' @examples
daily_spi <- function(Prod_data,  
                      trgt = NULL, 
                      ref_start = NULL,
                      ref_end = NULL, 
                      distribution = "Gamma", 
                      fit = "ub-pwm",
                      params = FALSE,
                      package = "SCI"){
  
  # Check Prod_data
  if(class(Prod_data) != "SpatRaster")
    stop("The object 'Prod_data' must be a SpatRaster")
  
  # Extract dates from Prod_data object
  dates <- terra::time(Prod_data)
  
  if(is.null(trgt))
    trgt <- dates[length(dates)]
  
  # Evaluate the 'trgt' element
  if(!is.null(trgt) & !as.Date(trgt) %in% dates)
    stop("The 'trgt' date should be included in 'Prod_data'.")
  
  # Check ref_start object
  if(class(ref_start) != "character" & !is.null(ref_start))
    stop("If object 'ref_start' is not set to NULL, it must be a character object that indicates the starting point of the reference period used for computing the SPEI. The format should be '%Y-%m'")
  
  # Check ref_end object
  if(class(ref_end) != "character" & !is.null(ref_end))
    stop("If object 'ref_end' is not set to NULL, it must be a character object that indicates the ending point of the reference period used for computing the SPEI. The format should be '%Y-%m'")
  
  # Check ref_start and ref_end
  if(is.null(ref_start) & !is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  if(!is.null(ref_start) & is.null(ref_end))
    stop("The objects 'ref_start' and 'ref_end' should be either both NULL or both character!")
  
  # Checking the 'package' object
  if(!package %in% c("SPEI", "SCI"))
    stop("The 'package' object must be either 'SPEI' or 'SCI'")
  
  # Apply daily SPI
  if(package == "SPEI"){
    
    idx <- terra::app(Prod_data, .spei_daily.spei, trgt = trgt, dates = dates, 
                      ref_start = ref_start, ref_end =ref_end, distribution = distribution, 
                      fit = fit, params = params)
    
  } else {
    
    idx <- terra::app(Prod_data, .spei_daily.sci, trgt = trgt, dates = dates, 
                      ref_start = ref_start, ref_end =ref_end, distribution = distribution, 
                      fit = fit, params = params)
    
  }
  
  ## set dates and return
  terra::time(idx)  <- rep(as.Date(trgt), terra::nlyr(idx))
  
  # Avoid NaNs and infinite values
  idx[is.nan(idx)]      <- NA
  idx[is.infinite(idx)] <- NA
  
  # Adding information in case that 'params' is set to true
  if(params){
    
    names <- switch(distribution,
                    "Gamma" = c('shape','rate'),
                    "PearsonIII" = c('mu','sigma','gamma'),
                    "log-Logistic" = c('xi','alpha','kappa')
    )
    
    
    names(idx) <- c("Values", names)
    
  } else {
    names(idx)  <- paste0(substr(trgt, 1, 7)) 
  }
  
  return(idx)
  
}


#' Standardised precipitation evapotranspiration index calculated for a specific day. 
#'  This function is a wrapper function of daily.spi
#'
#' @param Prod_data 'SpatRaster' object that contains spatially-distributed monthly data for a specific date (e.g., "%Y-02-28).
#'   This object contains the daily accumulated values of a specific day according to the scale provided. This object can be computed with
#'   the 'spi_agregate_daily' function.
#'   This 'SpatRaster' must only contain the days that corresponds to the specific selection.
#' @param trgt The day for which the function will be computed. The default is NULL, indicating that the last day of the data
#'  will be used to return the SPI values.
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
#'  @param params Should the parameters of the selected distributions be returned? Set to FALSE as the default.
#' @param package Either 'SCI' or 'SPEI'. Should the SCI or SPEI package be used in the implementation?
#'
#' @return
#' @export
#'
#' @examples
daily_spei <- function(Prod_data,  
                      trgt = NULL, 
                      ref_start = NULL,
                      ref_end = NULL, 
                      distribution = "log-Logistic", 
                      fit = "ub-pwm",
                      params = FALSE,
                      package = "SCI"){
  
 idx <- daily.spi(Prod_data, trgt = trgt, ref_start = ref_start, 
                  ref_end = ref_end, distribution = distribution, fit = fit, params = params)
  

  return(idx)
  
}




