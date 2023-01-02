################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: functions to calculate the SPI and SPEI for a specific date
################################################################################
##
## Creation date: 2023-01-01
##
################################################################################
################################################################################

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
calculate_params <- function(P_data,  
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
    
    # Calculating the probability of zero
    ref_period <- c(ref_start, ref_end)
    if(is.null(ref_period))
      ref_period <- "None"
    pze        <- terra::app(Prod_data, .calculate_pze, ref_period = ref_period, dates = dates)
    
    # Calculating the parameters for the respective day
    res[[i]] <- params_spi(Prod_data, trgt = target, ref_start = ref_start, 
                           ref_end = ref_end, distribution = distribution, fit = fit, package = package)
    
    cat(" Done! \n")
  }
  
  # Stacking the list
  res <-terra::rast(res)
  
  # Getting the names of the parameters
  param_names <- unique(names(res))
  
  sep_params <- list()
  for(i in 1:length(param_names)){
    
    sep_params[[i]] <- res[[which(names(res) %in% param_names[i])]]
    
  }
  
  # Setting the manes of the list
  names(sep_params) <- param_names
  
  # Setting additional information that will be required in the daily_spi or daily_spei function
  sep_params$probability_zero <- pze
  sep_params$distribution     <- distribution
  sep_params$package          <- package
  
  
  class(sep_params) <- "params_list"
  
  return(sep_params)
  
}


#' Kalman smoothing to daily SPI and SPEI parameters
#'
#' @param params_list List object that contain as many 'SparRaster' objects as parameters in the selected distribution.
#'    This list must contain the names of the parameters. Please see the 'calculate_params' function.
#' @param H Covariance matrix or array of disturbance terms \epsilon_tϵ of observation equations. See the KFAS package.
#' @param ... 
#'
#' @return The smoothed parameters in a list.
#' @export
#'
#' @examples
kalman_parameters <- function(params_list, H, ...){
  
  # Checking if the class of the object is 'params_list'
  if(!class(params_list) == "params_list")
    stop("The 'params_list' object must have a 'params_list' class. Please see the 'calculate_params' function.")
  
  # Checking the object H
  if(!is.numeric(H))
    stop("The parameter 'H' must be numeric")
  
  # Checking length of H
  if(!(length(H) == 1 | length(H) == length(params_list)))
    stop("The length of the parameter 'H' must be either 1 or equal to the number of parameters")
  
  # Checking the position of the distribution and package in the list
  pos_attributes <- which(names(params_list) %in% c("distribution", "package"))
  pos_pze        <- which(names(params_list) %in% "probability_zero") 
  
  # Storing the attributes in an object and excluding them from the 'params_list' object
  attributes     <- params_list[pos_attributes]
  pze            <- params_list[pos_pze]
  params_list    <- params_list[-c(pos_attributes, pos_pze)]
  
  # Applying the filter for every parameter
  parameters <- names(params_list)
  dates      <- terra::time(params_list[[1]])
  res        <- list()
  
  for(i in 1:length(parameters)){
    
    # Checking parameter H
    H_iter <- H[1]
    if(length(H) == length(params_list))
      H_iter <- H[i]
    
    res[[i]]              <- terra::app(params_list[[i]], .kalman_filter, H = H_iter, ...)
    terra::time(res[[i]]) <- dates
    names(res[[i]])       <- dates
    
  }
  
  # Setting the names to the smoothed parameters
  names(res) <- parameters
  
  # Passing the attributes to the final list
  res$probability_zero <- pze
  res$distribution <- attributes$distribution
  res$package      <- attributes$package
  
  class(res) <- "params_list"
  
  return(res)
  
}

#' Standardised precipitation index calculated for a specific day.
#'
#'@param param_list A list object that contains a 'SpatRaster' of 365 layers for each parameter of the selected distribution.
#' @param P_lyr 'SpatRaster' object that contains spatially-distributed monthly data for a specific date (e.g., "%Y-02-28).
#'   This object contains the daily accumulated values of a specific day according to the scale provided. This object can be computed with
#'   the 'spi_agregate_daily' function.
#'   This 'SpatRaster' must only contain the days that corresponds to the specific selection.
#'
#' @return This function returns one layer of SPI-n according to a selected day. The input object 'Prod_data'
#'    has to be created before, which can be achieved with the 'spi.agregate_daily' function. 
#' @export
#'
#' @examples
#' 
daily_spi <- function(param_list,
                      P_lyr){
  
  # Checking if the class of the object is 'params_list'
  if(!class(params_list) == "params_list")
    stop("The 'params_list' object must have a 'params_list' class. Please see the 'calculate_params' function.")
  
  # Checking the resolution of P_lyr
  if(!terra::compareGeom(P_lyr, param_list[[1]]))
    stop("The gridded parameters in 'params_list' and the object 'P_lyr' do not have the same raster geometry.")
  
  # Checking the position of the distribution and package in the list
  pos_attributes <- which(names(params_list) %in% c("distribution", "package"))
  pos_pze        <- which(names(params_list) %in% c("probability_zero"))
  
  # Storing the attributes in an object and excluding them from the 'params_list' object
  pze            <- params_list[pos_pze]
  attributes     <- params_list[pos_attributes]
  params_list    <- params_list[-c(pos_attributes, pos_pze)]
  
  # Converting the 'SpatRaster' object in matrices
  params_matrix <- list()
  for(i in 1:length(params_list))
    params_matrix[[i]] <- as.matrix(params_list[[i]])
  
  names(params_matrix) <- names(params_list)
  
  # Converting the P_lyr and pze to matrix
  P_matrix   <- as.matrix(P_lyr)
  pze_matrix <- as.matrix(pze)
  
  # Setting target day, dates, and used reference period
  trgt       <- terra::time(P_lyr)
  
  # Implementation of the SPEI package ('sbegueria' https://github.com/sbegueria/SPEI/blob/master/R/spei.R)
  if(attributes$package == "SPEI"){
    
    result <- .spei_daily.spei(params_matrix, P_matrix, pze_matrix,
                               trgt = trgt, distribution = attributes$distribution)
    
    
    # Implementation of the SCI package ( https://github.com/cran/SCI/blob/master/R/sci.r)
  } else {
    
    result <- .spei_daily.sci(params_matrix, P_matrix, pze_matrix,
                              trgt = trgt, distribution = attributes$distribution)
    
  }
  
  # converting it into a 'SpatRaster'
  final_layer <- param_list[[1]][[1]]
  terra::values(final_layer) <- result
  terra::time(final_layer)   <- trgt
  names(final_layer)         <- trgt
  
  return(final_layer)
}

#' Standardised precipitation evapotranspiration index calculated for a specific day. 
#'  This function is a wrapper function of daily.spi
#'
#'@param param_list A list object that contains a 'SpatRaster' of 365 layers for each parameter of the selected distribution.
#' @param P_lyr 'SpatRaster' object that contains spatially-distributed monthly data for a specific date (e.g., "%Y-02-28).
#'   This object contains the daily accumulated values of a specific day according to the scale provided. This object can be computed with
#'   the 'spi_agregate_daily' function.
#'   This 'SpatRaster' must only contain the days that corresponds to the specific selection.
#'
#' @return This function returns one layer of SPI-n according to a selected day. The input object 'Prod_data'
#'    has to be created before, which can be achieved with the 'spi.agregate_daily' function. 
#' @export
#'
#' @examples
#' 
daily_spi <- function(param_list,
                      P_lyr){
  
  idx <- daily_spi(param_list, P_lyr)
  
  
  return(idx)
  
}