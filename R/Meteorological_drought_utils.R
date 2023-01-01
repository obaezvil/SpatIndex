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
    distribution <- 'logis'
  
  if(distribution == "PearsonIII")
    distribution <- 'pe3'
  
  # Obtaining parameters for the reference period
  spi_ref_params <- SCI::fitSCI(x[pos_ini:pos_fin], first.mon = 1, 
                                distr = distribution, 
                                time.scale = scale, p0 = TRUE)
  
  # Apply the the parameters of the reference period to all data
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
    distribution <- 'logis'
  
  if(distribution == "PearsonIII")
    distribution <- 'pe3'
  
  # Obtaining parameters for the reference period
  spei_ref_params <- SCI::fitSCI(x[pos_ini:pos_fin], first.mon = 1, 
                                distr = distribution, 
                                time.scale = scale, p0 = FALSE)
  
  # Apply the the parameters of the reference period to all data
  spei <- SCI::transformSCI(x, first.mon = 1, obj = spei_ref_params)
  
  # Trasforming NaNs to NAs
  spei[is.nan(spei)] <- NA
  
  return(spei)
  
}

#' Utils function to calculate deciles
#'
#' @param x Numerical vector.
#'
#' @return Numerical vector with the corresponding deciles
#' @export
#'
#' @examples
.deciles <- function(x){
  
  # Creating zoo file
  x <- zoo::zoo(x, dates)
  
  # Extracting the month from the dates
  mt        <- as.numeric(substr(dates, 6, 7))
  mt_unique <- unique(mt)
  
  # Function to get quantile
  .get.decile <- function(a, deciles){
    
    euclidean <- sqrt((a - deciles) ^ 2)
    pos       <- which.min(euclidean)
    return(pos)
    
  }
  
  # Creating a new empty vector
  result <- c()
  
  # Creating a loop for iterating for each month
  for(i in 1:length(mt_unique)){
    
    # Storing the position of the month and subsetting the respective month
    pos_month <- which(mt %in% mt[i])
    x_month  <- x[pos_month]
    
    # Calculating deciles from a vector
    deciles <- quantile(x_month, probs = seq(0, 1, by = 0.1), na.rm = TRUE)
    
    # Obtaining mean value between deciles
    deciles <- zoo::rollapply(deciles, width = 2, mean, 
                              align = "left", partial = TRUE, na.rm = TRUE)[-11]
    
    # Apply function to all time series
    dec <- unlist(lapply(as.numeric(x_month), .get.decile, deciles))
    
    # Storing the values in the respective month
    result[pos_month] <- dec
  }
  
  # Replace with NAs over areas with full NA vectors
  if(length(result) != length(x))
    result <- rep(NA, length(x))
  
  return(result)
  
}

#' Utils function to calculate percent of normal 
#'
#' @param x Numerical vector.
#' @param dates Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_spei' function.
#'
#' @return Numerical vector with the corresponding percent of normal precipitation
#' @export
#'
#' @examples
.pni <- function(x, dates){
  
  # Creating zoo file
  x <- zoo::zoo(x, dates)
  
  # Aggregating mean monthly values
  monthly <- hydroTSM::monthlyfunction(x, FUN = mean)
  
  # Extracting the month from the dates
  mt <- as.numeric(substr(dates, 6, 7))
  
  # Function to compute the pni for each value
  .get.pni <- function(a, pos, monthly){
    
    r <- (as.numeric(a) / as.numeric(monthly)[pos]) * 100
    
    return(r)
    
  }
  
  # Apply function to all time series
  pni <- unlist(mapply(.get.pni, x, mt, MoreArgs = list(monthly)))
  
  # Replace with NAs over areas with full NA vectors
  if(length(pni) != length(x))
    pni <- rep(NA, length(x))
  
  return(pni)
  
}


#' Accumulations to be used by the daily_spi and daily_spei functions, 
#'  
#' @param Prod_data 'SpatRaster' object that contains spatially-distributed monthly precipitation data that will be used to calculate the SPI. 
#'  This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
#'  of the terra package.
#' @param scale Integer value that represents the time scale at which the SPI will be computed.
#' @param trgt The day for which the function will be computed. The default is NULL, indicating that the last day of the data
#'  will be selected.
#' @return This function returns a 'SpatRaster' of the specific date (e.g., "%Y-02-28") according to the defined scale, defined as a 30-day 
#'   accumulation period (e.g., scale = 3, accumulates the specific date plus the previous 89 days).
#' @export
#'
#' @examples
aggregate_days4spi <- function(Prod_data, 
                               scale, 
                               trgt = NULL){
  
  # Aggregating the time series according to the scale specified
  #   One month represents the previous 30 days
  no_days <- scale * 30
  
  # Extracting the dates from 'P_data'
  dates <- terra::time(Prod_data)
  
  # If the target date is not specified, the last day is selected as the target
  if(is.null(trgt))
    trgt <- dates[length(dates)]
  
  # extracting the desired period to compute the accumulations
  period <- substr(trgt, 6, 10)
  
  # Matching the dates of multiple years that are equal to 'period' and calculating the
  #   initial positions where the accumulation of days should begin (to calculate the monthly 
  #   P according tho the defined scale)
  pos_fin   <- grep(period, dates)
  pos_ini   <- (pos_fin - no_days) + 1
  dates_fin <- dates[pos_fin]
  
  # Avoiding negative indices. If the accumulation canot be calculated due to insufficient data
  #  the layer will be neglected in the following step
  pos_ini[pos_ini < 0] <- 0
  
  # If a year is not complete then we neglect it (this might happen in the first years depending on scale)
  periods          <- pos_fin - pos_ini
  incomplete_years <- which(periods != (no_days - 1))
  if(length(incomplete_years) > 0){
    pos_fin   <- pos_fin[-incomplete_years]
    pos_ini   <- pos_ini[-incomplete_years]
    dates_fin <- dates_fin[-incomplete_years]
  }
  
  # Starting an iterative process to calculate the accumulations per year
  accums <- list()
  for(i in 1:length(dates_fin)){
    
    positions <- pos_ini[i]:pos_fin[i]
    rst      <- Prod_data[[positions]]
    accums[[i]] <- sum(rst, na.rm = TRUE)
    
  }
  
  # Converting the data into a 'SpatRaster'
  accums              <- terra::rast(accums)
  terra::time(accums) <- dates_fin
    
  return(accums)
  
}



#' Utils function to calculate the parameters of the selected distribution function for the SPI or SPEI using a vector. This function is based on the code presented in the SPEI package
#'  from 'sbegueria' https://github.com/sbegueria/SPEI/blob/master/R/spei.R
#'
#' @param x Numerical vector.
#' @param trgt The day for which the function will be computed. The default is NULL, indicating that the last day of the data
#'  will be selected.
#' @param dates Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_spei' function.
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
#'
#' @return
#'
#' @examples
.get_params.spei <- function(x, 
                             trgt,
                             dates, 
                             ref_start,
                             ref_end,
                             distribution,
                             fit){
  
  # Checking distribution
  if(!distribution %in% c("Gamma", "log-Logistic", "PearsonIII"))
    stop("The accepted distributions are 'Gamma', 'log-Logistic', and 'PearsonIII'")
  
  # Checking fit
  if(!fit %in% c('ub-pwm', 'pp-pwm', 'max-lik'))
    stop("The accepted fit  are 'ub-pwm', 'pp-pwm', and 'max-lik'")
  
  coef <- switch(distribution,
                 "Gamma" = array(NA, c(2, 1, 1),
                                 list(par=c('alpha','beta'), colnames(x), NULL)),
                 "PearsonIII" = coef <- array(NA, c(3, 1, 1),
                                              list(par=c('mu','sigma','gamma'), colnames(x), NULL)),
                 "log-Logistic" = array(NA, c(3, 1, 1),
                                        list(par=c('xi','alpha','kappa'), colnames(x), NULL)),
                 "GEV" = array(NA, c(3, 1, 1),
                               list(par=c('xi','alpha','kappa'), colnames(x), NULL))
  )
  
  # converting into a ts
  acu <- zoo::zoo(x, as.Date(dates))
  
  # converting NaNs and Inf values to NAs
  acu[which(is.nan(acu))]      <- NA
  acu[which(is.infinite(acu))] <- NA
  
  # Trim data set to reference period for fitting (acu_ref)
  if(!is.null(ref_start) & !is.null(ref_end)){
    
    pos_ini <- min(which(as.Date(paste0(ref_start, "-01")) <= as.Date(dates)))
    pos_fin <- max(which(as.Date(paste0(ref_end, "-01")) >= as.Date(dates)))
    
    acu_ref <- acu[pos_ini:pos_fin]
    
  } else {
    
    acu_ref <- acu
    
  }
  
  # Start a object to store the results
  res <- c()
  ################################
  
  # Analysing the NA values
  f     <- acu_ref[!is.na(acu_ref)]
  ff    <- acu[!is.na(acu)]
  x_mon <- f
  
  # Probability of zero (pze)
  if(distribution != 'log-Logistic' & length(na.omit(x_mon)) > 0){
    pze   <- sum(x_mon==0) / length(x_mon)
    x_mon <- x_mon[x_mon > 0]
  }
  
  # Distribution parameters (f_params)
  #  Fit distribution parameters
  x_mon_sd <- sd(x_mon, na.rm=TRUE)
  
  #####################################################
  # Condition to evaluate the length of values in x_mon                                               
  if (length(x_mon) < 4) {
    
    val      <- NA
    f_params <- coef[,,1]
    
  } else if (is.na(x_mon_sd) || (x_mon_sd == 0)){
    
    val      <- NA
    f_params <- coef[,,1]
    
  } else if (length(na.omit(x_mon)) == 0){
    
    val      <- NA
    f_params <- coef[,,1]
    
  } else {
    
    
    # Calculate probability weighted moments based on `lmomco` or `TLMoments`
    pwm <- switch(fit,
                  'pp-pwm' = lmomco::pwm.pp(as.numeric(x_mon), -0.35, 0, nmom=3, sort=TRUE),
                  'ub-pwm' = TLMoments::PWM(as.numeric(x_mon), order=0:2)
    )
    
    # Check L-moments validity
    lmom <- lmomco::pwm2lmom(pwm)
    
    # `lmom` fortran functions need specific inputs L1, L2, T3
    # This is handled internally by `lmomco` with `lmorph`
    fortran_vec <- c(lmom$lambdas[1:2], lmom$ratios[3])
    
    # Calculate parameters based on distribution with `lmom`, then `lmomco`
    f_params <- switch(distribution,
                       'log-Logistic' = tryCatch(lmom::pelglo(fortran_vec),
                                                 error = function(e){ lmomco::parglo(lmom)$para }),
                       'Gamma' = tryCatch(lmom::pelgam(fortran_vec),
                                          error = function(e){ lmomco::pargam(lmom)$para }),
                       'PearsonIII' = tryCatch(lmom::pelpe3(fortran_vec),
                                               error = function(e){ lmomco::parpe3(lmom)$para })
    )
    
    # Adjust if user chose `log-Logistic` and `max-lik`
    if(distribution == 'log-Logistic' && fit=='max-lik'){
      f_params <- SPEI::parglo.maxlik(x.mon, f_params)$para
    }
    

    # Check L-moments validity
    if ( !lmomco::are.lmom.valid(lmom) || anyNA(lmom[[1]]) || any(is.nan(lmom[[1]])) ){ 
      f_params <- rep(NA, length(f_params))
    }
    
  }
  
  # Converting the parameter beta into the rate parameter (rate = 1 / Beta)
  if(distribution == "Gamma")
    f_params[2] <- 1 /f_params[2]
    
  return(f_params)
  
}

#' Utils function to calculate the parameters of the selected distribution function for the SPI or SPEI using a vector. This function is based on the code presented in the SCI package
#'  obtained from  https://github.com/cran/SCI/blob/master/R/sci.r
#'
#' @param x Numerical vector.
#' @param trgt The day for which the function will be computed. The default is NULL, indicating that the last day of the data
#'  will be selected.
#' @param dates Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_spei' function.
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
#' @param params Should the parameters of the selected distributions be returned? Set to FALSE as the default.
#'
#' @return
#'
#' @examples
.get_params.sci <- function(x, 
                            trgt,
                            dates, 
                            ref_start,
                            ref_end,
                            distribution,
                            fit,
                            params){
  
  # Saving the raw data in an object
  x_all <- x
  
  # Checking distribution
  if(!distribution %in% c("Gamma", "log-Logistic", "PearsonIII"))
    stop("The accepted distributions for the SCI package are 'Gamma', 'log-Logistic', and 'PearsonIII'")
  
  # Adjusting distributions for SCI package
  if(distribution == "Gamma")
    distribution <- tolower(distribution)
  
  if(distribution == "log-Logistic")
    distribution <- 'logis'
  
  if(distribution == "PearsonIII")
    distribution <- 'pe3'
  
  # model Probability of zero (precipitation) months is modeled with a mixed distribution as 
  #   D(x) = p_0 + (1-p_0)G(x)
  p0 = TRUE
  
  # Trim data set to reference period for fitting (acu_ref)
  if(!is.null(ref_start) & !is.null(ref_end)){
    
    pos_ini <- min(which(as.Date(paste0(ref_start, "-01")) <= as.Date(dates)))             
    pos_fin <- max(which(as.Date(paste0(ref_end, "-01")) >= as.Date(dates)))
    
    acu_ref <- x[pos_ini:pos_fin]
    
  } else {
    
    acu_ref <- x
    
  }
  
  # Converting the data to numeric
  x <- as.numeric(acu_ref)   
  
  # Setting attributes of the data
  nn         <- length(x)
  time.index <- 1:nn
  
  # Find the distibution parameters for the data
  x.fit                <- vector("list", 1)
  names(x.fit)         <- paste("M",1,sep="")
  x.fit.monitor        <- 0
  names(x.fit.monitor) <- names(x.fit)
  
  empty.fit <- try(SCI::dist.start(NA,distribution),silent=TRUE)
  if(class(empty.fit)=="try-error"){
    stop("distribution ", distribution," not defined for start.fun=",deparse(substitute(start.fun)),
         "\n or startfun not implemented correctly")
  } else {
    empty.fit <- unlist(empty.fit)
  }
  
  if(p0){
    empty.fit <- c(empty.fit,P0=NA)
  }
  
  # Setting an empty list for the mledist
  mledist.par <- list()
  mledist.par$distr <- distribution
  
  # For the specific month
  
  ## Find distribution parameter for each month
  mledist.par$data <- x
  mledist.par$data <- mledist.par$data[is.finite(mledist.par$data)] # condition to neglect non-finite values
  
  # If there are not finite values in the data
  if(length(mledist.par$data) == 0){
    
    x.fit[[1]]       <- empty.fit
    x.fit.monitor[1] <- 4
    
  } else if(all(mledist.par$data==mledist.par$data[1])){ ## if all values in calibration period are eaqual...
    
    x.fit[[1]]       <- empty.fit
    x.fit.monitor[1] <- 5
    
  } else {
    
    np0    <- sum(mledist.par$data==0)
    nn     <- length(mledist.par$data)
    p0.est <- np0/nn
    
    # Conditional if there are zeroes in the time series
    if(p0.est > 0){
      
      mledist.par$data <- mledist.par$data[mledist.par$data > 0]
      ## Catch that adds a single value close to zero if there are zeros in the
      ## fitting period.  This forces the distribution to tend towards zero,
      ## preventing a "gap" betwen 0 and data.
      mledist.par$data <- c(mledist.par$data, 0.01 * min(mledist.par$data, na.rm=TRUE))
      
    }
    
    mledist.par$start <- SCI::dist.start(x=mledist.par$data, distr = distribution)
    fail.value        <- mledist.par$start
    fail.value        <- unlist(fail.value)
    fail.value[] <- NA
    
    # Checking if there are NA values in the shape and rate
    if(any(is.na(unlist(mledist.par$start)))){
      
      x.fit[[1]]        <- fail.value
      x.fit.monitor[1] <- 1
      
    } else {
      
      sink <- capture.output(
        d.fit <- try(suppressWarnings(do.call(fitdistrplus::mledist, mledist.par)), silent=TRUE)
      )
      
      # Saving values
      if(class(d.fit)=="try-error"){
        
        x.fit[[1]]       <- fail.value
        x.fit.monitor[1] <- 2
        
      } else if(d.fit$convergence > 0){
        
        x.fit[[1]]       <- fail.value
        x.fit.monitor[1] <- 3
        
      } else {
        
        x.fit[[1]] <- d.fit$estimate
        
      }
      
      if(p0){
        
        x.fit[[1]] <- c(x.fit[[1]], P0 = p0.est)
        
      }
      
      # Condition to assess whether any value in x.fit is NA
      if(any(is.na(x.fit[[1]]))){
        
        x.fit[[1]][] <- NA
        
      }
      
    }
    
  }
  
  # Storing the fitted parameters
  f_params  <- do.call(cbind, x.fit)
  pos       <- which(rownames(f_params) == "P0")
  if(length(pos) > 0)
    f_params <- f_params[-pos,]
    
  return(f_params)
  
  
}

#' Utils function to compute the SPI or SPEI using a vector. This function is based on the code presented in the SPEI package
#'  from 'sbegueria' https://github.com/sbegueria/SPEI/blob/master/R/spei.R
#'
#' @param x Numerical vector.
#' @param trgt The day for which the function will be computed. The default is NULL, indicating that the last day of the data
#'  will be selected.
#' @param dates Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_spei' function.
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
#'
#' @return
#'
#' @examples
.spei_daily.spei <- function(x, 
                             trgt,
                             dates, 
                             ref_start,
                             ref_end,
                             distribution,
                             fit,
                             params){   #### RECONVERT BETA IN GAMMA
  
  # Checking distribution
  if(!distribution %in% c("Gamma", "log-Logistic", "PearsonIII"))
    stop("The accepted distributions are 'Gamma', 'log-Logistic', and 'PearsonIII'")
  
  # Checking fit
  if(!fit %in% c('ub-pwm', 'pp-pwm', 'max-lik'))
    stop("The accepted fit  are 'ub-pwm', 'pp-pwm', and 'max-lik'")
  
  coef <- switch(distribution,
                 "Gamma" = array(NA, c(2, 1, 1),
                                 list(par=c('alpha','beta'), colnames(x), NULL)),
                 "PearsonIII" = coef <- array(NA, c(3, 1, 1),
                                              list(par=c('mu','sigma','gamma'), colnames(x), NULL)),
                 "log-Logistic" = array(NA, c(3, 1, 1),
                                        list(par=c('xi','alpha','kappa'), colnames(x), NULL)),
                 "GEV" = array(NA, c(3, 1, 1),
                               list(par=c('xi','alpha','kappa'), colnames(x), NULL))
  )
  
  # converting into a ts
  acu <- zoo::zoo(x, as.Date(dates))
  
  # converting NaNs and Inf values to NAs
  acu[which(is.nan(acu))]      <- NA
  acu[which(is.infinite(acu))] <- NA
  
  # Trim data set to reference period for fitting (acu_ref)
  if(!is.null(ref_start) & !is.null(ref_end)){
    
    pos_ini <- min(which(as.Date(paste0(ref_start, "-01")) <= as.Date(dates)))
    pos_fin <- max(which(as.Date(paste0(ref_end, "-01")) >= as.Date(dates)))
    
    acu_ref <- acu[pos_ini:pos_fin]
    
  } else {
    
    acu_ref <- acu
    
  }
  
  # Start a object to store the results
  res <- c()
  ################################
  
  # Analysing the NA values
  f     <- acu_ref[!is.na(acu_ref)]
  ff    <- acu[!is.na(acu)]
  x_mon <- f
  
  # Probability of zero (pze)
  if(distribution != 'log-Logistic' & length(na.omit(x_mon)) > 0){
    pze   <- sum(x_mon==0) / length(x_mon)
    x_mon <- x_mon[x_mon > 0]
  }
  
  # Distribution parameters (f_params)
  #  Fit distribution parameters
  x_mon_sd <- sd(x_mon, na.rm=TRUE)
  
  #####################################################
  # Condition to evaluate the length of values in x_mon                                               
  if (length(x_mon) < 4) {
    
    val      <- NA
    f_params <- coef[,,1]
    
  } else if (is.na(x_mon_sd) || (x_mon_sd == 0)){
    
    val      <- NA
    f_params <- coef[,,1]
    
  } else if (length(na.omit(x_mon)) == 0){
    
    val      <- NA
    f_params <- coef[,,1]
    
  } else {
    
    
    # Calculate probability weighted moments based on `lmomco` or `TLMoments`
    pwm <- switch(fit,
                  'pp-pwm' = lmomco::pwm.pp(as.numeric(x_mon), -0.35, 0, nmom=3, sort=TRUE),
                  'ub-pwm' = TLMoments::PWM(as.numeric(x_mon), order=0:2)
    )
    
    # Check L-moments validity
    lmom <- lmomco::pwm2lmom(pwm)
    
    # `lmom` fortran functions need specific inputs L1, L2, T3
    # This is handled internally by `lmomco` with `lmorph`
    fortran_vec <- c(lmom$lambdas[1:2], lmom$ratios[3])
    
    # Calculate parameters based on distribution with `lmom`, then `lmomco`
    f_params <- switch(distribution,
                       'log-Logistic' = tryCatch(lmom::pelglo(fortran_vec),
                                                 error = function(e){ lmomco::parglo(lmom)$para }),
                       'Gamma' = tryCatch(lmom::pelgam(fortran_vec),
                                          error = function(e){ lmomco::pargam(lmom)$para }),
                       'PearsonIII' = tryCatch(lmom::pelpe3(fortran_vec),
                                               error = function(e){ lmomco::parpe3(lmom)$para })
    )
    
    # Adjust if user chose `log-Logistic` and `max-lik`
    if(distribution == 'log-Logistic' && fit=='max-lik'){
      f_params <- SPEI::parglo.maxlik(x.mon, f_params)$para
    }
    
    coef[,,1] <- f_params
    
    # Calculate CDF on `acu` using `f_params`
    cdf_res <- switch(distribution,
                      'log-Logistic' = lmom::cdfglo(acu, f_params),
                      'Gamma' = lmom::cdfgam(acu, f_params),
                      'PearsonIII' = lmom::cdfpe3(acu, f_params)
    )
    
    # Store the standardized values
    res <- qnorm(cdf_res)
    
    # Adjust for `pze` if distribution is Gamma or PearsonIII
    if(distribution == 'Gamma' | distribution == 'PearsonIII'){ 
      res <- qnorm(pze + (1-pze) * pnorm(res))
    }
    
    
    # Returning the result for the particular day
    if(!is.null(trgt)){
      
      val <- res[which(trgt == dates)]
      
      
    } else {
      
      val <- res[length(res)]
      
    }
    
    
    # Check L-moments validity
    if ( !lmomco::are.lmom.valid(lmom) || anyNA(lmom[[1]]) || any(is.nan(lmom[[1]])) ){ 
      val <- NA
    }
    
  }
  
  # Retrieving the parameters of the distribution
  if(params){
    
    # Converting the parameter beta into the rate parameter (rate = 1 / Beta)
    if(distribution == "Gamma")
      f_params[2] <- 1 /f_params[2]
    
    val <- c(as.numeric(val), f_params)
    
  }
  
  return(val)
  
}


#' Utils function to compute the SPI or SPEI using a vector. This function is based on the code presented in the SCI package
#'  obtained from  https://github.com/cran/SCI/blob/master/R/sci.r
#'
#' @param x Numerical vector.
#' @param trgt The day for which the function will be computed. The default is NULL, indicating that the last day of the data
#'  will be selected.
#' @param dates Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_spei' function.
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
#' @param params Should the parameters of the selected distributions be returned? Set to FALSE as the default.
#'
#' @return
#'
#' @examples
.spei_daily.sci <- function(x, 
                            trgt,
                            dates, 
                            ref_start,
                            ref_end,
                            distribution,
                            fit,
                            params){
  
  # Saving the raw data in an object
  x_all <- x
  
  # Checking distribution
  if(!distribution %in% c("Gamma", "log-Logistic", "PearsonIII"))
    stop("The accepted distributions for the SCI package are 'Gamma', 'log-Logistic', and 'PearsonIII'")
  
  # Adjusting distributions for SCI package
  if(distribution == "Gamma")
    distribution <- tolower(distribution)
  
  if(distribution == "log-Logistic")
    distribution <- 'logis'
  
  if(distribution == "PearsonIII")
    distribution <- 'pe3'
  
  # model Probability of zero (precipitation) months is modeled with a mixed distribution as 
  #   D(x) = p_0 + (1-p_0)G(x)
  p0 = TRUE
  
  # Trim data set to reference period for fitting (acu_ref)
  if(!is.null(ref_start) & !is.null(ref_end)){
    
    pos_ini <- min(which(as.Date(paste0(ref_start, "-01")) <= as.Date(dates)))             
    pos_fin <- max(which(as.Date(paste0(ref_end, "-01")) >= as.Date(dates)))
    
    acu_ref <- x[pos_ini:pos_fin]
    
  } else {
    
    acu_ref <- x
    
  }
  
  # Converting the data to numeric
  x <- as.numeric(acu_ref)   
  
  # Setting attributes of the data
  nn         <- length(x)
  time.index <- 1:nn
  
  # Find the distibution parameters for the data
  x.fit                <- vector("list", 1)
  names(x.fit)         <- paste("M",1,sep="")
  x.fit.monitor        <- 0
  names(x.fit.monitor) <- names(x.fit)
  
  empty.fit <- try(SCI::dist.start(NA,distribution),silent=TRUE)
  if(class(empty.fit)=="try-error"){
    stop("distribution ", distribution," not defined for start.fun=",deparse(substitute(start.fun)),
         "\n or startfun not implemented correctly")
  } else {
    empty.fit <- unlist(empty.fit)
  }
  
  if(p0){
    empty.fit <- c(empty.fit,P0=NA)
  }
  
  # Setting an empty list for the mledist
  mledist.par=list()
  mledist.par$distr <- distribution
  
  # For the specific month
  
  ## Find distribution parameter for each month
  mledist.par$data <- x
  mledist.par$data <- mledist.par$data[is.finite(mledist.par$data)] # condition to neglect non-finite values
  
  # If there are not finite values in the data
  if(length(mledist.par$data) == 0){
    
    x.fit[[1]]       <- empty.fit
    x.fit.monitor[1] <- 4
    
  } else if(all(mledist.par$data==mledist.par$data[1])){ ## if all values in calibration period are eaqual...
    
    x.fit[[1]]       <- empty.fit
    x.fit.monitor[1] <- 5
    
  } else {
    
    np0    <- sum(mledist.par$data==0)
    nn     <- length(mledist.par$data)
    p0.est <- np0/nn
    
    # Conditional if there are zeroes in the time series
    if(p0.est > 0){
      
      mledist.par$data <- mledist.par$data[mledist.par$data > 0]
      ## Catch that adds a single value close to zero if there are zeros in the
      ## fitting period.  This forces the distribution to tend towards zero,
      ## preventing a "gap" betwen 0 and data.
      mledist.par$data <- c(mledist.par$data, 0.01 * min(mledist.par$data, na.rm=TRUE))
      
    }
    
    mledist.par$start <- SCI::dist.start(x=mledist.par$data, distr = distribution)
    fail.value        <- mledist.par$start
    fail.value        <- unlist(fail.value)
    fail.value[] <- NA
    
    # Checking if there are NA values in the shape and rate
    if(any(is.na(unlist(mledist.par$start)))){
      
      x.fit[[1]]        <- fail.value
      x.fit.monitor[1] <- 1
      
    } else {
      
      sink <- capture.output(
        d.fit <- try(suppressWarnings(do.call(fitdistrplus::mledist, mledist.par)), silent=TRUE)
      )
      
      # Saving values
      if(class(d.fit)=="try-error"){
        
        x.fit[[1]]       <- fail.value
        x.fit.monitor[1] <- 2
        
      } else if(d.fit$convergence > 0){
        
        x.fit[[1]]       <- fail.value
        x.fit.monitor[1] <- 3
        
      } else {
        
        x.fit[[1]] <- d.fit$estimate
        
      }
      
      if(p0){
        
        x.fit[[1]] <- c(x.fit[[1]], P0 = p0.est)
        
      }
      
      # Condition to assess whether any value in x.fit is NA
      if(any(is.na(x.fit[[1]]))){
        
        x.fit[[1]][] <- NA
        
      }
      
    }
    
  }
  
  obj <- list(dist.para = do.call(cbind, x.fit),
              dist.para.flag = x.fit.monitor,
              time.scale = 1,
              distr = distribution,
              p0 = p0,
              p0.center.mass = FALSE,
              scaling = 1,
              call=match.call())
  
  class(obj) <- "fitSCI"
  f_params   <- obj$dist.para
  
  ####
  # Predicting the values based on the op object (probability based on fitted
  #   distribution)
  pdistr <- obj$distr
  pdistr <- match.fun(paste("p", pdistr, sep = ""))
  
  npar    <- nrow(obj$dist.para)
  PP0     <- obj$dist.para["P0",]
  distpar <- obj$dist.para[-npar,]
  
  # Getting the values for the particular moth
  xm <- as.numeric(x_all)
  
  if(any(is.na(distpar))){
    
    xm[] <- NA
    
  } else {
    
    xm <- try({
      
      xh <- do.call(pdistr, c(list(xm), distpar))
      xh <- PP0 + (1 - PP0) * xh
    },silent=TRUE)
    
    if(class(xm) == "try-error"){
      
      xm[] <- NA
      
    }
  }
  
  ## Transforming to normal values
  res <- qnorm(xm)
  res <- zoo::zoo(res, dates)
  
  # Returning the result for the particular day
  if(!is.null(trgt)){
    
    val <- res[which(trgt == dates)]
    
    
  } else {
    
    val <- res[length(res)]
    
  }
  
  # Retrieving the parameters of the distribution
  if(params){
    
    pos <- which(rownames(f_params) == "P0")
    if(length(pos) > 0)
      f_params <- f_params[-pos,]
    
    val <- c(as.numeric(val), as.numeric(f_params))
  
  }
  
  return(val)
  
  
}







