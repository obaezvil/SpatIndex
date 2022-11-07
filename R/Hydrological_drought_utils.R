################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: util functions for hydrological drought indices
################################################################################
##
## Creation date: 2022-08-23
##
################################################################################
################################################################################

#' Title
#'
#' @param x Zoo object with one or multiple monthly time series.
#' @param scale Integer value that represents the time scale at which the SSI will be computed.
#' @param ref_start optional value that represents the starting point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.
#' @param ref_end Optional value that represents the ending point of the reference period used for computing the index. 
#'  The date should be introduced as '\%Y-\%m'. For example: "1989-02".
#'  The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.
#' @param distribution Optional value indicating the name of the distribution function to be used for computing the SSI 
#'  (one of 'log-Logistic', 'GEV', 'Gamma', and 'PearsonIII').
#'
#' @return
#' @export
#'
#' @examples
.ssi <- function(x,
                 scale,
                 ref_start,
                 ref_end,
                 distribution){
  
  n_col <- ncol(data.frame(x))
  # Setting by default 'pos_ini' and 'pos_fin' as the first and last layer
  pos_ini <- 1
  if(n_col != 1){
    pos_fin <- nrow(x)
  } else {
    pos_fin <- length(x)
  }
  
  
  # Format ref_start if it is not NULL
  if(!is.null(ref_start)){
    
    pos_ini <- grep(ref_start, dates)
    
  }
  
  # Format ref_end if it is not NULL
  if(!is.null(ref_start)){
    
    pos_fin <- grep(ref_end, dates)
    
  }
  
  # Checking distribution
  if(!distribution %in% c("Gamma", "log-Logistic", "GEV", "PearsonIII"))
    stop("The accepted distributions for the SCI package are 'Gamma', 'log-Logistic', 'GEV', and 'PearsonIII'")
  
  # Adjusting distributions for SCI package
  if(distribution == "Gamma")
    distribution <- tolower(distribution)
  
  if(distribution == "log-Logistic")
    distribution <- 'logis'
  
  if(distribution == "GEV"){
    require(evd)
    distribution <- 'gev'
  }
    
  if(distribution == "PearsonIII")
    distribution <- 'pe3'
  
  # Creating an object to store the SSI data
  res <- x
  
  # Defining number of columns
  if(is.null(ncol(x))){
    col <- 1
  } else {
    col <- ncol(x)
  }
  
  # Iterating for every time series in x
  for(i in 1:col){
    
    # Subsetting the data
    val <- as.numeric(x[,i])
    
    # Obtaining parameters for the reference period
    ssi_ref_params <- SCI::fitSCI(val[pos_ini:pos_fin], first.mon = 1, 
                                   distr = distribution, 
                                   time.scale = scale, p0 = FALSE)
    
    # Apply the the parameters of the reference period to all data
    ssi <- SCI::transformSCI(val, first.mon = 1, obj = ssi_ref_params)
    
    # Trasforming NaNs to NAs
    ssi[is.nan(ssi)] <- NA
    
    # Populating the res object with the results
    if(n_col != 1){
      res[,i] <- ssi
    } else {
      res <- ssi
    }
    
  } # end iterative process
  
  return(res)
  
}

#' Utils function to calculate the Snow Water Equivalent Index (SWEI)
#'
#' @param x Numerical vector.
#' @param dates Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_swei' function.
#' @param scale Integer value that represents the time scale at which the ESSMI will be computed.
#' @param missing_ratio Ratio of missing data that is acceptable for the computation. Set to 0.2 by default (20\%).
#' @return Numerical vector with the corresponding to the SWEI
#'
#' @examples
.swei <- function(x, 
                   dates,
                   scale, 
                   missing_ratio){
  
  # Applying the accumulation according to scale
  x_zoo <- zoo::zoo(x, dates)
  x_zoo <- zoo::rollapply(x_zoo, scale, fill = NA, FUN = sum, align = "right")
  
  # Setting a threshold of 1mm 
  x_zoo[x_zoo < 1] <- 0
  
  # Obtaining the resulting unique steps
  steps   <- substr(zoo::index(x_zoo), 6, 10)
  u_steps <- unique(steps)
  
  # Calculating NAs in period
  nas <- length(which(is.na(x_zoo)))
  
  if(length(nas) < 1)
    nas <- 0
  
  missing <- nas / length(x_zoo)
  
  # Calculating unique values
  uniq_vals <- length(unique(x_zoo)) / length(x_zoo)
  
  # Conditional: if all values are NAs, return NAs (masked regions or oceans)
  if(missing_ratio < missing){
    
    swei <- rep(NA, length(dates))
    
  } else {
    
    # Creating an object to store the resulting data
    swei <- c()
    
    #######
    ####### Loop to evaluate each one of the time steps
    #######
    
    for(i in 1:length(u_steps)){
      
      # Extracting the data according to step 'i'
      pos_step   <- which(steps %in% u_steps[i])
      x_step     <- x_zoo[pos_step]
      
      # Dealing with zero values (see Huning et al 2020: https://www.pnas.org/doi/epdf/10.1073/pnas.1915921117)
      pos_zero <- which(x_step == 0)
      
      if(length(pos_zero) > 0){
        suppressWarnings({
          min_val          <- x_step[-pos_zero] - min(x_step[-pos_zero], na.rm = TRUE) / 100
          min_val          <- min(min_val, na.rm = TRUE)
      })
        x_step[pos_zero] <- min_val
      }
      
      # Storing the position of the NA values (if any)
      pos_nas <- which(is.na(x_step))
      
      # Removing NA values from vector
      if(length(pos_nas) > 0)
        x_step <- x_step[-pos_nas]
      
      # Ordering the vector in increasing form
      x_order <- order(x_step)
      x_sort  <- round(x_step[x_order], 2)
      
      # Calculating the unique values of x_sort
      x_sort_unique <- unique(x_sort)
      
      # Applying the empirical Gringorten plotting position
      p_ami <- (1:length(x_sort_unique)-0.44) / (length(x_sort_unique) + 0.12)
      
      # Matching equal values to x_sort
      pos_equal_vals <- match(x_sort, x_sort_unique)
      p_ami          <- p_ami[pos_equal_vals]
      
      # Applying the inverse standard normal distribution
      swei_ordered <- qnorm(p_ami)
      
      # Ordering the values according to the original data
      res <- swei_ordered[order(x_order)]
      
      # Including the NA values (if any)
      if(length(pos_nas) > 0)
        res <- R.utils::insert(res, pos_nas, NA)
      
      # Obtaining all values for that time step
      pos_dates   <- grep(u_steps[i], dates)

      
      # Storing the standardized values
      res[which(is.infinite(res))] <- NA
      res[which(is.nan(res))]      <- NA
      
      # Setting the values to NA if all values are zeros
      if("min_val" %in% ls()){

        if(length(min_val) < 1 | is.infinite(min_val))
          res <- rep(NA, length(res))

      }
      
      
      suppressWarnings(swei[pos_dates] <- res)
      
    } # End for
    
  } # End if else missing data
  
  
  return(swei)
  
}


