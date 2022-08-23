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
  
  # Setting by default 'pos_ini' and 'pos_fin' as the first and last layer
  pos_ini <- 1
  pos_fin <- nrow(x)
  
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
  
  if(distribution == "GEV")
    distribution <- 'gev'
  
  if(distribution == "PearsonIII")
    distribution <- 'pe3'
  
  # Creating an object to store the SSI data
  res <- x
  
  # Iterating for every time series in x
  for(i in 1:ncol(x)){
    
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
    res[,i] <- ssi
    
  } # end iterative process
  
  return(res)
  
}



