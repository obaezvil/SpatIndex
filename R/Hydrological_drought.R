################################################################################
################################################################################
##
##  Author: Oscar M. Baez-Villanueva
##  Collaborators:
##
################################################################################
## Objective: calculate hydrological drought indices
################################################################################
##
## Creation date: 2022-08-23
##
################################################################################
################################################################################

#' Wrapper function from the SCI package to calculate the Standardized Streamflow Index (SSI)
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
#'  #'
#' @return Standardised Streamflow Index.
#' @export
#'
#' @examples
ssi <- function(x, 
                scale, 
                ref_start = NULL, 
                ref_end = NULL, 
                distribution){
  
  # Check zoo data
  if(!zoo::is.zoo(x))
    stop("x must be a zoo object!")
  
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
  
  # Extract dates from the x object
  dates <- zoo::index(x)
  
  # Apply the SSI
  idx <- .ssi(x, scale = scale, ref_start = ref_start, ref_end = ref_end, distribution = distribution)

  return(idx)
  
}






