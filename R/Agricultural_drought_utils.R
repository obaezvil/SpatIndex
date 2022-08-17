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
## Creation date: 2022-08-11
##
################################################################################
################################################################################


#' Utils function to calculate the Vegetation Condition Index
#'
#' @param x Numerical vector.
#'
#' @return Numerical vector with the corresponding to the VCI
#' @export
#'
#' @examples
.vci <- function(x, ...){
  
  # Applying VCI formula
  vci <- ( x - min(x, ...) ) / ( max(x, ...) - min(x, ...) ) * 100
  
  # Replace with NAs over areas with full NA vectors
  if(length(vci) != length(x))
    vci <- rep(NA, length(x))
  
  return(vci)
  
}

#' Utils function to calculate the Temperature Condition Index
#'
#' @param x Numerical vector.
#'
#' @return Numerical vector with the corresponding to the TCI
#' @export
#'
#' @examples
.tci <- function(x, ...){
  
  # Applying VCI formula
  tci <- ( max(x, ...) - x ) / ( max(x, ...) - min(x, ...) ) * 100
  
  # Replace with NAs over areas with full NA vectors
  if(length(tci) != length(x))
    tci <- rep(NA, length(x))
  
  return(tci)
  
}




