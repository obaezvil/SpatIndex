% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPI_SPEI_at_specific_day_utils.R
\name{.get_params.sci}
\alias{.get_params.sci}
\title{Utils function to calculate the parameters of the selected distribution function for the SPI or SPEI using a vector. This function is based on the code presented in the SCI package
 obtained from  https://github.com/cran/SCI/blob/master/R/sci.r}
\usage{
.get_params.sci(x, trgt, dates, ref_start, ref_end, distribution, fit, params)
}
\arguments{
\item{x}{Numerical vector.}

\item{trgt}{The day for which the function will be computed. The default is NULL, indicating that the last day of the data
will be selected.}

\item{dates}{Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_spei' function.}

\item{ref_start}{optional value that represents the starting point of the reference period used for computing the index. 
The date should be introduced as '\%Y-\%m'. For example: "1989-02".
The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.}

\item{ref_end}{Optional value that represents the ending point of the reference period used for computing the index. 
The date should be introduced as '\%Y-\%m'. For example: "1989-02".
The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.}

\item{distribution}{Optional value indicating the name of the distribution function to be used for computing the SPI 
(one of 'log-Logistic', 'Gamma' and 'PearsonIII'). Defaults to 'log-Logistic' for SPEI.}

\item{fit}{Optional value indicating the name of the method used for computing the distribution function parameters 
(one of 'ub-pwm', 'pp-pwm' and 'max-lik'). Defaults to 'ub-pwm'.}

\item{params}{Should the parameters of the selected distributions be returned? Set to FALSE as the default.}
}
\description{
Utils function to calculate the parameters of the selected distribution function for the SPI or SPEI using a vector. This function is based on the code presented in the SCI package
 obtained from  https://github.com/cran/SCI/blob/master/R/sci.r
}
