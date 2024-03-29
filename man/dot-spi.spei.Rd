% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Meteorological_drought_utils.R
\name{.spi.spei}
\alias{.spi.spei}
\title{Utils function to calculate the SPI (SPEI package)}
\usage{
.spi.spei(x, dates, scale, ref_start, ref_end, distribution, fit, ...)
}
\arguments{
\item{x}{Numerical vector.}

\item{dates}{Vector of dates that is extracted from the 'SpatRaster' in the 'spatial_spi' function.}

\item{scale}{Integer value that represents the time scale at which the SPI will be computed.}

\item{ref_start}{optional value that represents the starting point of the reference period used for computing the index. 
The date should be introduced as '\%Y-\%m'. For example: "1989-02".
The default is NULL, which indicates that the first layer in the 'SpatRaster' will be used as starting point.}

\item{ref_end}{Optional value that represents the ending point of the reference period used for computing the index. 
The date should be introduced as '\%Y-\%m'. For example: "1989-02".
The default is NULL, which indicates that the last layer in the 'SpatRaster' will be used as ending point.}

\item{distribution}{Optional value indicating the name of the distribution function to be used for computing the SPI 
(one of 'log-Logistic', 'Gamma' and 'PearsonIII'). Defaults to 'Gamma' for SPI.}

\item{fit}{Optional value indicating the name of the method used for computing the distribution function parameters 
(one of 'ub-pwm', 'pp-pwm' and 'max-lik'). Defaults to 'ub-pwm'.}

\item{...}{Additional variables that can be used for the 'spi' function of the SPEI package.}
}
\description{
Utils function to calculate the SPI (SPEI package)
}
