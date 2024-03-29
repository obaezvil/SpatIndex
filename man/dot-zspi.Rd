% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Combined_drought_indices_utils.R
\name{.zspi}
\alias{.zspi}
\title{Utils function to calculate the boolean SPI (zSPI) to be used in the CDI}
\usage{
.zspi(SPI1_data, SPI3_data, threshold_spi1, threshold_spi3)
}
\arguments{
\item{SPI1_data}{'SpatRaster' object that contains the SPI-1 data that will be used to calculate the zSPI. 
This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
of the terra package.}

\item{SPI3_data}{'SpatRaster' object that contains the SPI-3 data that will be used to calculate the zSPI. 
This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
of the terra package.}

\item{threshold_spi1}{Threshold to construct the Boolean layer for SPI-1. Set to -2 by default.}

\item{threshold_spi3}{Threshold to construct the Boolean layer for SPI-3. Set to -1 by default.}
}
\value{
Returns the zSPI
}
\description{
Utils function to calculate the boolean SPI (zSPI) to be used in the CDI
}
