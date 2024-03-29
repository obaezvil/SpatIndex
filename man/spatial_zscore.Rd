% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Agricultural_drought.R
\name{spatial_zscore}
\alias{spatial_zscore}
\title{Anomaly Index (zscore)}
\usage{
spatial_zscore(Var_data)
}
\arguments{
\item{Var_data}{'SpatRaster' object that contains spatially-distributed data that will be used to calculate the zscore. 
This 'SpatRaster' must include the time that corresponds to the dates of the respective layers. They can be set with the function time
of the terra package.}
}
\value{
Spatially-distributed zfAPAR values.
}
\description{
Anomaly Index (zscore)
}
