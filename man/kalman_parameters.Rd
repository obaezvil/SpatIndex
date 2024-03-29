% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPI_SPEI_at_specific_day.R
\name{kalman_parameters}
\alias{kalman_parameters}
\title{Kalman smoothing to daily SPI and SPEI parameters (only works for gamma and log-logistic)}
\usage{
kalman_parameters(params_list, H, ...)
}
\arguments{
\item{params_list}{List object that contain as many 'SparRaster' objects as parameters in the selected distribution.
This list must contain the names of the parameters. Please see the 'calculate_params' function.}

\item{H}{Covariance matrix or array of disturbance terms \epsilon_tϵ of observation equations. See the KFAS package.}

\item{...}{}
}
\value{
The smoothed parameters in a list.
}
\description{
Kalman smoothing to daily SPI and SPEI parameters (only works for gamma and log-logistic)
}
