% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Climdex_utils.R
\name{.idx.rnnptot}
\alias{.idx.rnnptot}
\title{Rnnp (threshold should be applied: for R95pTOT thres = 0.95; for R99pTOT thres = 0.99)}
\usage{
.idx.rnnptot(rst, thres)
}
\arguments{
\item{rst}{'SpatRaster' object with daily precipitation values.}

\item{thres}{Threshold value that corresponds to nn.}
}
\value{
Contribution to total precipitation from very wet days (0.95) and extremely wet days (0.99)
}
\description{
Rnnp (threshold should be applied: for R95pTOT thres = 0.95; for R99pTOT thres = 0.99)
}
