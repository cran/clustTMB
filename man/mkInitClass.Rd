% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genInit.R
\name{mkInitClass}
\alias{mkInitClass}
\title{Apply classification method dependent on init.method}
\usage{
mkInitClass(n.g, n.i, n.j, control, y)
}
\arguments{
\item{n.g}{Number of clusters}

\item{n.i}{Number of observations}

\item{n.j}{Number of columns}

\item{control}{Classification settings from \code{\link[=init.options]{init.options()}}}

\item{y}{Observations}
}
\value{
classification vector
}
\description{
Apply classification method dependent on init.method
}
\examples{
data("faithful")
mkInitClass(2, nrow(faithful), ncol(faithful), init.options(), faithful)
}
