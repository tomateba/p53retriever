% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/p53density.R
\name{p53density}
\alias{p53density}
\title{Plot the density chart of responsive elements on a DNA sequence}
\usage{
p53density(p53.table, start.at, stop.at)
}
\arguments{
\item{p53.table}{A dataframe containing the responsive elements, as the one returned from the p53track function.}

\item{start.at}{A numeric value indicating the start coordinate of the plot}

\item{stop.at}{A numeric value indicating the stop coordinate of the plot}
}
\value{
an object of class \code{ggplot}
}
\description{
Plot the density chart of responsive elements on a DNA sequence
}
\examples{
data(RE_collection)
p53density(RE_collection,-10000,10000)
}

