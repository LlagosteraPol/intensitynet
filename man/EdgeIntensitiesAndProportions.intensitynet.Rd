% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intensitynet.R
\name{EdgeIntensitiesAndProportions.intensitynet}
\alias{EdgeIntensitiesAndProportions.intensitynet}
\title{Calculate all the edge intensities of the graph.}
\usage{
\method{EdgeIntensitiesAndProportions}{intensitynet}(obj)
}
\arguments{
\item{obj}{intensitynet object}
}
\value{
intensitynet class object where the graph contains all the edge intensities as an attribute
}
\description{
Calculate all the edge intensities of the graph. It's more fast than using iteratively the 
function EdgeIntensity for all edges.
}
