% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/intensitynet.R
\name{SetNetworkAttribute.intensitynet}
\alias{SetNetworkAttribute.intensitynet}
\title{Set attributes to the network edges or nodes}
\usage{
\method{SetNetworkAttribute}{intensitynet}(obj, where, name, value)
}
\arguments{
\item{obj}{intensitynet object}

\item{where}{'vertex' or 'edge', where to set the attribute}

\item{name}{name of the attribute}

\item{value}{vector containing the data for the attribute}
}
\value{
intensitynet object containing the network with the added attributes
}
\description{
Set attributes to the network edges or nodes
}
