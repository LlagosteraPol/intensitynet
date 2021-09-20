
rm(list = ls())

#Set working directory
setwd("R")

source("netintensityV2_Interface.R")

# Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
load("../Data/Castellon.RData")

# Node coordinates: Georeferenced coordinates from 'castellon' nodes
load("../Data/nodes.RData")

# Event (crime coordinates)
load("../Data/crimes.RData")

#subset of events
crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)


# If not specified, the default graph is: undirected-general
net <- Netintensity(Castellon, nodes$cx, nodes$cy, crim$X, crim$Y)

# gives a list with: list of mean node intensities, list of edge intensities and range of mean node intensities
data <- net$summary() 

# Plot degree histogram
net$plot("hist")

# Plot graph with edge intensities
net$plot("graph")