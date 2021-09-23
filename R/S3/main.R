require(spatstat)
require(igraph)
require(intervals)
require(aoos)
require(roxygen2)

source("S3/netTools.R")
source("S3/intensitynetUnd.R")

intensitynet <- function(adj_mtx, node_coords, events_mtx, graph_type = 'undirected'){
  
  if(class(adj_mtx) == "data.frame"){
    adj_mtx <- as.matrix(adj_mtx)
  }
  
  net_setup <- list(adj_mtx = adj_mtx, node_coords = node_coords, events_mtx = events_mtx, graph_type = graph_type)
  class(net_setup) <- "netTools"
  g <- init_graph(net_setup)
  
  #dist_mtx <- calculateDistancesMtx(node_coords)
  
  intnet <- list(graph = g, events = events_mtx, graph_type = graph_type, distances = NULL, intensities = NULL)
  attr(intnet, 'class') <- "intensitynet"
  
  # Select the proper class
  switch(graph_type, 
         'undirected' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetUnd")},
           'directed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetDir")},
              'mixed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetMix")})
  
  intnet # return
}

# -------- Network functions ----------
getDistancesMtx <- function(obj){
  UseMethod("getDistancesMtx")
}

getDistancesMtx.intensitynet <- function(obj){
  if(is.null(obj$distances)){
    calculateDistancesMtx(obj)
  }else{
    obj$distances
  }
}

# -------- Intensity functions ----------
pathIntensity <- function(obj){
  UseMethod("pathIntensity")
}

shortestPathIntensity <- function(obj){
  UseMethod("shortestPathIntensity")
}

nodeIntensity <- function(obj){
  UseMethod("nodeIntensity")
}

meanNodeIntensity <- function(obj){
  UseMethod("meanNodeIntensity")
}

edgeIntensity <- function(obj){
  UseMethod("edgeIntensity")
}

calculateAllIntensities <- function(obj){
  UseMethod("calculateAllIntensities")
}


