require(aoos)
require(contoureR)
require(igraph)
require(intergraph)
require(intervals)
require(roxygen2)
require(sna)
require(spatstat)
require(spdep)

source("S3/intensitynetDir.R")
source("S3/intensitynetMix.R")
source("S3/intensitynetUnd.R")
source("S3/netTools.R")

#' Constructor of the class intensitynet
#'
#' @name intensitynet
#'
#' @param adjacency_mtx Network adjacency matrix
#' @param node_coords Nodes latitude and longitude matrix
#' @param events_mtx Events latitude and longitude matrix
#' @param graph_type Network type: 'undirected' (default), 'directed' or 'mixed' 
#' 
#' @return intensitynet object containing: graph=<igraph>, events = <matrix>, graph_type = c('directed', 'undirected', 'mixed'), 
#' distances = <matrix>
#' 
intensitynet <- function(adjacency_mtx, node_coords, events_mtx, graph_type = 'undirected'){
  
  if(class(adjacency_mtx) == "data.frame"){
    adjacency_mtx <- as.matrix(adjacency_mtx)
  }
  
  if(class(node_coords) == "data.frame"){
    node_coords <- as.matrix(node_coords)
  }
  
  node_coords_obj <- list(node_coords = node_coords)
  class(node_coords_obj) <- "netTools"
  dist_mtx <- CalculateDistancesMtx(node_coords_obj)
  
  net_setup <- list(adjacency_mtx = adjacency_mtx, 
                    node_coords = node_coords, 
                    events = events_mtx, 
                    distances = dist_mtx, 
                    graph_type = graph_type)
  class(net_setup) <- "netTools"
  g <- InitGraph(net_setup)
  
  
  
  intnet <- list(graph = g, events = events_mtx, graph_type = graph_type, distances = dist_mtx)
  attr(intnet, 'class') <- "intensitynet"
  
  # Select the proper class
  switch(graph_type, 
         'undirected' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetUnd")},
           'directed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetDir")},
              'mixed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetMix")})
  
  intnet # return
}

#intnet_all <- CalculateEventIntensities(intnet)
plot <- function(plot, ...){
  UseMethod("plot")
}

#' Plot intensitynet object
#'
#' @name plot.intensitynet 
#'
#' @param obj intensitynet object
#' 
plot.intensitynet <- function(obj, ...){
  geoplot_obj <- list(graph=obj$graph, distances_mtx = obj$distances)
  class(geoplot_obj) <- "netTools"
  GeoreferencedPlot(geoplot_obj)
}

# -------- Network functions ----------
EventCorrelation <- function(obj, edge_id1, edge_id2){
  UseMethod("EventCorrelation")
}

NodeLocalCorrelation <- function(obj, mode){
  UseMethod("NodeLocalCorrelation")
}

# -------- Intensity functions ----------
PathIntensity <- function(obj, path_nodes){
  UseMethod("PathIntensity")
}

ShortestPathIntensity <- function(obj,  node_id1, node_id2, weighted = FALSE){
  UseMethod("ShortestPathIntensity")
}

CalculateEventIntensities <- function(obj){
  UseMethod("CalculateEventIntensities")
}

MeanNodeIntensity <- function(obj, node_id){
  UseMethod("MeanNodeIntensity")
}

EdgeIntensity <- function(obj, node_id1, node_id2, z){
  UseMethod("EdgeIntensity")
}

#' If not calculated, calculates the intesnity of the edge with nodes; node_id1, node_id2. 
#' If the edge already contains an intensity, gives it directly.
#'
#' @name EdgeIntensity.intensitynetUnd
#' 
#' @param node_id1 First node ID of the edge
#' @param node_id2 Second node ID of the edge
#' 
#' @return edge_intensity - Intensity of the edge
#'
#TODO: Set function as non-visible
EdgeIntensity.intensitynet= function(obj,  node_id1, node_id2, z=50){
  
  if(node_id1 == node_id2){
    stop("Both vertices cannot be the same.")
  }
  
  if(z <= 0){
    print("Warning: 'z' cannot be equal or less than 0, using default.")
    z <- 50
  }
  
  g <- obj$graph
  distances_mtx <- obj$distances
  events_mtx <- obj$events
  
  # Note that the igraph library already handle the error when one of the node id's 
  # are not part of the graph. Also gives the proper information about it.
  edge_id <- get.edge.ids(g, c(node_id1, node_id2))
  
  # If the intensity of this edge was previously calculated, then return it
  if(edge_id != 0 & !is.null(edge_attr(g, "intensity", index=edge_id))){
    if(length(is.na(vertex_attr(g, "intensity", edge_id)))==0){
      return(edge_attr(g, 'intensity', index=edge_id))
    }
  }
  
  # Distance between the node and its neighbor
  res <- tryCatch(
    {
      abs(distances_mtx[node_id1, node_id2]) # Distance between the node and its neighbor 
    },
    # If the nodes are not part of the graph, give the proper information of the error.
    error=function(cond) {
      neighbors_list <- neighbors(g, node_id1)
      if(! V(g)[node_id2] %in% neighbors_list){
        message("Second vertice (node_id2) it's not a neighbor of first vertice (node_id1)")
      }else{
        message(cond)
      }
    }
  )    
  
  node1 <- c(vertex_attr(g, "xcoord", node_id1), vertex_attr(g, "ycoord", node_id1))
  node2 <- c(vertex_attr(g, "xcoord", node_id2), vertex_attr(g, "ycoord", node_id2))
  
  indicator <- 0
  # Counting events
  for(row in 1:nrow(events_mtx)) {
    dist_obj <- list(p1= node1, p2= node2, ep=c(events_mtx[row, 1], events_mtx[row, 2]))
    class(dist_obj) <- 'netTools'
    d <- PointToLine(dist_obj)
    
    # If the event is at a perpendicular distance less or equal 'z' from the line connecting
    # both given points (the road), then is counted as an event of that road
    if(d <= z){
      indicator <- indicator + 1
    } 
  }
  edge_intensity <- indicator/res
  
  edge_intensity
}


#' Calculates the intensity of the given path
#'
#' @name PathIntensity.intensitynet
#'
#' @param obj intensitynet object
#' @param path_nodes vector containing the node ID's of the path
#' 
#' @return intensity of the path
#' 
PathIntensity.intensitynet <- function(obj, path_nodes){
  edge_counts <- list()
  path_intensity <- 0
  
  prev <- NULL
  for(node_id in path_nodes){
    if(is.null(prev)){
      prev <- node_id
      next
    }
    
    path_intensity <- path_intensity + Reduce('+', EdgeIntensity(obj, prev, node_id))
    
    prev <- node_id
  }
  # Divide the intensity of the edges by their number (In a path -> N edges = N vertices - 1)
  path_intensity <- path_intensity / (length(path_nodes) - 1)
}


#' Calculates the shortest path between two vertices and calculates its intensity
#'
#' @name ShortestPathIntensity.intensitynet
#'
#' @param obj intensitynet object
#' @param node_id1 starting node
#' @param node_id2 ending node
#' @param weighted TRUE or FALSE (default), tell if the distances must be taken into account 
#' 
#' @return intensity of the path the shortest path and the path
#' 
ShortestPathIntensity.intensitynet <- function(obj,  node_id1, node_id2, weighted = FALSE){
  g <- obj$graph
  
  if(weighted){
    path <- ShortestDistance(node_id1, node_id2)$path
  }else{
    path <- unlist(get.shortest.paths(g, node_id1, node_id2)$vpath)
  }
  
  return(list(intensity = PathIntensity(path), path = path))
}
