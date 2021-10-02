require(spatstat)
require(igraph)
require(intervals)
require(aoos)
require(roxygen2)
library(contoureR)

source("S3/netTools.R")
source("S3/intensitynetUnd.R")

intensitynet <- function(adjacency_mtx, node_coords, events_mtx, graph_type = 'undirected'){
  
  if(class(adjacency_mtx) == "data.frame"){
    adjacency_mtx <- as.matrix(adjacency_mtx)
  }
  
  if(class(node_coords) == "data.frame"){
    node_coords <- as.matrix(node_coords)
  }
  
  node_coords_obj <- list(node_coords = node_coords)
  class(node_coords_obj) <- "netTools"
  dist_mtx <- calculateDistancesMtx(node_coords_obj)
  
  net_setup <- list(adjacency_mtx = adjacency_mtx, 
                    node_coords = node_coords, 
                    events = events_mtx, 
                    distances = dist_mtx, 
                    graph_type = graph_type)
  class(net_setup) <- "netTools"
  g <- init_graph(net_setup)
  
  
  
  intnet <- list(graph = g, events = events_mtx, graph_type = graph_type, distances = dist_mtx)
  attr(intnet, 'class') <- "intensitynet"
  
  # Select the proper class
  switch(graph_type, 
         'undirected' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetUnd")},
           'directed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetDir")},
              'mixed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetMix")})
  
  intnet # return
}

#intnet_all <- calculateEventIntensities(intnet)
plot <- function(plot, ...){
  UseMethod("plot")
}

plot.intensitynet <- function(obj, ...){
  geoplot_obj <- list(graph=obj$graph, distances_mtx = obj$distances)
  class(geoplot_obj) <- "netTools"
  georeferencedPlot(geoplot_obj)
}

# -------- Network functions ----------


# -------- Intensity functions ----------
pathIntensity <- function(obj, path_nodes){
  UseMethod("pathIntensity")
}

shortestPathIntensity <- function(obj,  node_id1, node_id2, weighted = FALSE){
  UseMethod("shortestPathIntensity")
}

calculateEventIntensities <- function(obj){
  UseMethod("calculateEventIntensities")
}

meanNodeIntensity.intensitynetUnd <- function(obj, node_id){
  UseMethod("meanNodeIntensity.intensitynetUnd")
}

edgeIntensity <- function(obj, node_id1, node_id2, z){
  UseMethod("edgeIntensity")
}

#' If not calculated, calculates the intesnity of the edge with nodes; node_id1, node_id2 and
#' input into the edge attribute of the graph. If the edge already contains an intensity,
#' gives it directly.
#'
#' @name edgeIntensity.intensitynetUnd
#' 
#' @param node_id1 First node ID of the edge
#' @param node_id2 Second node ID of the edge
#' 
#' @return edge_intensity - Intensity of the edge
#'
#TODO: Set function as non-visible
edgeIntensity.intensitynet= function(obj,  node_id1, node_id2, z=50){
  
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
    d <- perpenidularDistance(dist_obj)
    
    # If the event is at a perpendicular distance less or equal 'z' from the line connecting
    # both given points (the road), then is counted as an event of that road
    if(d <= z){
      indicator <- indicator + 1
    } 
  }
  edge_intensity <- indicator/res
  
  edge_intensity
}

pathIntensity.intensitynet <- function(obj, path_nodes){
  edge_counts <- list()
  path_intensity <- 0
  
  prev <- NULL
  for(node_id in path_nodes){
    if(is.null(prev)){
      prev <- node_id
      next
    }
    
    path_intensity <- path_intensity + Reduce('+', edgeIntensity(obj, prev, node_id))
    
    prev <- node_id
  }
  # Divide the intensity of the edges by their number (In a path -> N edges = N vertices - 1)
  path_intensity <- path_intensity / (length(path_nodes) - 1)
}

shortestPathIntensity.intensitynet <- function(obj,  node_id1, node_id2, weighted = FALSE){
  g <- obj$graph
  
  if(weighted){
    path <- shortestDistance(node_id1, node_id2)$path
  }
  else{
    path <- unlist(get.shortest.paths(g, node_id1, node_id2)$vpath)
  }
  
  return(list(intensity = pathIntensity(path), path = path))
}


