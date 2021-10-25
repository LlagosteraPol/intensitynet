require(aoos)
require(contoureR)
require(igraph)
require(intergraph)
require(intervals)
require(ggplot2)
require(ggraph)
require(roxygen2)
require(sna)
require(spatstat)
require(spdep)
require(visNetwork)

source("./intensitynetDir.R", local = TRUE)
source("./intensitynetMix.R", local = TRUE)
source("./intensitynetUnd.R", local = TRUE)
source("./netTools.R", local = TRUE)

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
                    distances_mtx = dist_mtx, 
                    graph_type = graph_type)
  class(net_setup) <- "netTools"
  g <- InitGraph(net_setup)
  
  
  
  intnet <- list(graph = g, events = events_mtx, graph_type = graph_type, distances_mtx = dist_mtx)
  attr(intnet, 'class') <- "intensitynet"
  
  # Select the proper class
  switch(graph_type, 
         'undirected' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetUnd")},
           'directed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetDir")},
              'mixed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetMix")})
  
  intnet # return
}

# -------- Network functions ----------
NodeGeneralCorrelation <- function(obj, dep_type, lag_max, intensity){
  UseMethod("NodeGeneralCorrelation")
}

NodeLocalCorrelation <- function(obj, dep_type = 'moran_i', intensity){
  UseMethod("NodeLocalCorrelation")
}

plot <- function(obj, vertex_intensity='none', edge_intensity='none', xy_axes=TRUE, enable_grid=FALSE, ...){
  UseMethod("plot")
}

gplot <- function(obj, heatmap='none', ...){
  UseMethod("gplot")
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

SetNetworkAttribute <- function(obj, where, name, value){
  UseMethod("SetNetworkAttribute")
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
EdgeIntensity.intensitynet <- function(obj,  node_id1, node_id2, z=5){
  
  if(node_id1 == node_id2){
    stop("The two vertices cannot be the same.")
  }
  
  if(z <= 0){
    print("Warning: 'z' cannot be equal or less than 0, using default.")
    z <- 15
  }
  
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
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
    ep <- c(events_mtx[row, 1], events_mtx[row, 2])
    dist_obj <- list(p1 = node1, p2 = node2, ep = ep)
    class(dist_obj) <- 'netTools'
    d <- PointToLine(dist_obj)
    
    # If the event is at a distance less or equal 'z' from the edge connecting
    # both given points (the road), then is counted as an event of that road
    if(min(node1[1], node2[1]) - z <= ep[1] & ep[1] <= max(node1[1], node2[1]) + z & 
       min(node1[2], node2[2]) - z <= ep[2] & ep[2] <= max(node1[2], node2[2]) + z & 
       d <= z){
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
  
  path_intensity
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


#' Gives general node correlation of the network (choosing from: normal correlation, covariance, 
#' moran-i or geary-g)
#' 
#' @name NodeGeneralCorrelation.intensitynet
#'
#' @param obj intensitynet object
#' @param dep_type the type of dependence statistic to be computed ("correlation", "covariance",
#' "moran", "geary").
#' @param lag_max Maximum geodesic lag at which to compute dependence
#' @param intensity vector containing the intensity values that the heatmaps
#' 
#' @return A vector containing the dependence statistics (ascending from order 0). 
#' 
NodeGeneralCorrelation.intensitynet <- function(obj, dep_type, lag_max, intensity){
  g <- obj$graph
  g_sna <- intergraph::asNetwork(g)
  nacf(g_sna, intensity, type = dep_type, mode = "graph", lag.max = lag_max)
}


#' Gives node local moran-i or geary-g correlations
#' 
#' @name NodeLocalCorrelation.intensitynet
#'
#' @param obj intensitynet object
#' @param dep_type the type of dependence statistic to be computed ('moran_i' or 'geary_g'),
#' default = 'moran_i.
#' @param intensity vector containing the intensity values that the heatmaps
#' 
#' @return An intnet object wich contains a igraph network with the selected correlation added
#' 
NodeLocalCorrelation.intensitynet <- function(obj, dep_type = 'moran_i', intensity){
  g <- obj$graph
  adj_mtx <- as_adj(graph = g)
  adj_listw <- mat2listw(adj_mtx)
  nb <- adj_listw$neighbours
  
  if(dep_type=='geary_g'){
    b_listw <- nb2listw(nb, style="B", zero.policy=TRUE) 
    locg <- localG(x = intensity, listw = b_listw)
    intnet <- SetNetworkAttribute(obj = obj, where = 'vertex', name = "geary_g", value = locg)
    return(list(correlation = locg, intnet = intnet))
  } else{
    w_listw <- nb2listw(nb, style="W", zero.policy=TRUE) 
    locmoran <- localmoran(x = intensity, listw = w_listw, zero.policy=TRUE)
    intnet <- SetNetworkAttribute(obj = obj, where = 'vertex', name = 'moran_i', value = locmoran[, 'Ii'])
    return(list(correlation = locmoran, intnet = intnet))
  } 
}


#' Plot the network and if specified, the moran_i or geary_g heatmap.
#'
#' @name gplot.intensitynet
#'
#' @param obj intensitynet object
#' @param intensity vector containing the intensity values that the heatmaps
#' will use. Default value = NULL
#' @param heatmap local 'moran_i' or 'geary_g'
#' @param ... extra arguments for the class gplot
#' 
gplot.intensitynet  <- function(obj, intensity = NULL, heatmap='none', ...){
  g <- obj$graph
  adj_mtx <- as_adj(graph = g)
  adj_listw <- mat2listw(adj_mtx)
  nb <- adj_listw$neighbours
  
  if(is.null(intensity)){
    intensity <- if(!is.null(vertex_attr(g)$intensity)) vertex_attr(g)$intensity else vertex_attr(g)$intensity_in
  }
  
  node_coords <- data.frame(xcoord = vertex_attr(g)$xcoord, ycoord = vertex_attr(g)$ycoord)
  rownames(node_coords) <- sprintf("V%s",seq(1:nrow(node_coords)))
  
  if(heatmap == 'moran_i'){
    w_listw <- nb2listw(nb, style="W", zero.policy=TRUE) 
    locmoran <- localmoran(x = intensity, listw = w_listw, zero.policy=TRUE, na.action = na.omit)
    
    # Calculate deviations
    node_int_deviation <- intensity - mean(intensity)  
    locmoran_deviation <- locmoran[, 'Ii'] - mean(locmoran[, 'Ii'])
    
    # create a new variable identifying the moran plot quadrant for each observation, dismissing the non-significant ones
    quad_sig <- NA
    significance <- 0.5
    
    # non-significant observations
    quad_sig[(locmoran[, 5] > significance)] <- 0 # "insignificant"  
    # low-low quadrant
    quad_sig[(node_int_deviation < 0 & locmoran_deviation < 0) & (locmoran[, 5] <= significance)] <- 1 # "low-low"
    # low-high quadrant
    quad_sig[(node_int_deviation < 0 & locmoran_deviation > 0) & (locmoran[, 5] <= significance)] <- 2 # "low-high"
    # high-low quadrant
    quad_sig[(node_int_deviation > 0 & locmoran_deviation < 0) & (locmoran[, 5] <= significance)] <- 3 # "high-low"
    # high-high quadrant
    quad_sig[(node_int_deviation > 0 & locmoran_deviation > 0) & (locmoran[, 5] <= significance)] <- 4 # "high-high"
    
    
    data_df <- data.frame(intensity = intensity , 
                          xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord, 
                          heatmap = quad_sig)
  }else if(heatmap=='geary_g'){
    b_listw <- nb2listw(nb, style="B", zero.policy=TRUE) 
    # local net G
    locg_all <- localG(x = intensity, listw = b_listw)
    locg <- unlist(as.list(round(locg_all, 1)))
    
    data_df <- data.frame(intensity = intensity, 
                          xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord, 
                          heatmap = locg)
    
  }else{
    data_df <- data.frame(intensity = intensity, 
                          xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord, 
                          heatmap = NA)
  }
  geoplot_obj <- list(graph=g, data_df = data_df, mode=heatmap)
  class(geoplot_obj) <- "netTools"
  
  GeoreferencedGgplot2(geoplot_obj, ...)
}

#' Set attributes to the network edges or nodes
#'
#' @name SetNetworkAttribute.intensitynet
#'
#' @param obj intensitynet object
#' @param where 'vertex' or 'edge', where to set the attribute
#' @param name name of the attribute
#' @param value vector containing the data for the attribute
#' 
#' @return intensitynet object containing the network with the added attributes
#' 
SetNetworkAttribute.intensitynet <- function(obj, where, name, value){
  g <- obj$graph
  
  if(where == 'edge') g <- g %>% set_edge_attr(name = name, value = value)
  else g <- g %>% set_vertex_attr(name = name, value = value)
  
  intnet <- list(graph = g, events = obj$events, graph_type = obj$graph_type, distances_mtx = obj$distances_mtx)
  attr(intnet, 'class') <- class(obj)
  intnet
}