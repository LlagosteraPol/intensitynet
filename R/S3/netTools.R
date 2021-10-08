#TODO: Declare all these functions as non-visible at namespace


InitGraph <- function(obj){
  UseMethod("InitGraph")
}

#' Creates an igraph network with the given data
#'
#' @name InitGraph.netTools 
#'
#' @param obj netTools object
#' 
#' @return igraph network
#' 
InitGraph.netTools <- function(obj){
  weighted_mtx = obj$adjacency_mtx * obj$distances
  if(obj$graph_type == 'undirected') g <- graph_from_adjacency_matrix(weighted_mtx, mode = obj$graph_type, weighted=TRUE)
  else  g <- graph_from_adjacency_matrix(weighted_mtx, mode = 'directed', weighted=TRUE)
  
  net_coords <- list(graph = g, node_coords = obj$node_coords)
  class(net_coords) <- "netTools"
  g <- SetNetCoords(net_coords)
  
  g # return
}

SetNetCoords <- function(obj){
  UseMethod("SetNetCoords")
}


#' Set igraph network node coordinates as its attributes
#'
#' @name InitGraph.netTools 
#'
#' @param obj netTools object
#' 
#' @return igraph network with the given coordinates as the attributes of the nodes
#' 
SetNetCoords.netTools = function(obj){
  x_coord_node <- obj$node_coords[, 1]
  y_coord_node <- obj$node_coords[, 2]
  
  # TODO: change x and y to long (longitude) and lat (latitude), beware of the coordinate system type (e.g WGS84)
  g <- obj$g %>% set_vertex_attr(name = "xcoord", value = x_coord_node) %>% 
    set_vertex_attr(name = "ycoord", value = y_coord_node)
  g
}

CalculateDistancesMtx <- function(obj){
  UseMethod("CalculateDistancesMtx")
}


#' Calculates the distances between all pair of nodes from the given network
#'
#' @name CalculateDistancesMtx.netTools 
#'
#' @param obj netTools object
#' 
#' @return distances matrix
#' 
CalculateDistancesMtx.netTools <- function(obj){
  x_coord_node <- obj$node_coords[, 1]
  y_coord_node <- obj$node_coords[, 2] 
  
  distances_mtx <- pairdist(ppp(x_coord_node,
                                y_coord_node,
                                xrange=c(min(as.numeric(x_coord_node)), max(as.numeric(x_coord_node))),
                                yrange=c(min(as.numeric(y_coord_node)), max(as.numeric(y_coord_node)))))
  rownames(distances_mtx) <- colnames(distances_mtx) <- sprintf("V%s",seq(1:ncol(distances_mtx)))
  distances_mtx
}

SetEdgeIntensity <- function(obj){
  UseMethod("SetEdgeIntensity")
}


#' Sets the given intensites as an edge attribute to the given igraph network
#'
#' @name SetEdgeIntensity.netTools 
#'
#' @param obj netTools object
#' 
#' @return igraph network with the given intensities as an attributes of the edges
#' 
#TODO: Declare non-visible at namespace
SetEdgeIntensity.netTools <- function(obj){
  g <- obj$graph
  node_id1 <- obj$node_id1
  node_id2 <- obj$node_id2
  value <- obj$value
  
  edge_id <- get.edge.ids(g, c(node_id1, node_id2))
  g <- g %>% set_edge_attr(name = "intensity", index = edge_id, value = value)
  g
}

SetNodeIntensity <- function(obj){
  UseMethod("SetNodeIntensity")
}


#' Sets the given intensites as a node attribute to the given igraph network
#'
#' @name SetNodeIntensity.netTools 
#'
#' @param obj netTools object
#' 
#' @return igraph network with the given intensities as an attributes of the nodes
SetNodeIntensity.netTools = function(obj){
  g <- obj$graph
  node_id <- obj$node_id
  intensity <- obj$intensity
  
  g <- g %>% set_vertex_attr(name = "intensity", index = node_id, value = intensity)
  g
}


ShortestDistance <- function(obj){
  UseMethod("ShortestDistance")
}


#' Calculates the shortest distance path between two nodes 
#'
#' @name ShortestDistance.netTools 
#'
#' @param obj netTools object
#' 
#' @return distance of the path and the nodes of the path
#' 
ShortestDistance.netTools = function(obj){
  g <- obj$graph
  node_id1 <- obj$node_id1
  node_id2 <- obj$node_id2
  distances_mtx <- obj$distances_mtx
  
  weighted_path <- unlist(get.shortest.paths(g, node_id1, node_id2)$vpath)
  if(!is.null(distances_mtx)){
    weight_sum <- sum(E(g, path = unlist(weighted_path))$weight)
  }
  else{
    weight_sum <- length(weighted_path)
  }
  list(weight = weight_sum, path = weighted_path)  
}

GeoreferencedPlot <- function(obj){
  UseMethod("GeoreferencedPlot")
}


#' Plot the given network using its node coordinates
#'
#' @name GeoreferencedPlot.netTools 
#'
#' @param obj netTools object
#' 
GeoreferencedPlot.netTools = function(obj){
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  
  if(!is.null(distances_mtx)){
    norm_coords = layout.norm(matrix(cbind(vertex_attr(g)$xcoord, vertex_attr(g)$ycoord), ncol=2))
    plot(g, 
         layout = norm_coords, 
         vertex.label=NA, 
         vertex.size=2, 
         window=TRUE, 
         axes=TRUE, 
         edge.label = edge_attr(g)$intensity, 
         edge.label.cex = 0.5)
  }
  else{
    plot(g, 
         vertex.label=NA, 
         vertex.size=2,
         vertex.size2=2)
  }
}


PointToLine <- function(obj){
  UseMethod("PointToLine")
}


#' Return the perpendicular distance between an event and the edge between two nodes.
#'
#' @name PointToLine.netTools  
#'
#' @param obj netTools object -> list(p1:c(coordx, coordy), p2:c(coordx, coordy), e:c(coordx, coordy))
#' 
#' @return the perpendicular distance
#' 
PointToLine.netTools <- function(obj){
  p1 <- obj$p1
  p2 <- obj$p2
  ep <- obj$ep
  
  v1 <- p1 - p2
  v2 <- p2 - ep
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  
  d
}