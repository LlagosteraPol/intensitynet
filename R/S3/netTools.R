#TODO: Declare all these functions as non-visible at namespace


init_graph <- function(obj){
  UseMethod("init_graph")
}

init_graph.netTools <- function(obj){
  g <- graph_from_adjacency_matrix(obj$adj_mtx, mode = obj$graph_type)
  
  net_coords <- list(graph = g, node_coords = obj$node_coords)
  class(net_coords) <- "netTools"
  g <- setNetCoords(net_coords)
  
  g # return
}

setNetCoords <- function(obj){
  UseMethod("setNetCoords")
}

setNetCoords.netTools = function(obj){
  x_coord_node <- obj$node_coords[1]
  y_coord_node <- obj$node_coords[2]
  
  # TODO: change x and y to long (longitude) and lat (latitude), beware of the coordinate system type (e.g WGS84)
  g <- obj$g %>% set_vertex_attr(name = "xcoord", value = x_coord_node) %>% 
    set_vertex_attr(name = "ycoord", value = y_coord_node)
  g
}

calculateDistancesMtx <- function(obj){
  UseMethod("calculateDistancesMtx")
}

calculateDistancesMtx.netTools <- function(node_coords){
  x_coord_node <- node_coords$cx # TODO: get first column (not by name)
  y_coord_node <- node_coords$cy # TODO: get second column (not by name)
  
  distances_mtx <- pairdist(ppp(x_coord_node,
                                y_coord_node,
                                xrange=c(min(as.numeric(x_coord_node)), max(as.numeric(x_coord_node))),
                                yrange=c(min(as.numeric(y_coord_node)), max(as.numeric(y_coord_node)))))
  distances_mtx
}

setEdgeIntensity <- function(obj){
  UseMethod("setEdgeIntensity")
}

#TODO: Declare non-visible at namespace
setEdgeIntensity.netTools <- function(g, node_id1, node_id2, value){
  edge_id <- get.edge.ids(g, c(node_id1, node_id2))
  g <- g %>% set_edge_attr(name = "intensity", index = edge_id, value = value)
  g
}