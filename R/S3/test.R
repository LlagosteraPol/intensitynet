
# to throw error: stop(sprintf("The object must be of class 'intensitynet', current class is %s", class(tt)[1]))


rm(list = ls())

source("S3/main.R")


# Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
load("../Data/Castellon.RData")

# Node coordinates: Georeferenced coordinates from 'castellon' nodes
load("../Data/nodes.RData")

# Event (crime coordinates)
load("../Data/crimes.RData")


#subset of events
crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)


init_graph <- function(obj){
  g <- graph_from_adjacency_matrix(obj$adj_mtx, mode = obj$graph_type)
  
  net_coords <- list(graph = g, node_coords = obj$node_coords)
  class(net_coords) <- "netTools"
  g <- setNetCoords(net_coords)
  
  g # return
}

setNetCoords <- function(obj){
  UseMethod("setNetCoords")
}

#TODO: Declare non-visible at namespace
setNetCoords.netTools = function(obj){
  x_coord_node <- obj$node_coords[1]
  y_coord_node <- obj$node_coords[2]
  
  # TODO: change x and y to long (longitude) and lat (latitude), beware of the coordinate system type (e.g WGS84)
  g <- obj$g %>% set_vertex_attr(name = "xcoord", value = x_coord_node) %>% 
    set_vertex_attr(name = "ycoord", value = y_coord_node)
  g
}

net_setup <- list(adj_mtx = as.matrix(Castellon), node_coords = nodes, events_mtx = crim, graph_type='undirected')
g <- init_graph(net_setup)