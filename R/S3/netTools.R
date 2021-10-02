#TODO: Declare all these functions as non-visible at namespace


init_graph <- function(obj){
  UseMethod("init_graph")
}

init_graph.netTools <- function(obj){
  weighted_mtx = obj$adjacency_mtx * obj$distances
  g <- graph_from_adjacency_matrix(weighted_mtx, mode = obj$graph_type, weighted=TRUE)
  
  net_coords <- list(graph = g, node_coords = obj$node_coords)
  class(net_coords) <- "netTools"
  g <- setNetCoords(net_coords)
  
  g # return
}

setNetCoords <- function(obj){
  UseMethod("setNetCoords")
}

setNetCoords.netTools = function(obj){
  x_coord_node <- obj$node_coords[, 1]
  y_coord_node <- obj$node_coords[, 2]
  
  # TODO: change x and y to long (longitude) and lat (latitude), beware of the coordinate system type (e.g WGS84)
  g <- obj$g %>% set_vertex_attr(name = "xcoord", value = x_coord_node) %>% 
    set_vertex_attr(name = "ycoord", value = y_coord_node)
  g
}

calculateDistancesMtx <- function(obj){
  UseMethod("calculateDistancesMtx")
}

calculateDistancesMtx.netTools <- function(obj){
  x_coord_node <- obj$node_coords[, 1]
  y_coord_node <- obj$node_coords[, 2] 
  
  distances_mtx <- pairdist(ppp(x_coord_node,
                                y_coord_node,
                                xrange=c(min(as.numeric(x_coord_node)), max(as.numeric(x_coord_node))),
                                yrange=c(min(as.numeric(y_coord_node)), max(as.numeric(y_coord_node)))))
  rownames(distances_mtx) <- colnames(distances_mtx) <- sprintf("V%s",seq(1:ncol(distances_mtx)))
  distances_mtx
}

setEdgeIntensity <- function(obj){
  UseMethod("setEdgeIntensity")
}

#TODO: Declare non-visible at namespace
setEdgeIntensity.netTools <- function(obj){
  g <- obj$graph
  node_id1 <- obj$node_id1
  node_id2 <- obj$node_id2
  value <- obj$value
  
  edge_id <- get.edge.ids(g, c(node_id1, node_id2))
  g <- g %>% set_edge_attr(name = "intensity", index = edge_id, value = value)
  g
}

setNodeIntensity <- function(obj){
  UseMethod("setNodeIntensity")
}

setNodeIntensity.netTools = function(obj){
  g <- obj$graph
  node_id <- obj$node_id
  intensity <- obj$intensity
  
  g <- g %>% set_vertex_attr(name = "intensity", index = node_id, value = intensity)
  g
}


shortestDistance <- function(obj){
  UseMethod("shortestDistance")
}

shortestDistance.netTools = function(obj){
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
  list(path = weighted_path, weight = weight_sum)  
}

georeferencedPlot <- function(obj){
  UseMethod("georeferencedPlot")
}

georeferencedPlot.netTools = function(obj){
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

clockwise <- function(obj){
  UseMethod("clockwise")
}

clockwise.netTools <- function(obj) {
  df <- obj$df
  
  x.coords <- c(df[[1]], df[[1]][1])
  y.coords <- c(df[[2]], df[[2]][1])
  
  double.area <- sum(sapply(2:length(x.coords), function(i) {
    (x.coords[i] - x.coords[i-1])*(y.coords[i] + y.coords[i-1])
  }))
  
  double.area > 0
} 

perpenidularDistance <- function(obj){
  UseMethod("perpenidularDistance")
}

perpenidularDistance.netTools <- function(obj){
  p1 <- obj$p1
  p2 <- obj$p2
  ep <- obj$ep
  
  v1 <- p1 - p2
  v2 <- p2 - ep
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  
  d
}