#TODO: Declare all these functions as non-visible at namespace


InitGraph <- function(obj){
  UseMethod("InitGraph")
}

#' Creates an igraph network with the given data
#'
#' @name InitGraph.netTools 
#'
#' @param obj netTools object -> list(adjacency_mtx: graph adjacency matrix, distances: distances between every pair of nodes,
#' graph_type: directed or undirected, node_coords: node coordinates matrix)
#' 
#' @return igraph network
#' 
InitGraph.netTools <- function(obj){
  adjacency_mtx <- obj$adjacency_mtx
  distances_mtx <- obj$distances_mtx
  graph_type <- obj$graph_type
  node_coords <- obj$node_coords
    
  weighted_mtx = adjacency_mtx * distances_mtx
  if(graph_type == 'undirected') g <- graph_from_adjacency_matrix(weighted_mtx, mode = graph_type, weighted=TRUE)
  else  g <- graph_from_adjacency_matrix(weighted_mtx, mode = 'directed', weighted=TRUE)
  
  net_coords <- list(graph = g, node_coords = node_coords)
  class(net_coords) <- "netTools"
  g <- SetNetCoords(net_coords)
  
  # Delete isolated vertices
  igraph::delete.vertices(g, igraph::degree(g)==0)
  
  g # return
}

SetNetCoords <- function(obj){
  UseMethod("SetNetCoords")
}


#' Set igraph network node coordinates as its attributes
#'
#' @name InitGraph.netTools 
#'
#' @param obj netTools object -> list(graph: igraph, list(): with the node coordinates 'x' and 'y') 
#' 
#' @return igraph network with the given coordinates as the attributes of the nodes
#' 
SetNetCoords.netTools = function(obj){
  g <- obj$graph
  x_coord_node <- obj$node_coords[, 1]
  y_coord_node <- obj$node_coords[, 2]
  
  # TODO: change x and y to long (longitude) and lat (latitude), beware of the coordinate system type (e.g WGS84)
  g <- g %>% set_vertex_attr(name = "xcoord", value = x_coord_node) %>% 
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
#' @param obj netTools object -> list(): with the node coordinates 'x' and 'y'
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
#' @param obj netTools object -> list(graph: igraph, node_id1: node id, node_id2: node id, intensity: edge intensity)
#' 
#' @return igraph network with the given intensities as an attributes of the edges
#' 
#TODO: Declare non-visible at namespace
SetEdgeIntensity.netTools <- function(obj){
  g <- obj$graph
  node_id1 <- obj$node_id1
  node_id2 <- obj$node_id2
  intensity <- obj$intensity
  
  edge_id <- get.edge.ids(g, c(node_id1, node_id2))
  g <- g %>% set_edge_attr(name = "intensity", index = edge_id, value = intensity)
  g
}

SetNodeIntensity <- function(obj){
  UseMethod("SetNodeIntensity")
}


#' Sets the given intensites as a node attribute to the given igraph network
#'
#' @name SetNodeIntensity.netTools 
#'
#' @param obj netTools object -> list(graph: igraph, node_id: node id, intensity: node intensity)
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
#' @param obj netTools object -> list(graph: igraph, node_id1: node id, node_id2: node id, distances_mtx: distances matrix))
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

GeoreferencedPlot <- function(obj, node_label='none', edge_label='none', ...){
  UseMethod("GeoreferencedPlot")
}


#' Plot the given network using its node coordinates
#'
#' @name GeoreferencedPlot.netTools 
#'
#' @param obj netTools object -> list(graph: igraph, distances_mtx: distances matrix))
#' 
GeoreferencedPlot.netTools = function(obj, vertex_intensity='', edge_intensity='', xy_axes=TRUE, enable_grid=FALSE, ...){
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  arguments <- list(...)
  
  if(!is.null(distances_mtx)){
    node_coords <- matrix(cbind(vertex_attr(g)$xcoord, vertex_attr(g)$ycoord), ncol=2)
    
    min_x <- min(node_coords[,1])
    max_x <- max(node_coords[,1])
    min_y <- min(node_coords[,2])
    max_y <- max(node_coords[,2])
    
    x_dist <- max_x - min_x
    y_dist <- max_y - min_y

    plot(g, 
           layout=node_coords, 
           vertex.label = vertex_intensity, 
           vertex.label.cex = if(exists('vertex.label.cex', where=arguments)) arguments[['vertex.label.cex']] else 0.3,  
           vertex.size= if(exists('vertex.size', where=arguments)) arguments[['vertex.size']] else 2, 
           edge.label = edge_intensity, 
           edge.label.cex = if(exists('edge.label.cex', where=arguments)) arguments[['edge.label.cex']] else 0.3,
           edge.arrow.size = if(exists('edge.arrow.size', where=arguments)) arguments[['edge.arrow.size']] else 0.1,
           ...)
    
    # Square encapsulating the plot
    rect(-1.05,-1.05,1.05,1.05)
    
    # X and Y coordinates
    if(xy_axes){
      mtext(expression(bold("x-coordinate")), at=0, line = -28, cex=0.70)
      mtext(floor(min_x + (1 * (x_dist/6))), at=-0.67, line = -27, cex=0.70)
      mtext(floor(min_x + (2 * (x_dist/6))), at=-0.34, line = -27, cex=0.70)
      mtext(floor(min_x + (3 * (x_dist/6))), at=0, line = -27, cex=0.70)
      mtext(floor(min_x + (4 * (x_dist/6))), at=0.34, line = -27, cex=0.70)
      mtext(floor(min_x + (5 * (x_dist/6))), at=0.67, line = -27, cex=0.70)
      
      mtext(expression(bold("y-coordinate")), at=0, line=0, cex=0.70, side = 2)
      mtext(floor(min_y + (1 * (y_dist/6))), at=-0.67, line = -1, cex=0.70, side = 2)
      mtext(floor(min_y + (2 * (y_dist/6))), at=-0.34, line = -1, cex=0.70, side = 2)
      mtext(floor(min_y + (3 * (y_dist/6))), at=0, line = -1, cex=0.70, side = 2)
      mtext(floor(min_y + (4 * (y_dist/6))), at=0.34, line = -1, cex=0.70, side = 2)
      mtext(floor(min_y + (5 * (y_dist/6))), at=0.67, line = -1, cex=0.70, side = 2)
    }
    
    #grid (if specified)
    if(enable_grid){
      grid_col <- rgb(0,0,0,alpha=0.2)
      # X
      lines(c(-1.05,1.05), c(-1,-1), col = grid_col)
      lines(c(-1.05,1.05), c(-1 + (1/3), -1 + (1/3)), col = grid_col)
      lines(c(-1.05,1.05), c(-1 + (2/3), -1 + (2/3)), col = grid_col)
      lines(c(-1.05,1.05), c(-1 + (3/3), -1 + (3/3)), col = grid_col)
      lines(c(-1.05,1.05), c(1 - (1/3), 1 - (1/3)), col = grid_col)
      lines(c(-1.05,1.05), c(1 - (2/3), 1 - (2/3)), col = grid_col)
      lines(c(-1.05,1.05), c(1 - (3/3), 1 - (3/3)), col = grid_col)
      lines(c(-1.05,1.05), c(1,1), col = grid_col)
      
      # Y
      lines(c(-1,-1), c(-1.05,1.05), col = grid_col)
      lines(c(-1 + (1/3), -1 + (1/3)), c(-1.05,1.05), col = grid_col)
      lines(c(-1 + (2/3), -1 + (2/3)), c(-1.05,1.05), col = grid_col)
      lines(c(-1 + (3/3), -1 + (3/3)), c(-1.05,1.05), col = grid_col)
      lines(c(1 - (1/3), 1 - (1/3)), c(-1.05,1.05), col = grid_col)
      lines(c(1 - (2/3), 1 - (2/3)), c(-1.05,1.05), col = grid_col)
      lines(c(1 - (3/3), 1 - (3/3)), c(-1.05,1.05), col = grid_col)
      lines(c(1,1), c(-1.05,1.05), col = grid_col)
    }
  }
  else{
    plot(g, 
         vertex.label=NA, 
         vertex.size=2,
         vertex.size2=2)
  }
}


GeoreferencedGgplot2 <- function(obj, ...){
  UseMethod("GeoreferencedGgplot2")
}

GeoreferencedGgplot2.netTools = function(obj, ...){
  arguments <- list(...)
  
  g <- obj$graph
  data_df <- obj$data_df
  mode <- obj$mode
 
  node_coords <- data.frame(xcoord = vertex_attr(g)$xcoord, ycoord = vertex_attr(g)$ycoord)
  rownames(node_coords) <- sprintf("V%s",seq(1:nrow(node_coords)))
  #get edges, which are pairs of node IDs
  edgelist <- get.edgelist(g)
  #convert to a four column edge data frame with source and destination coordinates
  edges <- data.frame(node_coords[edgelist[,1],], node_coords[edgelist[,2],])
  colnames(edges) <- c("xcoord1","ycoord1","xcoord2","ycoord2")
  
  if(is.null(data_df$intensity) || is.na(data_df$heatmap)){
    ggplot(data_df, aes(xcoord,ycoord), ...) + 
      geom_point(shape=19, size=1.5) +
      geom_segment(aes(x=xcoord1, y=ycoord1, xend = xcoord2, yend = ycoord2), 
                   data=edges, 
                   size = 0.5, 
                   colour="grey") +
      scale_y_continuous(name="y-coordinate") + 
      scale_x_continuous(name="x-coordinate") + theme_bw()

  }else{
    if(mode=='moran_i') {
      ggplot(data_df, aes(xcoord,ycoord), ...) + 
        geom_point(aes(colour=as.factor(heatmap)), shape=19, size=1.5) +
        geom_tile(aes(fill=as.factor(heatmap)), show.legend = FALSE) + 
        scale_color_manual(values=c("gray","skyblue", "yellow", "darkorange", "red4"), 
                           name="", breaks=c(0,1,2,3,4), labels=c("insignificant","low-low","low-high","high-low","high-high")) +
        geom_segment(aes(x=xcoord1, y=ycoord1, xend = xcoord2, yend = ycoord2), 
                     data=edges, 
                     size = 0.5, 
                     colour="grey") +
        scale_y_continuous(name="y-coordinate") + 
        scale_x_continuous(name="x-coordinate") + theme_bw()
    }else if(mode=='geary_g'){
      ggplot(data_df, aes(xcoord,ycoord), ...) +  
        geom_point(alpha = 0) + 
        geom_tile()+ 
        geom_text(aes(label=heatmap),hjust=0, vjust=0, size=3, check_overlap = T) +
        scale_colour_grey(guide='none') + 
        geom_segment(aes(x = xcoord1, y = ycoord1, xend = xcoord2, yend = ycoord2), 
                     data = edges, 
                     size = 0.5, 
                     colour = "grey") +
        scale_y_continuous(name="y-coordinate") + 
        scale_x_continuous(name="x-coordinate") + theme_bw() 
    }
  }
}


PointToLine <- function(obj){
  UseMethod("PointToLine")
}


#' Return the perpendicular distance between an event and the line formed by two nodes.
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


Undirected2RandomDirectedAdjMtx <- function(obj){
  UseMethod("Undirected2RandomDirectedAdjMtx")
}

#' Creates a directed adjacency matrix from an Undirected one with random directions (in out edges) 
#' but with the same connections between nodes.
#'
#' @param obj netTools object -> list(mtx: matrix)
#' 
#' @return directed adjacency matrix with random directions
#' 
Undirected2RandomDirectedAdjMtx.netTools  <- function(obj){
  mtx <- obj$mtx
  prob <- 0.25
  
  for(row in 1:nrow(mtx)) {
    for(col in row:ncol(mtx)) {
      if(mtx[row, col] != 0){
        if(runif(1) <= prob) mtx[row, col] <- 0
      }
    }
  }
  mtx
}