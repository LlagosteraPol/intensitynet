#----------------------------------------UseMethod's----------------------------------------

InitGraph <- function(obj){
  UseMethod("InitGraph")
}


CalculateDistancesMtx <- function(obj){
  UseMethod("CalculateDistancesMtx")
}


SetNetCoords <- function(obj){
  UseMethod("SetNetCoords")
}


SetEdgeIntensity <- function(obj){
  UseMethod("SetEdgeIntensity")
}


SetNodeIntensity <- function(obj){
  UseMethod("SetNodeIntensity")
}


GeoreferencedPlot <- function(obj, vertex_labels='', edge_labels='', xy_axes=TRUE, enable_grid=FALSE, ...){
  UseMethod("GeoreferencedPlot")
}


GeoreferencedGgplot2 <- function(obj, ...){
  UseMethod("GeoreferencedGgplot2")
}


PointToLine <- function(obj){
  UseMethod("PointToLine")
}


PointToSegment <- function(obj){
  UseMethod("PointToSegment")
}


Undirected2RandomDirectedAdjMtx <- function(obj){
  UseMethod("Undirected2RandomDirectedAdjMtx")
}
#-------------------------------------------------------------------------------------------


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
  if(is.null(colnames(weighted_mtx))){
    colnames(weighted_mtx) <- sprintf("V%s", seq(1:nrow(weighted_mtx)))
  } 
  if(graph_type == 'undirected'){
    g <- igraph::graph_from_adjacency_matrix(weighted_mtx, mode = graph_type, weighted = TRUE)
  } else {
    g <- igraph::graph_from_adjacency_matrix(weighted_mtx, mode = 'directed', weighted = TRUE)
  } 
  
  net_coords <- list(graph = g, node_coords = node_coords)
  class(net_coords) <- "netTools"
  g <- SetNetCoords(net_coords)
  
  # Delete isolated vertices
  igraph::delete.vertices(g, igraph::degree(g)==0)
  
  g # return
}


#' Calculates the distances between all pairs of nodes from the given network
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
  
  distances_mtx <- spatstat.geom::pairdist(
    spatstat.geom::ppp(x_coord_node,
                       y_coord_node,
                       xrange = c(min(as.numeric(x_coord_node)), max(as.numeric(x_coord_node))),
                       yrange = c(min(as.numeric(y_coord_node)), max(as.numeric(y_coord_node))))
  )
  rownames(distances_mtx) <- colnames(distances_mtx) <- sprintf("V%s", seq(1:ncol(distances_mtx)))
  distances_mtx
}


#' Set igraph network node coordinates as its attributes
#'
#' @name InitGraph.netTools 
#'
#' @param obj netTools object -> list(graph: igraph, list(): with the node coordinates 'x' and 'y') 
#' 
#' @return igraph network with the given coordinates as the attributes of the nodes
#' 
SetNetCoords.netTools <- function(obj){
  g <- obj$graph
  x_coord_node <- obj$node_coords[, 1]
  y_coord_node <- obj$node_coords[, 2]
  
  # g <- g %>% 
  #      igraph::set_vertex_attr(name = "xcoord", value = x_coord_node) %>% 
  #      igraph::set_vertex_attr(name = "ycoord", value = y_coord_node)
  g <- igraph::set_vertex_attr(g, name = "xcoord", value = x_coord_node)
  g <- igraph::set_vertex_attr(g, name = "ycoord", value = y_coord_node)
  g
}


#' Sets the given intensities as an edge attribute to the given igraph network
#'
#' @name SetEdgeIntensity.netTools 
#'
#' @param obj netTools object -> list(graph: igraph, node_id1: node id, node_id2: node id, intensity: edge intensity)
#' 
#' @return igraph network with the given intensities as attributes of the edges
#' 
SetEdgeIntensity.netTools <- function(obj){
  g <- obj$graph
  node_id1 <- obj$node_id1
  node_id2 <- obj$node_id2
  intensity <- obj$intensity
  
  edge_id <- igraph::get.edge.ids(g, c(node_id1, node_id2))
  #g <- g %>% igraph::set_edge_attr(name = "intensity", index = edge_id, value = intensity)
  g <- igraph::set_edge_attr(g, name = "intensity", index = edge_id, value = intensity)
  g
}


#' Sets the given intensities as a node attribute to the given igraph network
#'
#' @name SetNodeIntensity.netTools 
#'
#' @param obj netTools object -> list(graph: igraph, node_id: node id, intensity: node intensity)
#' 
#' @return igraph network with the given intensities as attributes of the nodes
#' 
SetNodeIntensity.netTools <- function(obj){
  g <- obj$graph
  node_id <- obj$node_id
  intensity <- obj$intensity
  
  #g <- g %>% igraph::set_vertex_attr(name = "intensity", index = node_id, value = intensity)
  g <- igraph::set_vertex_attr(g, name = "intensity", index = node_id, value = intensity)
  g
}


#' Plot the given network using its node coordinates
#'
#' @name GeoreferencedPlot.netTools 
#'
#' @param obj netTools object -> list(graph: igraph, distances_mtx: distances matrix))
#' @param vertex_labels list -> labels for the vertices
#' @param edge_labels list -> labels for the edges
#' @param xy_axes show the x and y axes
#' @param enable_grid draw a background grid
#' @param ... extra arguments for the plot
#' 
GeoreferencedPlot.netTools <- function(obj, vertex_labels='', edge_labels='', xy_axes=TRUE, enable_grid=FALSE, ...){
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  arguments <- list(...)
  
  if(!is.null(distances_mtx)){
    node_coords <- matrix(cbind(igraph::vertex_attr(g)$xcoord, igraph::vertex_attr(g)$ycoord), ncol=2)
    
    min_x <- min(node_coords[,1])
    max_x <- max(node_coords[,1])
    min_y <- min(node_coords[,2])
    max_y <- max(node_coords[,2])
    
    x_dist <- max_x - min_x
    y_dist <- max_y - min_y
    
    igraph::plot.igraph(g,
                        layout = node_coords, 
                        vertex.label = vertex_labels, 
                        vertex.label.cex = if(exists('vertex.label.cex', where=arguments)) arguments[['vertex.label.cex']] else 0.3,  
                        vertex.size = if(exists('vertex.size', where=arguments)) arguments[['vertex.size']] else 2, 
                        edge.label = edge_labels, 
                        edge.label.cex = if(exists('edge.label.cex', where=arguments)) arguments[['edge.label.cex']] else 0.3,
                        edge.arrow.size = if(exists('edge.arrow.size', where=arguments)) arguments[['edge.arrow.size']] else 0.1,
                        ...)
    
    # Square encapsulating the plot
    graphics::rect(-1.05, -1.05, 1.05, 1.05)
    
    # X and Y coordinates
    if(xy_axes){
      graphics::mtext(expression(bold("x-coordinate")), at = 0, line = -28, cex = 0.70)
      graphics::mtext(floor(min_x + (1 * (x_dist/6))), at = -0.67, line = -27, cex = 0.70)
      graphics::mtext(floor(min_x + (2 * (x_dist/6))), at = -0.34, line = -27, cex = 0.70)
      graphics::mtext(floor(min_x + (3 * (x_dist/6))), at = 0, line = -27, cex = 0.70)
      graphics::mtext(floor(min_x + (4 * (x_dist/6))), at = 0.34, line = -27, cex = 0.70)
      graphics::mtext(floor(min_x + (5 * (x_dist/6))), at = 0.67, line = -27, cex = 0.70)
      
      graphics::mtext(expression(bold("y-coordinate")), at = 0, line = 0, cex = 0.70, side = 2)
      graphics::mtext(floor(min_y + (1 * (y_dist/6))), at = -0.67, line = -1, cex = 0.70, side = 2)
      graphics::mtext(floor(min_y + (2 * (y_dist/6))), at = -0.34, line = -1, cex = 0.70, side = 2)
      graphics::mtext(floor(min_y + (3 * (y_dist/6))), at = 0, line = -1, cex = 0.70, side = 2)
      graphics::mtext(floor(min_y + (4 * (y_dist/6))), at = 0.34, line = -1, cex = 0.70, side = 2)
      graphics::mtext(floor(min_y + (5 * (y_dist/6))), at = 0.67, line = -1, cex = 0.70, side = 2)
    }
    
    #grid (if specified)
    if(enable_grid){
      grid_col <- grDevices::rgb(0,0,0,alpha=0.2)
      # X
      graphics::lines(c(-1.05,1.05), c(-1,-1), col = grid_col)
      graphics::lines(c(-1.05,1.05), c(-1 + (1/3), -1 + (1/3)), col = grid_col)
      graphics::lines(c(-1.05,1.05), c(-1 + (2/3), -1 + (2/3)), col = grid_col)
      graphics::lines(c(-1.05,1.05), c(-1 + (3/3), -1 + (3/3)), col = grid_col)
      graphics::lines(c(-1.05,1.05), c(1 - (1/3), 1 - (1/3)), col = grid_col)
      graphics::lines(c(-1.05,1.05), c(1 - (2/3), 1 - (2/3)), col = grid_col)
      graphics::lines(c(-1.05,1.05), c(1 - (3/3), 1 - (3/3)), col = grid_col)
      graphics::lines(c(-1.05,1.05), c(1,1), col = grid_col)
      
      # Y
      graphics::lines(c(-1,-1), c(-1.05,1.05), col = grid_col)
      graphics::lines(c(-1 + (1/3), -1 + (1/3)), c(-1.05,1.05), col = grid_col)
      graphics::lines(c(-1 + (2/3), -1 + (2/3)), c(-1.05,1.05), col = grid_col)
      graphics::lines(c(-1 + (3/3), -1 + (3/3)), c(-1.05,1.05), col = grid_col)
      graphics::lines(c(1 - (1/3), 1 - (1/3)), c(-1.05,1.05), col = grid_col)
      graphics::lines(c(1 - (2/3), 1 - (2/3)), c(-1.05,1.05), col = grid_col)
      graphics::lines(c(1 - (3/3), 1 - (3/3)), c(-1.05,1.05), col = grid_col)
      graphics::lines(c(1,1), c(-1.05,1.05), col = grid_col)
    }
  }
  else{
    igraph::plot.igraph(g, 
                        vertex.label = NA, 
                        vertex.size = 2,
                        vertex.size2 = 2)
  }
}


#' This function uses 'ggplot' to plot heatmaps of a network
#' 
#' @name GeoreferencedGgplot2.netTools
#' 
#' @param obj netTools object -> list(graph: igraph, data_df: dataframe(intensity: intensity of the nodes, 
#' xcoord: x coordinates of the nodes, ycoord: y coordinates of the nodes, heattype: data wich the heatmap will refer), 
#' mode: ('moran', 'getis' or 'intensity'))
#' @param ... extra arguments for the ggplot
#' 
GeoreferencedGgplot2.netTools <- function(obj, ...){
  arguments <- list(...)
  
  g <- obj$graph
  data_df <- obj$data_df
  mode <- obj$mode
  highlighted_df <- data_df[as.numeric(obj$net_vertices),]
  
  node_coords <- data.frame(xcoord = igraph::vertex_attr(g)$xcoord, ycoord = igraph::vertex_attr(g)$ycoord)
  rownames(node_coords) <- igraph::vertex_attr(g)$name
  #get edges, which are pairs of node IDs
  edgelist <- igraph::get.edgelist(g)
  #convert to a four column edge data frame with source and destination coordinates
  edges_df <- data.frame(node_coords[edgelist[,1],], node_coords[edgelist[,2],])
  colnames(edges_df) <- c("xcoord1","ycoord1","xcoord2","ycoord2")
  
  #if(is.null(data_df$intensity) || is.na(data_df$heattype)){
  if(mode == 'moran') {
    # ggplot2::ggplot(data_df, ggplot2::aes(xcoord, ycoord), ...) + 
    #   ggplot2::ggplot2::geom_point(shape = 19,
    #              size = 1.5,
    #              colour="gray") +
    #   ggplot2::ggplot2::geom_point(data = highlighted_df,
    #              shape = 19,
    #              size = 1.5,
    #              aes(xcoord, ycoord, colour = value)) +
    #   viridis::scale_color_viridis() +
    #   ggplot2::geom_tile(ggplot2::aes( fill = as.factor(value) ), 
    #             show.legend = FALSE) + 
    #   ggplot2::labs( title = 'Moran-i Heatmap\n',
    #         color = 'Correlation') +
    #   ggplot2::geom_segment(ggplot2::aes(x = xcoord1, y = ycoord1, xend = xcoord2, yend = ycoord2), 
    #                data = edges, 
    #                size = 0.5, 
    #                colour = "grey") +
    #   ggplot2::scale_y_continuous(name = "y-coordinate") + 
    #   ggplot2::scale_x_continuous(name = "x-coordinate") + 
    #   ggplot2::theme_bw() +
    #   theme( plot.title = ggplot2::element_text(size = 14, 
    #                                    face = "bold", 
    #                                    hjust = 0.5) )
    ggplot2::ggplot(data_df, ggplot2::aes_string(x = 'xcoord', y = 'ycoord'), ...) +
      ggplot2::geom_point(ggplot2::aes_string( colour = 'as.factor(value)' ),
                          shape = 19,
                          size = 1.5) +
      ggplot2::geom_tile(ggplot2::aes_string( fill = 'as.factor(value)' ),
                         show.legend = FALSE) +
      ggplot2::labs( title = 'Moran-i Heatmap\n' ) +
      ggplot2::scale_color_manual(values = c("black", "gray", "skyblue", "yellow", "darkorange", "red4"),
                                  name = "", breaks=c(1,2,3,4,5,6),
                                  labels = c("Not contemplated","insignificant","low-low","low-high","high-low","high-high") ) +
      ggplot2::geom_segment(ggplot2::aes_string(x = 'xcoord1', y = 'ycoord1', 
                                                xend = 'xcoord2', yend = 'ycoord2'),
                            data = edges_df,
                            size = 0.5,
                            colour = "grey") +
      ggplot2::scale_y_continuous(name = "y-coordinate") +
      ggplot2::scale_x_continuous(name = "x-coordinate") +
      ggplot2::theme_bw() +
      ggplot2::theme( plot.title = ggplot2::element_text(size = 14,
                                                         face = "bold",
                                                         hjust = 0.5) )
  }else if(mode == 'geary'){
    ggplot2::ggplot(data_df, ggplot2::aes_string('xcoord', 'ycoord'), ...) +
      ggplot2::geom_point( ggplot2::aes_string( colour = 'as.factor(value)' ),
                           shape = 19,
                           size = 1.5 ) +
      ggplot2::geom_tile( ggplot2::aes_string( fill = 'as.factor(value)' ),
                          show.legend = FALSE ) +
      ggplot2::labs(title = 'Geary-c Heatmap\n') +
      ggplot2::scale_color_manual(values = c("black", "green", "gray", "red"),
                                  name = "",
                                  breaks = c(1,2,3, 4),
                                  labels = c("Not contemplated", "positive auto.","no auto.","negative auto.") ) +
      ggplot2::geom_segment(ggplot2::aes_string(x = 'xcoord1', y = 'ycoord1', 
                                                xend = 'xcoord2', yend = 'ycoord2'),
                            data = edges_df,
                            size = 0.5,
                            colour = "grey") +
      ggplot2::scale_y_continuous( name = "y-coordinate" ) +
      ggplot2::scale_x_continuous( name = "x-coordinate" ) +
      ggplot2::theme_bw() +
      ggplot2::theme( plot.title = ggplot2::element_text( size = 14,
                                                          face = "bold",
                                                          hjust = 0.5 ) )
    # ggplot2::ggplot(data_df, ggplot2::aes(xcoord, ycoord), ...) +
    #   ggplot2::geom_point(shape = 19,
    #                       size = 1.5,
    #                       colour="gray") +
    #   ggplot2::geom_point(data = highlighted_df,
    #                       shape = 19,
    #                       size = 1.5,
    #                       ggplot2::aes(xcoord, ycoord, colour = value)) +
    #   viridis::scale_color_viridis() +
    #   ggplot2::labs(title = 'Geary-c Heatmap\n',
    #                 color = 'Correlation') +
    #   ggplot2::geom_segment(ggplot2::aes(x = xcoord1, y = ycoord1, xend = xcoord2, yend = ycoord2),
    #                         data = edges_df,
    #                         size = 0.5,
    #                         colour="grey") +
    #   ggplot2::scale_y_continuous(name = "y-coordinate") +
    #   ggplot2::scale_x_continuous(name = "x-coordinate") +
    #   ggplot2::theme_bw() +
    #   ggplot2::theme(legend.title = ggplot2::element_text(face = "bold"),
    #                  plot.title = ggplot2::element_text( size = 14,
    #                                                      face = "bold",
    #                                                      hjust = 0.5) )
  }else if(mode == 'getis'){
    #TODO: implement
    
  }else if( mode == 'v_intensity' ){
    ggplot2::ggplot(data_df, ggplot2::aes_string('xcoord', 'ycoord'), ...) +
      ggplot2::geom_point(shape = 19,
                          size = 1.5,
                          colour="gray") +
      ggplot2::geom_point(data = highlighted_df,
                          shape = 19,
                          size = 1.5,
                          ggplot2::aes_string(x = 'xcoord', y = 'ycoord', colour = 'value')) +
      viridis::scale_color_viridis() +
      ggplot2::labs(title = 'Vertex Intensity Heatmap\n',
                    color = 'Norm. intensity') +
      ggplot2::geom_segment(ggplot2::aes_string(x = 'xcoord1', y = 'ycoord1', 
                                                xend = 'xcoord2', yend = 'ycoord2'),
                            data = edges_df,
                            size = 0.5,
                            colour="grey") +
      ggplot2::scale_y_continuous(name = "y-coordinate") +
      ggplot2::scale_x_continuous(name = "x-coordinate") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title = ggplot2::element_text(face = "bold"),
                     plot.title = ggplot2::element_text( size = 14,
                                                         face = "bold",
                                                         hjust = 0.5) )
  }else if (mode == 'e_intensity'){
    if(length(obj$net_vertices) == length(igraph::V(g))){
      edge_int <- igraph::edge_attr(g, 'intensity')
      norm_int <- (edge_int - min(edge_int)) / (max(edge_int) - min(edge_int))
      
      ggplot2::ggplot(data_df, ggplot2::aes_string(x = 'xcoord', y = 'ycoord'), ...) +
        ggplot2::geom_point(shape = 19,
                            size = 1.5,
                            colour="gray") +
        viridis::scale_color_viridis() +
        ggplot2::labs(title = 'Edge Intensity Heatmap\n',
                      color = 'Norm. intensity') +
        ggplot2::geom_segment(ggplot2::aes_string(x = 'xcoord1', y = 'ycoord1', 
                                                  xend = 'xcoord2', yend = 'ycoord2',
                                                  colour = 'norm_int'),
                              data = edges_df,
                              size = 0.5) +
        ggplot2::scale_y_continuous(name = "y-coordinate") +
        ggplot2::scale_x_continuous(name = "x-coordinate") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.title = ggplot2::element_text(face = "bold"),
                       plot.title = ggplot2::element_text( size = 14,
                                                           face = "bold",
                                                           hjust = 0.5) )
    }else{
      sub_g <- igraph::induced_subgraph(graph = g, vids = obj$net_vertices)
      sub_edges <- igraph::ends(sub_g, es = igraph::E(sub_g))
      
      highlighted_edges <- igraph::get.edge.ids(g, c(t(sub_edges)))
      edge_int <- igraph::edge_attr(g, 'intensity', highlighted_edges)
      norm_int <- (edge_int - min(edge_int)) / (max(edge_int) - min(edge_int))
      
      #convert to a four column edge data frame with source and destination coordinates
      sub_edges_df <- data.frame(node_coords[sub_edges[,1],], node_coords[sub_edges[,2],])
      colnames(sub_edges_df) <- c("xcoord1","ycoord1","xcoord2","ycoord2")
      
      ggplot2::ggplot(data_df, ggplot2::aes_string(x = 'xcoord', y = 'ycoord'), ...) +
        ggplot2::geom_point(shape = 19,
                            size = 1.5,
                            colour="gray") +
        viridis::scale_color_viridis() +
        ggplot2::labs(title = 'Edge Intensity Heatmap\n',
                      color = 'Norm. intensity') +
        ggplot2::geom_segment(ggplot2::aes_string(x = 'xcoord1', y = 'ycoord1', 
                                                  xend = 'xcoord2', yend = 'ycoord2'),
                              data = edges_df,
                              size = 0.5,
                              colour = 'grey') +
        ggplot2::geom_segment(ggplot2::aes_string(x = 'xcoord1', y = 'ycoord1', 
                                                  xend = 'xcoord2', yend = 'ycoord2', 
                                                  colour = 'norm_int'),
                              data = sub_edges_df,
                              size = 0.5) +
        ggplot2::scale_y_continuous(name = "y-coordinate") +
        ggplot2::scale_x_continuous(name = "x-coordinate") +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.title = ggplot2::element_text(face = "bold"),
                       plot.title = ggplot2::element_text( size = 14,
                                                           face = "bold",
                                                           hjust = 0.5) )
    }
  }else if(mode == 'intensity'){
    sub_g <- igraph::induced_subgraph(graph = g, vids = obj$net_vertices)
    sub_edges <- igraph::ends(sub_g, es = igraph::E(sub_g))
    
    highlighted_edges <- igraph::get.edge.ids(g, c(t(sub_edges)))
    edge_int <- igraph::edge_attr(g, 'intensity', highlighted_edges)
    norm_e_int <- (edge_int - min(edge_int)) / (max(edge_int) - min(edge_int))
    
    #convert to a four column edge data frame with source and destination coordinates
    sub_edges_df <- data.frame(node_coords[sub_edges[,1],], node_coords[sub_edges[,2],])
    colnames(sub_edges_df) <- c("xcoord1","ycoord1","xcoord2","ycoord2")
    
    ggplot2::ggplot(data_df, ggplot2::aes_string(x = 'xcoord', y = 'ycoord'), ...) +
      ggplot2::geom_point(shape = 19,
                          size = 1.5,
                          colour="gray") +
      ggplot2::geom_point(data = highlighted_df,
                          shape = 19,
                          size = 1.5,
                          ggplot2::aes_string(x = 'xcoord', y = 'ycoord', colour = 'value')) +
      viridis::scale_color_viridis() +
      ggplot2::labs(title = 'Intensity Heatmap\n',
                    color = 'Norm. intensity') +
      ggplot2::geom_segment(ggplot2::aes_string(x = 'xcoord1', y = 'ycoord1', 
                                                xend = 'xcoord2', yend = 'ycoord2'),
                            data = edges_df,
                            size = 0.5,
                            colour = 'grey') +
      ggplot2::geom_segment(ggplot2::aes_string(x = 'xcoord1', y = 'ycoord1', 
                                                xend = 'xcoord2', yend = 'ycoord2', 
                                                colour = 'norm_e_int'),
                            data = sub_edges_df,
                            size = 0.5) +
      ggplot2::scale_y_continuous(name = "y-coordinate") +
      ggplot2::scale_x_continuous(name = "x-coordinate") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.title = ggplot2::element_text(face = "bold"),
                     plot.title = ggplot2::element_text( size = 14,
                                                         face = "bold",
                                                         hjust = 0.5) )
  }else{
    ggplot2::ggplot(data_df, ggplot2::aes_string(x = 'xcoord', y = 'ycoord'), ...) + 
      ggplot2::geom_point(shape = 19, 
                          size = 1.5) +
      ggplot2::geom_segment(ggplot2::aes_string(x = 'xcoord1', y = 'ycoord1', 
                                                xend = 'xcoord2', yend = 'ycoord2'), 
                            data = edges_df, 
                            size = 0.5, 
                            colour = "grey") +
      ggplot2::scale_y_continuous(name = "y-coordinate") + 
      ggplot2::scale_x_continuous(name = "x-coordinate") + 
      ggplot2::theme_bw()
  }
}


#' Return the distance between an event and the line (not segment) formed by two nodes.
#'
#' @name PointToLine.netTools  
#'
#' @param obj netTools object -> list(p1:c(coordx, coordy), p2:c(coordx, coordy), e:c(coordx, coordy))
#' 
#' @return the distance to the line
#' 
PointToLine.netTools <- function(obj){
  p1 <- obj$p1
  p2 <- obj$p2
  ep <- obj$ep
  
  v1 <- p1 - p2
  v2 <- p2 - ep
  m <- cbind(v1,v2)
  d <- abs(det(m)) / sqrt(sum(v1 * v1))
  
  d
}


#' Return the shortest distance between an event and the segment formed by two nodes.
#'
#' @name PointToSegment.netTools  
#'
#' @param obj netTools object -> list(p1:c(coordx, coordy), p2:c(coordx, coordy), e:c(coordx, coordy))
#' 
#' @return distance to the segment
#' 
PointToSegment <- function(obj) {
  #start_time <- Sys.time() # debug only
  p1 <- obj$p1
  p2 <- obj$p2
  ep <- obj$ep
  A <- ep[1] - p1[1]
  B <- ep[2] - p1[2]
  C <- p2[1] - p1[1]
  D <- p2[2] - p1[2]
  
  dot <- A * C + B * D
  len_sq <- C * C + D * D
  param <- -1
  if (len_sq != 0){
    param <- dot / len_sq # in case of 0 length line
  } 
  
  if (param < 0) {
    xx <- p1[1]
    yy <- p1[2]
  }
  else if (param > 1) {
    xx <- p2[1]
    yy <- p2[2]
  }
  else {
    xx <- p1[1] + param * C
    yy <- p1[2] + param * D
  }
  
  dx <- ep[1] - xx
  dy <- ep[2] - yy
  #cat(paste0("PointToSegment time: ", Sys.time() - start_time, "\n")) # debug only
  return(sqrt(dx * dx + dy * dy))
}


#' Creates a directed adjacency matrix from an Undirected one with random directions (in-out edges) 
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
        if(stats::runif(1) <= prob) mtx[row, col] <- 0
      }
    }
  }
  mtx
}