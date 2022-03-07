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


GeoreferencedPlot <- function(obj, ...){
  UseMethod("GeoreferencedPlot")
}


GeoreferencedGgplot2 <- function(obj, ...){
  UseMethod("GeoreferencedGgplot2")
}


PointToLine <- function(obj){
  UseMethod("PointToLine")
}

PointToSegment_deprecated <- function(obj){
  UseMethod("PointToSegment_deprecated")
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
#' @param obj netTools object -> list(intnet: intensitynet object, vertex_labels: list of labels for the vertices,
#' edge_labels: list of labels for the edges, xy_axes: boolean to show or not the x and y axes, 
#' enable_grid: boolean to draw or not a background grid, show_events: boolean to show or not the events as orange squares)
#' @param ... extra arguments for the plot
#' 
GeoreferencedPlot.netTools <- function(obj, ...){
  g <- obj$intnet$graph
  distances_mtx <- obj$intnet$distances_mtx
  
  arguments <- list(...)
  
  if(!is.null(distances_mtx)){
    node_coords <- matrix(cbind(igraph::vertex_attr(g)$xcoord, igraph::vertex_attr(g)$ycoord), ncol=2)
    
    min_x <- min(node_coords[,1])
    max_x <- max(node_coords[,1])
    min_y <- min(node_coords[,2])
    max_y <- max(node_coords[,2])
    x_range <- c(min_x, max_x)
    y_range <- c(min_y, max_y)
    
    x_dist <- max_x - min_x
    y_dist <- max_y - min_y
    
    n_lines <- 6
    margin <- (x_dist + y_dist)/2 * 0.05
    
    igraph::plot.igraph(g,
                        layout = node_coords,
                        rescale = FALSE,
                        xlim = x_range,
                        ylim = y_range,
                        vertex.color = 'blue',
                        vertex.label = obj$vertex_labels, 
                        vertex.label.cex = if(exists('vertex.label.cex', where=arguments)) arguments[['vertex.label.cex']] else 0.3,  
                        vertex.size = if(exists('vertex.size', where=arguments)) arguments[['vertex.size']] else 0.8 * max(x_range, y_range),
                        edge.label = obj$edge_labels, 
                        edge.label.cex = if(exists('edge.label.cex', where=arguments)) arguments[['edge.label.cex']] else 0.3,
                        edge.arrow.size = if(exists('edge.arrow.size', where=arguments)) arguments[['edge.arrow.size']] else 0.1,
                        ...)
    
    # X and Y coordinates
    if(obj$xy_axes){
      # Square encapsulating the plot
      graphics::rect(min_x - margin, min_y - margin, max_x + margin, max_y + margin)
      
      # X coordinates
      graphics::text(x = sum(x_range) / 2, y = min_y - margin * 3, label = expression(bold("x-coordinate")), adj = 0.5)
      
      # Y coordinates
      graphics::text(x =  min_x - margin * 3, y = sum(y_range) / 2, label = expression(bold("y-coordinate")), srt = 90, adj = 0.5)
      
      
      for(i in 0:(n_lines)){
        # X
        graphics::text(x = min_x + i * x_dist / n_lines, 
                       y =  min_y - margin * 2, 
                       label = floor(min_x + i * x_dist / n_lines))
        # Y
        graphics::text(x = min_x - margin * 2,
                       y =  min_y + i * y_dist / n_lines, 
                       label = floor(min_y + i * y_dist / n_lines),
                       srt = 90)
      }
    }
    
    #grid (if specified)
    if(obj$enable_grid){
      grid_col <- grDevices::rgb(0,0,0,alpha=0.2)
      
      for(i in 0:n_lines){
        # X
        graphics::lines(c(min_x - margin, max_x + margin), 
                        c(min_y + i * y_dist / n_lines, min_y + i * y_dist / n_lines), 
                        col = grid_col)
        # Y
        graphics::lines(c(min_x + i * x_dist / n_lines, min_x + i * x_dist / n_lines), 
                        c(min_y - margin, max_y + margin), 
                        col = grid_col) 
      } 
    }
    
    if(obj$show_events){
      tmp_g <- igraph::make_empty_graph(n = length(obj$intnet$events), directed = FALSE)
      igraph::plot.igraph(tmp_g, 
                          layout = obj$intnet$events,
                          rescale = FALSE,
                          xlim = c(min_x, max_x),
                          ylim = c(min_y, max_y),
                          vertex.color = 'orange',
                          vertex.label = '', 
                          vertex.label.cex = 0.3,  
                          vertex.size = if(exists('vertex.size', where=arguments)) arguments[['vertex.size']] else 0.8 * max(x_range, y_range), 
                          vertex.shape = "square",
                          add = TRUE
      )
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
#' @param obj netTools object -> list(intnet: intensitynet object, 
#' data_df: dataframe(intensity: intensity of the nodes, xcoord: x coordinates of the nodes, 
#' ycoord: y coordinates of the nodes, heattype: data which the heatmap will refer,
#' mode: ('moran', 'getis' or 'intensity'), show_events: boolean to show or not the events as orange squares)
#' @param ... extra arguments for the ggplot
#' 
GeoreferencedGgplot2.netTools <- function(obj, ...){
  arguments <- list(...)
  
  g <- obj$intnet$graph
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
    hplot <- ggplot2::ggplot(data_df, ggplot2::aes_string(x = 'xcoord', y = 'ycoord'), ...) +
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
    hplot <- ggplot2::ggplot(data_df, ggplot2::aes_string('xcoord', 'ycoord'), ...) +
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
    hplot <- ggplot2::ggplot(data_df, ggplot2::aes_string('xcoord', 'ycoord'), ...) +
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
      
      hplot <- ggplot2::ggplot(data_df, ggplot2::aes_string(x = 'xcoord', y = 'ycoord'), ...) +
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
      
      hplot <- ggplot2::ggplot(data_df, ggplot2::aes_string(x = 'xcoord', y = 'ycoord'), ...) +
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
    
    hplot <- ggplot2::ggplot(data_df, ggplot2::aes_string(x = 'xcoord', y = 'ycoord'), ...) +
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
    hplot <- ggplot2::ggplot(data_df, ggplot2::aes_string(x = 'xcoord', y = 'ycoord'), ...) + 
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
  
  if(obj$show_events){
    hplot + ggplot2::geom_point(data = as.data.frame(obj$intnet$events),
                                mapping = ggplot2::aes(x = xcoord, y = ycoord),
                                shape = 22, fill = 'orange', color = 'orange')
  }else{
    hplot
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
#' @name PointToSegment_deprecated.netTools  
#'
#' @param obj netTools object -> list(p1:c(coordx, coordy), p2:c(coordx, coordy), e:c(coordx, coordy))
#' 
#' @return distance to the segment
#' 
PointToSegment_deprecated <- function(obj) {
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
  return(sqrt(dx * dx + dy * dy))
}


#' Return the shortest distance between an event and a set of segments.
#'
#' @name PointToSegment.netTools  
#'
#' @param obj netTools object -> list(p1:matrix(coordx, coordy), p2:matrix(coordx, coordy), e:matrix(coordx, coordy))
#' 
#' @return distance vector to each segment
#' 
PointToSegment <- function(obj) {
  p1 <- obj$p1
  p2 <- obj$p2
  ep <- obj$ep
  
  if(!is.matrix(p1)){
    if(is.data.frame(p1)) p1 <- data.matrix(p1)
    else p1 <- matrix(p1, ncol = 2) 
  } 
  if(!is.matrix(p2)){
    if(is.data.frame(p2)) p2 <- data.matrix(p2)
    else p2 <- matrix(p2, ncol = 2) 
  } 
  if(!is.matrix(ep)){
    if(is.data.frame(ep)) ep <- data.matrix(ep)
    else ep <- matrix(ep, ncol = 2) 
  } 
  
  A <- ep[,1] - p1[,1]
  B <- ep[,2] - p1[,2]
  C <- p2[,1] - p1[,1]
  D <- p2[,2] - p1[,2]
  
  dot <- A * C + B * D
  len_sq <- C * C + D * D
  
  param <- ifelse(len_sq != 0, dot / len_sq, -1)
  
  
  xx <- ifelse(param < 0, p1[,1], ifelse(param > 1, p2[,1], p1[,1] + param * C))
  yy <- ifelse(param < 0, p1[,2], ifelse(param > 1, p2[,2], p1[,2] + param * D))
  
  
  dx <- ep[,1] - xx
  dy <- ep[,2] - yy
  
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