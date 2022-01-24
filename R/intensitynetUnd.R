#' Calculates the mean intensity of the given node (intensity of all the edges of the node/number of edges of the node)
#' 
#' @name nodeIntensity.intensitynetUnd
#' 
#' @param obj intensitynetUnd object
#' @param node_id ID of the node
#' 
#' @return mean intensity of the given node
#' 
MeanNodeIntensity.intensitynetUnd = function(obj, node_id){
  g <- obj$graph
  
  # If the intensity is already calculated, return it
  if(!is.null(igraph::vertex_attr(g, 'intensity', index=node_id))){
    if(!is.na(igraph::vertex_attr(g, "intensity", index=node_id))[1]){
      return(igraph::vertex_attr(g, 'intensity', index=node_id))
    } 
  }
  
  if(igraph::degree(g, node_id) > 0){
    neighbors_list <- igraph::neighbors(g, node_id)
    
    ev_mat <- matrix(0, ncol = length(neighbors_list)) 
    colnames(ev_mat) <- as.vector(neighbors_list) 
    rownames(ev_mat) <- node_id
    
    for (neighbor_id in neighbors_list){
      ev_mat[as.character(node_id), as.character(neighbor_id)] <- EdgeIntensity(obj, 
                                                                                igraph::V(g)[node_id]$name, 
                                                                                igraph::V(g)[neighbor_id]$name)
    }
    
    mean_intensity <- Reduce('+', ev_mat) / igraph::degree(g, node_id)
    
    mean_intensity
  }
}


#' Calculates edgewise and mean nodewise intensities for Undirected networks
#' 
#' @name CalculateEventIntensities.intensitynetUnd
#' 
#' @param obj intensitynetUnd object
#' 
#' @return intensitynetUnd object with a graph containing all the intensities as attributes of its nodes and edges
#' 
#' @export
CalculateEventIntensities.intensitynetUnd = function(obj){
  g <- obj$graph
  counts <- c()
  
  if(length(obj$events) == 0){
    warning("No events, cannot calculate any intensity.")
    return(obj)
  }
  
  tmp_obj <- AllEdgeIntensities.intensitynet(obj)
  g <- tmp_obj$graph
  
  pb = utils::txtProgressBar(min = 0, max = igraph::gorder(g), initial = 0) 
  message("Calculating node intensities...")
  
  # check if the intensities was previously calculated, if not, calculate them
  for(node_id in igraph::V(g)){
    
    utils::setTxtProgressBar(pb,node_id)
    
    if(is.null(igraph::vertex_attr(g, 'intensity', node_id))){
      if(igraph::degree(g, node_id) > 0){
        #Adds result of Nodewise mean intensity function to 'counts'
        counts[[node_id]]  <- MeanNodeIntensity(tmp_obj, node_id)
      }else{
        # Counts for isolated nodes or NA values
        counts[[node_id]] <- 0
      }
    }else if(is.na(igraph::vertex_attr(g, 'intensity', node_id))[1]){
      counts[[node_id]] <- 0
    }else{
      counts[[node_id]] <- igraph::vertex_attr(g, 'intensity', node_id)
    }
  }
  close(pb)
  
  #g <- g %>% igraph::set_vertex_attr(name = "intensity", value = as.matrix(counts))
  g <- igraph::set_vertex_attr(g, name = "intensity", value = as.matrix(counts))
  
  intnet <- list(graph = g, events = obj$events, graph_type = obj$graph_type, distances_mtx = obj$distances_mtx)
  attr(intnet, 'class') <- c("intensitynet", "intensitynetUnd")
  return(intnet)
}


#' Plot intensitynet object
#'
#' @name plot.intensitynetUnd
#'
#' @param x intensitynet object
#' @param vertex_labels list -> labels for the vertices
#' @param edge_labels list -> labels for the edges
#' @param xy_axes show the x and y axes
#' @param enable_grid draw a background grid
#' @param ... extra arguments for the plot
#' 
#' @return No return value, same as graphics::plot.
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' plot(und_intnet_chicago) # basic plot
#' plot(und_intnet_chicago, enable_grid = TRUE) # with grid
#' plot(und_intnet_chicago, xy_axes = FALSE) # without axes
#' 
#' @export
plot.intensitynetUnd <- function(x, vertex_labels = 'none', edge_labels = 'none', 
                                 xy_axes = TRUE, enable_grid = FALSE, ...){
  g <- x$graph
  
  v_label <- switch(vertex_labels, 
                    none = {''}, 
                    intensity = {round(igraph::vertex_attr(g)$intensity, 4)},
                    '')
  
  e_label <- switch(edge_labels, 
                    none = {''}, 
                    intensity = {round(igraph::edge_attr(g)$intensity, 4)},
                    '')
  
  geoplot_obj <- list(graph = g, distances_mtx = x$distances_mtx)
  class(geoplot_obj) <- "netTools"
  
  GeoreferencedPlot(geoplot_obj, 
                    vertex_labels = v_label, 
                    edge_labels = e_label, 
                    xy_axes = xy_axes, 
                    enable_grid = enable_grid, 
                    ...)
}