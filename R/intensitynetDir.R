#' Given a node, calculates its mean intensities regarding in and out edges associated with the node.
#' 
#' @name nodeIntensity.intensitynetDir
#' 
#' @param obj intensitynetDir object
#' @param node_id ID of the node
#' 
#' @return mean intensities of the given node for in and out edges
#' 
MeanNodeIntensity.intensitynetDir= function(obj, node_id){
  g <- obj$graph
  
  # If the intensities are already calculated, return them
  if(!is.null(igraph::vertex_attr(g, 'intensity_in', index = node_id)) &
     !is.null(igraph::vertex_attr(g, 'intensity_out', index = node_id))) {
    
    if(!is.na(igraph::vertex_attr(g, "intensity_in", index = node_id))[1] &
       !is.na(igraph::vertex_attr(g, "intensity_out", index = node_id))[1]) {
      
      return( list(in_int  = igraph::vertex_attr(g, 'intensity_in', index=node_id),
                   out_int = igraph::vertex_attr(g, 'intensity_out', index=node_id)))
    }
  }
  
  if(igraph::degree(g, node_id) > 0){
    in_neighbors  <- igraph::neighbors(g, node_id, mode = 'in')
    out_neighbors <- igraph::neighbors(g, node_id, mode = 'out')
    
    if(length(in_neighbors) > 0){
      in_mat <- matrix(0, ncol = length(in_neighbors)) 
      colnames(in_mat) <- as.vector(in_neighbors) 
      rownames(in_mat) <- node_id
      
      for (neighbor_id in in_neighbors){
        in_mat[as.character(node_id), 
               as.character(neighbor_id)] <- EdgeIntensity(obj = obj, 
                                                           node_id1 = igraph::V(g)[node_id]$name, 
                                                           node_id2 = igraph::V(g)[neighbor_id]$name,
                                                           z = obj$event_correction)
      }
      in_intensity <- Reduce('+', in_mat) / length(in_neighbors)
    }else{
      in_intensity <- 0
    }
    
    if(length(out_neighbors) > 0){
      out_mat <- matrix(0, ncol = length(out_neighbors)) 
      colnames(out_mat) <- as.vector(out_neighbors) 
      rownames(out_mat) <- node_id
      
      for (neighbor_id in out_neighbors){
        out_mat[as.character(node_id), 
                as.character(neighbor_id)] <- EdgeIntensity(obj = obj, 
                                                            node_id1 = igraph::V(g)[node_id]$name, 
                                                            node_id2 = igraph::V(g)[neighbor_id]$name,
                                                            z = obj$event_correction)
      }
      
      out_intensity <- Reduce('+', out_mat) / length(out_neighbors)
    }else{
      out_intensity <- 0
    }
    
    list(in_int = in_intensity, out_int = out_intensity)
  }
}


#' Calculates edgewise and mean nodewise intensities for Directed networks
#' 
#' @name CalculateEventIntensities.intensitynetDir
#' 
#' @param obj intensitynetDir object
#' 
#' @return intensitynetDir object with a graph containing all the intensities as attributes of its nodes and edges
#' 
#' @export
CalculateEventIntensities.intensitynetDir = function(obj){
  g <- obj$graph
  intensities <- obj$intensities
  in_counts <- c()
  out_counts <- c()
  
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
    
    if(is.null(igraph::vertex_attr(g, 'intensity_in', node_id)) || is.null(igraph::vertex_attr(g, 'intensity_out', node_id))){
      if(igraph::degree(g, node_id) > 0){
        #Adds result of Nodewise mean intensity function to 'counts'
        intensities <- MeanNodeIntensity(tmp_obj, node_id)
        in_counts[[node_id]]  <- intensities$in_int
        out_counts[[node_id]] <- intensities$out_int
      }else{
        # Counts for isolated nodes or NA values
        in_counts[[node_id]]  <- 0
        out_counts[[node_id]] <- 0
      }
    }else if(is.na(igraph::vertex_attr(g, 'intensity_in', node_id))[1] ||
             is.na(igraph::vertex_attr(g, 'intensity_out', node_id))[1]){
      
      if(is.na(igraph::vertex_attr(g, 'intensity_in', node_id))[1]) in_counts[[node_id]]   <- 0
      if(is.na(igraph::vertex_attr(g, 'intensity_out', node_id))[1]) out_counts[[node_id]] <- 0
      
    }else{
      in_counts[[node_id]]  <- igraph::vertex_attr(g, 'intensity_in', node_id)
      out_counts[[node_id]] <- igraph::vertex_attr(g, 'intensity_out', node_id)
    }
  }
  close(pb)
  
  # g <- g %>% igraph::set_vertex_attr(name = "intensity_in", value = as.matrix(in_counts)) %>% 
  #            igraph::set_vertex_attr(name = "intensity_out", value = as.matrix(out_counts))
  g <- igraph::set_vertex_attr(g, name = "intensity_all", value = as.matrix(in_counts))
  g <- igraph::set_vertex_attr(g, name = "intensity_out", value = as.matrix(out_counts))
  
  intnet <- list(graph = g, events = obj$events, graph_type = obj$graph_type, distances_mtx = obj$distances_mtx)
  attr(intnet, 'class') <- c("intensitynet", "intensitynetDir")
  return(intnet)
}


#' Plot intensitynet object
#'
#' @name plot.intensitynetDir
#'
#' @param x intensitynet object
#' @param vertex_labels list -> labels for the vertices
#' @param edge_labels list -> labels for the edges
#' @param xy_axes show the x and y axes
#' @param enable_grid draw a background grid
#' @param show_events option to show the events as orange squares, FALSE by default
#' @param ... extra arguments for the plot
#' 
#' @return No return value, same as graphics::plot.
#' 
#' @examples
#' 
#' data("dir_intnet_chicago")
#' plot(dir_intnet_chicago) # basic plot
#' plot(dir_intnet_chicago, enable_grid = TRUE) # with grid
#' plot(dir_intnet_chicago, xy_axes = FALSE) # without axes
#' 
#' @export
plot.intensitynetDir <- function(x, vertex_labels='none', edge_labels='none', 
                                 xy_axes=TRUE, enable_grid=FALSE, show_events = FALSE, ...){
  g <- x$graph
  
  v_label <- switch(vertex_labels, 
                    none = {''}, 
                    intensity_in = {round(igraph::vertex_attr(g)$intensity_in, 4)},
                    intensity_out = {round(igraph::vertex_attr(g)$intensity_out, 4)},
                    '')
  
  e_label <- switch(edge_labels, 
                    none = {''}, 
                    intensity = {round(igraph::edge_attr(g)$intensity, 4)},
                    '')
  
  geoplot_obj <- list(intnet = x, 
                      vertex_labels = v_label, 
                      edge_labels = e_label, 
                      xy_axes = xy_axes, 
                      enable_grid = enable_grid, 
                      show_events = show_events)
  class(geoplot_obj) <- "netTools"
  
  GeoreferencedPlot(geoplot_obj, ...)
}