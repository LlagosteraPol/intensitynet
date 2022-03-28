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
      
      message("Warning: Node intensities were already calculated in a previous instance, returning the same intensity.")
      return(igraph::vertex_attr(g, 'intensity', index=node_id))
    } 
  }
  
  if(igraph::degree(g, node_id) > 0){
    neighbors_list <- igraph::neighbors(g, node_id)
    
    ev_mat <- matrix(0, ncol = length(neighbors_list)) 
    colnames(ev_mat) <- as.vector(neighbors_list) 
    rownames(ev_mat) <- node_id
    
    for (neighbor_id in neighbors_list){
      ev_mat[as.character(node_id), as.character(neighbor_id)] <- EdgeIntensity(obj = obj, 
                                                                                node_id1 = igraph::V(g)[node_id]$name, 
                                                                                node_id2 = igraph::V(g)[neighbor_id]$name)
    }
    
    mean_intensity <- Reduce('+', ev_mat) / igraph::degree(g, node_id)
    
    mean_intensity
  }
}


#' Calculates edgewise and mean nodewise intensities for for Undirected networks and, for each edge, the proportions of
#' all event covariates.
#' 
#' @name RelateEventsToNetwork.intensitynetUnd
#' 
#' @param obj intensitynetUnd object
#' 
#' @return proper intensitynetUnd object with a graph containing the nodewise intensity in the node 
#' attributes and the edgewise intensities and event covariate proportions as edge attributes.
#' 
#' @export
RelateEventsToNetwork.intensitynetUnd = function(obj){
  g <- obj$graph
  counts <- c()
  
  if(length(obj$events) == 0){
    warning("No events, cannot calculate any intensity.")
    return(obj)
  }
  
  tmp_obj <- EdgeIntensitiesAndProportions.intensitynet(obj)
  g <- tmp_obj$graph
  
  message("\nCalculating node intensities...")
  pb = utils::txtProgressBar(min = 0, max = igraph::gorder(g), initial = 0, style = 3) 
  
  # check if the intensities was previously calculated, if not, calculate them
  v_count <- 0
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
      v_count <- v_count + 1
      counts[[node_id]] <- igraph::vertex_attr(g, 'intensity', node_id)
    }
  }
  close(pb)
  # If the intensity of all edges is already calculated return the object
  if(v_count == length(igraph::V(g))){
    message("Warning: Intensities were already calculated in a previous instance, returning the same object.")
    return(obj)
  } 
  
  
  #g <- g %>% igraph::set_vertex_attr(name = "intensity", value = as.matrix(counts))
  g <- igraph::set_vertex_attr(g, name = "intensity", value = as.matrix(counts))
  
  intnet <- list(graph = g, 
                 events = obj$events, 
                 graph_type = obj$graph_type, 
                 distances_mtx = obj$distances_mtx,
                 event_correction = obj$event_correction)
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
#' @param show_events option to show the events as orange squares, FALSE by default
#' @param path vector with the nodes of the path to be highlighted. Default NULL
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
#' plot(und_intnet_chicago, path = c("V1","V2","V24","V25","V26","V48")) # highlight a path
#' 
#' @export
plot.intensitynetUnd <- function(x, vertex_labels = 'none', edge_labels = 'none', 
                                 xy_axes = TRUE, enable_grid = FALSE, show_events = FALSE, path = NULL, ...){
  g <- x$graph
  
  if(!is.null(path) && length(path) == 1){
    stop("A path must contain more than one vertex")
  }
  
  v_label <- switch(vertex_labels, 
                    none = {''}, 
                    intensity = {round(igraph::vertex_attr(g)$intensity, 4)},
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
                      show_events = show_events,
                      path = path)
  class(geoplot_obj) <- "netTools"
  
  GeoreferencedPlot(geoplot_obj, ...)
}