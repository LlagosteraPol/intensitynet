#' Given a node, calculates its mean intensities depending on the edges associated with the node, those intensities are: 
#' in, out (for directed edges), undirected and total intensity.
#' 
#' @name nodeIntensity.intensitynetMix
#' 
#' @param obj intensitynetMix object
#' @param node_id ID of the node
#' 
#' @return mean intensities of the given node for undirected edges, in and out directed and total intensity.
#' 
MeanNodeIntensity.intensitynetMix = function(obj, node_id){
  g <- obj$graph
  
  # If the intensities are already calculated, return them
  if(!is.null(igraph::vertex_attr(g, 'intensity_und', index = node_id)) &
     !is.null(igraph::vertex_attr(g, 'intensity_in', index = node_id))  &
     !is.null(igraph::vertex_attr(g, 'intensity_out', index = node_id)) &
     !is.null(igraph::vertex_attr(g, 'intensity_all', index = node_id))) {
    
    if(!is.na(igraph::vertex_attr(g, "intensity_und", index = node_id))[1] &
       !is.na(igraph::vertex_attr(g, "intensity_in", index = node_id))[1]  &
       !is.na(igraph::vertex_attr(g, "intensity_out", index = node_id))[1] &
       !is.na(igraph::vertex_attr(g, "intensity_all", index = node_id))[1]) {
      
      message("Warning: Node intensities were already calculated in a previous instance, returning the same intensity.")
      return( list( und_int  = igraph::vertex_attr(g, 'intensity_und', index=node_id),
                    in_int   = igraph::vertex_attr(g, 'intensity_in', index=node_id),
                    out_int  = igraph::vertex_attr(g, 'intensity_out', index=node_id),
                    all_int  = igraph::vertex_attr(g, 'intensity_all', index=node_id) ) )
    }
  }
  
  if(igraph::degree(g, node_id) > 0){
    in_neighbors_tmp  <- as.vector( igraph::neighbors(g, node_id, mode = 'in') )
    out_neighbors_tmp <- as.vector( igraph::neighbors(g, node_id, mode = 'out') )
    
    und_neighbors <- intersect(in_neighbors_tmp, out_neighbors_tmp)
    in_neighbors  <- setdiff(in_neighbors_tmp, out_neighbors_tmp) 
    out_neighbors  <- setdiff(out_neighbors_tmp, in_neighbors_tmp) 
    
    if(length(und_neighbors) > 0){
      und_mat <- matrix(0, ncol = length(und_neighbors)) 
      colnames(und_mat) <- und_neighbors
      rownames(und_mat) <- node_id
      
      for (neighbor_id in und_neighbors){
        und_mat[as.character(node_id), as.character(neighbor_id)] <- EdgeIntensity(obj = obj, 
                                                                                   node_id1 = igraph::V(g)[node_id]$name, 
                                                                                   node_id2 = igraph::V(g)[neighbor_id]$name)
      }
      und_intensity <- Reduce('+', und_mat) / length(und_neighbors)
    }else{
      und_intensity <- 0
    }
    
    if(length(in_neighbors) > 0){
      in_mat <- matrix(0, ncol = length(in_neighbors)) 
      colnames(in_mat) <- in_neighbors
      rownames(in_mat) <- node_id
      
      for (neighbor_id in in_neighbors){
        in_mat[as.character(node_id), as.character(neighbor_id)] <- EdgeIntensity(obj = obj, 
                                                                                  node_id1 = igraph::V(g)[node_id]$name, 
                                                                                  node_id2 = igraph::V(g)[neighbor_id]$name)
      }
      in_intensity <- Reduce('+', in_mat) / length(in_neighbors)
    }else{
      in_intensity <- 0
    }
    
    if(length(out_neighbors) > 0){
      out_mat <- matrix(0, ncol = length(out_neighbors)) 
      colnames(out_mat) <- out_neighbors 
      rownames(out_mat) <- node_id
      
      for (neighbor_id in out_neighbors){
        out_mat[as.character(node_id), as.character(neighbor_id)] <- EdgeIntensity(obj = obj, 
                                                                                   node_id1 = igraph::V(g)[node_id]$name, 
                                                                                   node_id2 = igraph::V(g)[neighbor_id]$name)
      }
      out_intensity <- Reduce('+', out_mat) / length(out_neighbors)
    }else{
      out_intensity <- 0
    }
    
    all_intensity <- (und_intensity + in_intensity + out_intensity) / 
      (length(und_neighbors) + length(in_neighbors) + length(out_neighbors))
    
    list(und_int = und_intensity, 
         in_int  = in_intensity, 
         out_int = out_intensity,
         all_int = all_intensity)
  }
}


#' Calculates edgewise and mean nodewise intensities for Mixed networks and, for each edge, the proportions of
#' all event covariates.
#' 
#' @name RelateEventsToNetwork.intensitynetMix
#' 
#' @param obj intensitynetMix object
#' 
#' @return proper intensitynetMix object with a graph containing the nodewise intensity in the node 
#' attributes and the edgewise intensities and event covariate proportions as edge attributes.
#' 
#' @export
RelateEventsToNetwork.intensitynetMix = function(obj){
  g <- obj$graph
  intensities <- obj$intensities
  und_counts <- c()
  in_counts <- c()
  out_counts <- c()
  all_counts <- c()
  
  if(length(obj$events) == 0){
    warning("No events, cannot calculate any intensity.")
    return(obj)
  }
  
  tmp_obj <- EdgeIntensitiesAndProportions.intensitynet(obj)
  g <- tmp_obj$graph
  
  message("Calculating node intensities...")
  pb = utils::txtProgressBar(min = 0, max = igraph::gorder(g), initial = 0, style = 3) 
  
  # check if the intensities was previously calculated, if not, calculate them
  v_count <- 0
  for(node_id in igraph::V(g)){
    utils::setTxtProgressBar(pb,node_id)
    
    if(is.null(igraph::vertex_attr(g, 'intensity_in', node_id)) || is.null(igraph::vertex_attr(g, 'intensity_out', node_id))){
      if(igraph::degree(g, node_id) > 0){
        #Adds result of Nodewise mean intensity function to 'counts'
        intensities <- MeanNodeIntensity(tmp_obj, node_id)
        und_counts[[node_id]]  <- intensities$und_int
        in_counts[[node_id]]  <- intensities$in_int
        out_counts[[node_id]] <- intensities$out_int
        all_counts[[node_id]] <- intensities$all_int
      }else{
        # Counts for isolated nodes or NA values
        und_counts[[node_id]]  <- 0
        in_counts[[node_id]]   <- 0
        out_counts[[node_id]]  <- 0
        all_counts[[node_id]]  <- 0
      }
    }else if(is.na(igraph::vertex_attr(g, 'intensity_und', node_id))[1] ||
             is.na(igraph::vertex_attr(g, 'intensity_in', node_id))[1]  ||
             is.na(igraph::vertex_attr(g, 'intensity_out', node_id))[1] ||
             is.na(igraph::vertex_attr(g, 'intensity_all', node_id))[1]){
      
      if(is.na(igraph::vertex_attr(g, 'intensity_und', node_id))[1]) in_counts[[node_id]] <- 0
      if(is.na(igraph::vertex_attr(g, 'intensity_in', node_id))[1]) out_counts[[node_id]] <- 0
      if(is.na(igraph::vertex_attr(g, 'intensity_out', node_id))[1]) out_counts[[node_id]] <- 0
      if(is.na(igraph::vertex_attr(g, 'intensity_all', node_id))[1]) out_counts[[node_id]] <- 0
      
    }else{
      v_count <- v_count + 1
      und_counts[[node_id]] <- igraph::vertex_attr(g, 'intensity_und', node_id)
      in_counts[[node_id]]  <- igraph::vertex_attr(g, 'intensity_in', node_id)
      out_counts[[node_id]] <- igraph::vertex_attr(g, 'intensity_out', node_id)
      in_counts[[node_id]]  <- igraph::vertex_attr(g, 'intensity_all', node_id)
    }
  }
  close(pb)
  # If the intensity of all edges is already calculated return the object
  if(v_count == length(igraph::V(g))){
    message("Warning: Intensities were already calculated in a previous instance, returning the same object.")
    return(obj)
  } 
  
  # g <- g %>% igraph::set_igraph::vertex_attr(name = "intensity_und", value = as.matrix(und_counts)) %>% 
  #            igraph::set_vertex_attr(name = "intensity_in", value = as.matrix(in_counts))   %>% 
  #            igraph::set_vertex_attr(name = "intensity_out", value = as.matrix(out_counts)) %>% 
  #            igraph::set_vertex_attr(name = "intensity_all", value = as.matrix(all_counts))
  g <- igraph::set_vertex_attr(g, name = "intensity_und", value = as.matrix(und_counts))
  g <- igraph::set_vertex_attr(g, name = "intensity_in", value = as.matrix(in_counts)) 
  g <- igraph::set_vertex_attr(g, name = "intensity_out", value = as.matrix(out_counts))
  g <- igraph::set_vertex_attr(g, name = "intensity_all", value = as.matrix(all_counts))
  
  intnet <- list(graph = g, 
                 events = obj$events, 
                 graph_type = obj$graph_type, 
                 distances_mtx = obj$distances_mtx,
                 event_correction = obj$event_correction)
  attr(intnet, 'class') <- c("intensitynet", "intensitynetMix")
  return(intnet)
}


#' Plot intensitynet object
#'
#' @name plot.intensitynetMix
#'
#' @param x intensitynet object
#' @param vertex_labels list -> labels for the vertices
#' @param edge_labels list -> labels for the edges
#' @param xy_axes show the x and y axes
#' @param enable_grid draw a background grid
#' @param path vector with the nodes of the path to be highlighted. Default NULL
#' @param show_events option to show the events as orange squares, FALSE by default
#' @param alpha optional argument to set the transparency of the events (show_events = TRUE). The range is from 0.1 (transparent) to 1 (opaque). Default: alpha = 1
#' @param ... extra arguments for the plot
#' 
#' @return No return value, same as graphics::plot.
#' 
#' @examples
#' 
#' data("mix_intnet_chicago")
#' plot(mix_intnet_chicago) # basic plot
#' plot(mix_intnet_chicago, enable_grid = TRUE) # with grid
#' plot(mix_intnet_chicago, xy_axes = FALSE) # without axes
#' plot(mix_intnet_chicago, path = c("V1","V2","V24","V25","V26","V48")) # highlight a path
#' 
#' @export
plot.intensitynetMix <- function(x, vertex_labels='none', edge_labels='none', 
                                 xy_axes=TRUE, enable_grid=FALSE, show_events = FALSE, path = NULL, alpha = 1, ...){
  g <- x$graph
  
  if(!is.null(path) && length(path) == 1){
    stop("A path must contain more than one vertex")
  }
  
  v_label <- switch(vertex_labels, 
                    none = {''}, 
                    intensity_und = {round(igraph::vertex_attr(g)$intensity_und, 4)},
                    intensity_in = {round(igraph::vertex_attr(g)$intensity_in, 4)},
                    intensity_out = {round(igraph::vertex_attr(g)$intensity_out, 4)},
                    intensity_all = {round(igraph::vertex_attr(g)$intensity_all, 4)},
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
                      path = path,
                      alpha = alpha)
  class(geoplot_obj) <- "netTools"
  
  GeoreferencedPlot(geoplot_obj, ...)
}