#' Calculates the mean intensity of the given node (for directed networks)
#' 
#' @description 
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
      
      message("Warning: Node intensities were already calculated in a previous instance, returning the same intensity.")
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
                                                           node_id2 = igraph::V(g)[neighbor_id]$name)
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
                                                            node_id2 = igraph::V(g)[neighbor_id]$name)
      }
      
      out_intensity <- Reduce('+', out_mat) / length(out_neighbors)
    }else{
      out_intensity <- 0
    }
    
    list(in_int = in_intensity, out_int = out_intensity)
  }
}


#' Calculates intensity statistics for the given intensitynet object
#' 
#' @description 
#' Calculates edgewise and mean nodewise intensities for Directed networks and, for each edge, the proportions of
#' all event covariates.
#' 
#' @name RelateEventsToNetwork.intensitynetDir
#' 
#' @param obj intensitynetDir object
#' 
#' @return proper intensitynetDir object with a graph containing the nodewise intensity in the node 
#' attributes and the edgewise intensities and event covariate proportions as edge attributes.
#' 
#' @export
RelateEventsToNetwork.intensitynetDir = function(obj){
  g <- obj$graph
  intensities <- obj$intensities
  in_counts <- c()
  out_counts <- c()
  
  if(length(obj$events) == 0){
    stop("Error: No events, cannot calculate any intensity.")
  }
  
  tmp_obj <- EdgeIntensitiesAndProportions.intensitynet(obj)
  g <- tmp_obj$graph
  
  message("\nCalculating node intensities...")
  pb = utils::txtProgressBar(min = 0, max = igraph::gorder(g), initial = 0, style = 3) 
  
  # check if the intensities was previously calculated, if not, calculate them
  v_count <- 0
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
      v_count <- v_count + 1
      in_counts[[node_id]]  <- igraph::vertex_attr(g, 'intensity_in', node_id)
      out_counts[[node_id]] <- igraph::vertex_attr(g, 'intensity_out', node_id)
    }
  }
  close(pb)
  # If the intensity of all edges is already calculated return the object
  if(v_count == length(igraph::V(g))){
    message("Warning: Intensities were already calculated in a previous instance, returning the same object.")
    return(obj)
  } 
  
  # g <- g %>% igraph::set_vertex_attr(name = "intensity_in", value = as.matrix(in_counts)) %>% 
  #            igraph::set_vertex_attr(name = "intensity_out", value = as.matrix(out_counts))
  g <- igraph::set_vertex_attr(g, name = "intensity_in", value = as.numeric(in_counts))
  g <- igraph::set_vertex_attr(g, name = "intensity_out", value = as.numeric(out_counts))
  
  intnet <- list(graph = g, 
                 events = obj$events, 
                 graph_type = obj$graph_type, 
                 distances_mtx = obj$distances_mtx,
                 event_correction = obj$event_correction,
                 events_related = TRUE)
  attr(intnet, 'class') <- c("intensitynetDir", "intensitynet")
  return(intnet)
}