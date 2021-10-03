#' Calculates the intensity of the given node and input into the node attribute of the graph.
#' 
#' @name nodeIntensity.intensitynetUnd
#' 
#' @param node_id ID of the node
#' 
#' @return list with the total intensity of the node ('intensity') and 
#' the intensity respect each neighbor ('detailed') 
#' 
#TODO: Set function as non-visible
meanNodeIntensity.intensitynetUnd = function(obj, node_id){
  g <- obj$graph
  
  # If the intensity is already calculated, return it
  if(!is.null(vertex_attr(g, 'intensity', index=node_id))){
    if(length(is.na(vertex_attr(g, "intensity", index=node_id)))==0){
      return(vertex_attr(g, 'intensity', index=node_id))
    } 
  }
  
  if(degree(g, node_id) > 0){
    neighbors_list <- neighbors(g, node_id)
    
    ev_mat <- matrix(0, ncol = length(neighbors_list)) 
    colnames(ev_mat) <- as.vector(neighbors_list) 
    rownames(ev_mat) <- node_id
    
    for (neighbor_id in neighbors_list){
      ev_mat[as.character(node_id), as.character(neighbor_id)] <- edgeIntensity(obj, V(g)[node_id]$name
                                                                                   , V(g)[neighbor_id]$name)
    }
    
    total_intensity <- Reduce('+', ev_mat)
    
    mean_intensity <- total_intensity/degree(g, node_id)
    mean_intensity
  }
}


#' Calculates edgewise and mean nodewise intensity function for Undirected networks
#' 
#' @name calculateEventIntensities.intensitynetUnd
#' 
#' @return 
#' 
calculateEventIntensities.intensitynetUnd = function(obj){
  g <- obj$graph
  intensities <- obj$intensities
  edge_counts <- c()
  counts <- c()
  
  pb = txtProgressBar(min = 0, max = gsize(g), initial = 0) 
  cat("Calculating edge intensities...\n")
  for(edge_id in E(g)){
    setTxtProgressBar(pb,edge_id)
    if(is.null(edge_attr(g, 'intensity', edge_id))){
      #Adds result of Edgewise intenisty function to 'edge_counts'
      edge_counts[[edge_id]] <- edgeIntensity(obj, ends(g, edge_id)[1], ends(g, edge_id)[2])
    }else if(length(is.na(edge_attr(g, 'intensity', edge_id)))!=0){
      edge_counts[[edge_id]] <- 0
    }else{
      edge_counts[[edge_id]] <- edge_attr(g, 'intensity', edge_id)
    }
  }
  close(pb)
  
  g <- g %>% set_edge_attr(name = "intensity", value = as.matrix(edge_counts))
  
  # Encapsulate Edge intensities to pass them to 'meanNodeIntensity' function to prevent its re-calculation
  tmp_obj <- list(graph = g, events = obj$events, graph_type = obj$graph_type, distances = obj$distances)
  attr(tmp_obj, 'class') <- c("intensitynet", "intensitynetUnd")
  
  pb = txtProgressBar(min = 0, max = gorder(g), initial = 0) 
  cat("Calculating node intensities...\n")
  # check if the intensities was previously calculated, if not, calculate them
  for(node_id in V(g)){
    setTxtProgressBar(pb,node_id)
    
    if(is.null(vertex_attr(g, 'intensity', node_id))){
      if(degree(g, node_id) > 0){
        #Adds result of Nodewise mean intenisty function to 'counts'
        counts[[node_id]]  <- meanNodeIntensity(tmp_obj, node_id)
      }else{
        # Counts for isolated nodes or NA values
        counts[[node_id]] <- 0
      }
    }else if(length(is.na(vertex_attr(g, 'intensity', node_id)))!=0){
      counts[[node_id]] <- 0
    }else{
      counts[[node_id]] <- vertex_attr(g, 'intensity', node_id)
    }
  }
  close(pb)
    
  g <- g %>% set_vertex_attr(name = "intensity", value = as.matrix(counts))

  intnet <- list(graph = g, events = obj$events, graph_type = obj$graph_type, distances = distances)
  attr(intnet, 'class') <- c("intensitynet", "intensitynetUnd")
  return(intnet)
}