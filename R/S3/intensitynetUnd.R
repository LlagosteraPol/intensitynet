#' Calculates the mean intensity of the given node (intensity of all the edges of the node / number of edges of the node)
#' 
#' @name nodeIntensity.intensitynetUnd
#' 
#' @param obj intensitynetUnd object
#' @param node_id ID of the node
#' 
#' @return mean intensity of the given node
#' 
#TODO: Set function as non-visible
MeanNodeIntensity.intensitynetUnd = function(obj, node_id){
  g <- obj$graph
  
  # If the intensity is already calculated, return it
  if(!is.null(vertex_attr(g, 'intensity', index=node_id))){
    if(length(is.na(vertex_attr(g, "intensity", index=node_id)))==0){
      return(vertex_attr(g, 'intensity', index=node_id))
    } 
  }
  
  if(igraph::degree(g, node_id) > 0){
    neighbors_list <- neighbors(g, node_id)
    
    ev_mat <- matrix(0, ncol = length(neighbors_list)) 
    colnames(ev_mat) <- as.vector(neighbors_list) 
    rownames(ev_mat) <- node_id
    
    for (neighbor_id in neighbors_list){
      ev_mat[as.character(node_id), as.character(neighbor_id)] <- EdgeIntensity(obj, V(g)[node_id]$name
                                                                                   , V(g)[neighbor_id]$name)
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
CalculateEventIntensities.intensitynetUnd = function(obj){
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
      edge_counts[[edge_id]] <- EdgeIntensity(obj, ends(g, edge_id)[1], ends(g, edge_id)[2])
    }else if(length(is.na(edge_attr(g, 'intensity', edge_id)))!=0){
      edge_counts[[edge_id]] <- 0
    }else{
      edge_counts[[edge_id]] <- edge_attr(g, 'intensity', edge_id)
    }
  }
  close(pb)
  
  # g <- g %>% set_edge_attr(name = "intensity", value = as.matrix(edge_counts))
  # 
  # # Encapsulate Edge intensities to pass them to 'MeanNodeIntensity' function to prevent its re-calculation
  # tmp_obj <- list(graph = g, events = obj$events, graph_type = obj$graph_type, distances_mtx = obj$distances_mtx)
  # attr(tmp_obj, 'class') <- c("intensitynet", "intensitynetUnd")
  
  # Encapsulate Edge intensities to pass them to 'MeanNodeIntensity' function to prevent its re-calculation
  tmp_obj <- SetNetworkAttribute(obj = obj, where = 'edge', name = 'intensity', value = as.matrix(edge_counts))
  g <- tmp_obj$graph
  
  pb = txtProgressBar(min = 0, max = gorder(g), initial = 0) 
  cat("Calculating node intensities...\n")
  # check if the intensities was previously calculated, if not, calculate them
  for(node_id in V(g)){
    setTxtProgressBar(pb,node_id)
    
    if(is.null(vertex_attr(g, 'intensity', node_id))){
      if(igraph::degree(g, node_id) > 0){
        #Adds result of Nodewise mean intensity function to 'counts'
        counts[[node_id]]  <- MeanNodeIntensity(tmp_obj, node_id)
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

  intnet <- list(graph = g, events = obj$events, graph_type = obj$graph_type, distances_mtx = obj$distances_mtx)
  attr(intnet, 'class') <- c("intensitynet", "intensitynetUnd")
  return(intnet)
}


#' Plot intensitynet object
#'
#' @name plot.intensitynetUnd
#'
#' @param obj intensitynet object
#' 
plot.intensitynetUnd <- function(obj, vertex_intensity='none', edge_intensity='none', xy_axes=TRUE, enable_grid=FALSE, ...){
  g <- obj$graph
  
  v_label <- switch(vertex_intensity, 
                    none = {''}, 
                    intensity = {round(vertex_attr(g)$intensity, 4)},
                    '')
  
  e_label <- switch(edge_intensity, 
                    none = {''}, 
                    intensity = {round(edge_attr(g)$intensity, 4)},
                    '')
  
  geoplot_obj <- list(graph=g, distances_mtx = obj$distances_mtx)
  class(geoplot_obj) <- "netTools"
  
  GeoreferencedPlot(geoplot_obj, vertex_intensity=v_label, edge_intensity=e_label, xy_axes=xy_axes, enable_grid=enable_grid, ...)
}