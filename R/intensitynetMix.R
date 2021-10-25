#' Given a node, calculates its mean intensites depending of the edges associated with the node, those intensities are: 
#' in, out (for directed edges), undirected and total intensity.
#'
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
  if(!is.null(vertex_attr(g, 'intensity_und', index=node_id)) &
     !is.null(vertex_attr(g, 'intensity_in', index=node_id)) &
     !is.null(vertex_attr(g, 'intensity_out', index=node_id)) &
     !is.null(vertex_attr(g, 'intensity_all', index=node_id))){
    if(length(is.na(vertex_attr(g, "intensity_und", index=node_id)))==0 &
       length(is.na(vertex_attr(g, "intensity_in", index=node_id)))==0 &
       length(is.na(vertex_attr(g, "intensity_out", index=node_id)))==0 &
       length(is.na(vertex_attr(g, "intensity_all", index=node_id)))==0){
      return( list(und_int  = vertex_attr(g, 'intensity_und', index=node_id),
                   in_int   = vertex_attr(g, 'intensity_in', index=node_id),
                   out_int  = vertex_attr(g, 'intensity_out', index=node_id),
                   all_int  = vertex_attr(g, 'intensity_all', index=node_id)))
    }
  }
  
  if(igraph::degree(g, node_id) > 0){
    in_neighbors_tmp  <- as.vector(neighbors(g, node_id, mode = 'in'))
    out_neighbors_tmp <- as.vector(neighbors(g, node_id, mode = 'out'))
    
    und_neighbors <- intersect(in_neighbors_tmp, out_neighbors_tmp)
    in_neighbors  <- setdiff(in_neighbors_tmp, out_neighbors_tmp) 
    out_neighbors  <- setdiff(out_neighbors_tmp, in_neighbors_tmp) 
    
    if(length(und_neighbors) > 0){
      und_mat <- matrix(0, ncol = length(und_neighbors)) 
      colnames(und_mat) <- und_neighbors
      rownames(und_mat) <- node_id
      
      for (neighbor_id in und_neighbors){
        und_mat[as.character(node_id), as.character(neighbor_id)] <- EdgeIntensity(obj, V(g)[node_id]$name
                                                                                   , V(g)[neighbor_id]$name)
      }
      und_intensity <- Reduce('+', und_mat)/length(und_neighbors)
    }else{
      und_intensity <- 0
    }
    
    if(length(in_neighbors) > 0){
      in_mat <- matrix(0, ncol = length(in_neighbors)) 
      colnames(in_mat) <- in_neighbors
      rownames(in_mat) <- node_id
      
      for (neighbor_id in in_neighbors){
        in_mat[as.character(node_id), as.character(neighbor_id)] <- EdgeIntensity(obj, V(g)[node_id]$name
                                                                                  , V(g)[neighbor_id]$name)
      }
      in_intensity <- Reduce('+', in_mat)/length(in_neighbors)
    }else{
      in_intensity <- 0
    }
      
    if(length(out_neighbors) > 0){
      out_mat <- matrix(0, ncol = length(out_neighbors)) 
      colnames(out_mat) <- out_neighbors 
      rownames(out_mat) <- node_id
      
      for (neighbor_id in out_neighbors){
        out_mat[as.character(node_id), as.character(neighbor_id)] <- EdgeIntensity(obj, V(g)[node_id]$name
                                                                                   , V(g)[neighbor_id]$name)
      }
      out_intensity <- Reduce('+', out_mat)/length(out_neighbors)
    }else{
      out_intensity <- 0
    }
 
    all_intensity <- und_intensity + in_intensity + out_intensity / 
                    (length(und_neighbors) + length(in_neighbors) + length(out_neighbors))
    
    list(und_int = und_intensity, 
         in_int  = in_intensity, 
         out_int = out_intensity,
         all_int = all_intensity)
  }
}


#' Calculates edgewise and mean nodewise intensities for Mixed networks
#' 
#' @name CalculateEventIntensities.intensitynetMix
#' 
#' @param obj intensitynetMix object
#' 
#' @return intensitynetMix object with a graph containing all the intensities as attributes of its nodes and edges
#' 
CalculateEventIntensities.intensitynetMix = function(obj){
  g <- obj$graph
  intensities <- obj$intensities
  edge_counts <- c()
  und_counts <- c()
  in_counts <- c()
  out_counts <- c()
  all_counts <- c()
  
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
  
  # Encapsulate Edge intensities to pass them to 'MeanNodeIntensity' function to prevent its re-calculation
  tmp_obj <- SetNetworkAttribute(obj = obj, where = 'edge', name = 'intensity', value = as.matrix(edge_counts))
  g <- tmp_obj$graph
  
  pb = txtProgressBar(min = 0, max = gorder(g), initial = 0) 
  cat("Calculating node intensities...\n")
  # check if the intensities was previously calculated, if not, calculate them
  for(node_id in V(g)){
    setTxtProgressBar(pb,node_id)
    
    if(is.null(vertex_attr(g, 'intensity_in', node_id)) || is.null(vertex_attr(g, 'intensity_out', node_id))){
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
    }else if(length(is.na(vertex_attr(g, 'intensity_und', node_id)))!=0 ||
             length(is.na(vertex_attr(g, 'intensity_in', node_id)))!=0 ||
             length(is.na(vertex_attr(g, 'intensity_out', node_id)))!=0 ||
             length(is.na(vertex_attr(g, 'intensity_all', node_id)))!=0){
      if(length(is.na(vertex_attr(g, 'intensity_und', node_id)))!=0) in_counts[[node_id]]   <- 0
      if(length(is.na(vertex_attr(g, 'intensity_in', node_id)))!=0) out_counts[[node_id]] <- 0
      if(length(is.na(vertex_attr(g, 'intensity_out', node_id)))!=0) out_counts[[node_id]] <- 0
      if(length(is.na(vertex_attr(g, 'intensity_all', node_id)))!=0) out_counts[[node_id]] <- 0
    }else{
      und_counts[[node_id]] <- vertex_attr(g, 'intensity_und', node_id)
      in_counts[[node_id]]  <- vertex_attr(g, 'intensity_in', node_id)
      out_counts[[node_id]] <- vertex_attr(g, 'intensity_out', node_id)
      in_counts[[node_id]]  <- vertex_attr(g, 'intensity_all', node_id)
    }
  }
  close(pb)
  
  g <- g %>% set_vertex_attr(name = "intensity_und", value = as.matrix(und_counts)) %>% 
             set_vertex_attr(name = "intensity_in", value = as.matrix(in_counts))   %>% 
             set_vertex_attr(name = "intensity_out", value = as.matrix(out_counts))  %>% 
             set_vertex_attr(name = "intensity_all", value = as.matrix(all_counts))
  
  intnet <- list(graph = g, events = obj$events, graph_type = obj$graph_type, distances_mtx = obj$distances_mtx)
  attr(intnet, 'class') <- c("intensitynet", "intensitynetMix")
  return(intnet)
}

#' Plot intensitynet object
#'
#' @name plot.intensitynetMix
#'
#' @param obj intensitynet object
#' 
plot.intensitynetMix <- function(obj, vertex_intensity='none', edge_intensity='none', xy_axes=TRUE, enable_grid=FALSE, ...){
  g <- obj$graph
  
  v_label <- switch(vertex_intensity, 
                    none = {''}, 
                    intensity_und = {round(vertex_attr(g)$intensity_und, 4)},
                    intensity_in = {round(vertex_attr(g)$intensity_in, 4)},
                    intensity_out = {round(vertex_attr(g)$intensity_out, 4)},
                    intensity_all = {round(vertex_attr(g)$intensity_all, 4)},
                    '')
  
  e_label <- switch(edge_intensity, 
                    none = {''}, 
                    intensity = {round(edge_attr(g)$intensity, 4)},
                    '')
  
  geoplot_obj <- list(graph=g, distances_mtx = obj$distances_mtx)
  class(geoplot_obj) <- "netTools"
  
  GeoreferencedPlot(geoplot_obj, vertex_intensity=v_label, edge_intensity=e_label, xy_axes=xy_axes, enable_grid=enable_grid, ...)
}