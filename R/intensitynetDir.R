#' Given a node, calculates its mean intensities regarding in and out edges edges associated with the node.
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
  if(!is.null(vertex_attr(g, 'intensity_in', index=node_id)) &
     !is.null(vertex_attr(g, 'intensity_out', index=node_id))){
    if(length(is.na(vertex_attr(g, "intensity_in", index=node_id)))==0 &
       length(is.na(vertex_attr(g, "intensity_out", index=node_id)))==0){
        return( list(in_int  = vertex_attr(g, 'intensity_in', index=node_id),
                     out_int = vertex_attr(g, 'intensity_out', index=node_id)))
    }
  }
  
  if(igraph::degree(g, node_id) > 0){
    in_neighbors  <- neighbors(g, node_id, mode = 'in')
    out_neighbors <- neighbors(g, node_id, mode = 'out')
    
    if(length(in_neighbors) > 0){
      in_mat <- matrix(0, ncol = length(in_neighbors)) 
      colnames(in_mat) <- as.vector(in_neighbors) 
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
      colnames(out_mat) <- as.vector(out_neighbors) 
      rownames(out_mat) <- node_id
      
      for (neighbor_id in out_neighbors){
        out_mat[as.character(node_id), as.character(neighbor_id)] <- EdgeIntensity(obj, V(g)[node_id]$name
                                                                                   , V(g)[neighbor_id]$name)
      }
      
      out_intensity <- Reduce('+', out_mat)/length(out_neighbors)
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
CalculateEventIntensities.intensitynetDir = function(obj){
  g <- obj$graph
  intensities <- obj$intensities
  edge_counts <- c()
  in_counts <- c()
  out_counts <- c()
  
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
        in_counts[[node_id]]  <- intensities$in_int
        out_counts[[node_id]] <- intensities$out_int
      }else{
        # Counts for isolated nodes or NA values
        in_counts[[node_id]]  <- 0
        out_counts[[node_id]] <- 0
      }
    }else if(length(is.na(vertex_attr(g, 'intensity_in', node_id)))!=0 ||
             length(is.na(vertex_attr(g, 'intensity_out', node_id)))!=0){
      if(length(is.na(vertex_attr(g, 'intensity_in', node_id)))!=0) in_counts[[node_id]]   <- 0
      if(length(is.na(vertex_attr(g, 'intensity_out', node_id)))!=0) out_counts[[node_id]] <- 0
    }else{
      in_counts[[node_id]]  <- vertex_attr(g, 'intensity_in', node_id)
      out_counts[[node_id]] <- vertex_attr(g, 'intensity_out', node_id)
    }
  }
  close(pb)
  
  g <- g %>% set_vertex_attr(name = "intensity_in", value = as.matrix(in_counts)) %>% 
             set_vertex_attr(name = "intensity_out", value = as.matrix(out_counts))
  
  intnet <- list(graph = g, events = obj$events, graph_type = obj$graph_type, distances_mtx = obj$distances_mtx)
  attr(intnet, 'class') <- c("intensitynet", "intensitynetDir")
  return(intnet)
}


#' Plot intensitynet object
#'
#' @name plot.intensitynetDir
#'
#' @param obj intensitynet object
#' 
plot.intensitynetDir <- function(obj, vertex_intensity='none', edge_intensity='none', xy_axes=TRUE, enable_grid=FALSE, ...){
  g <- obj$graph
  
  v_label <- switch(vertex_intensity, 
                    none = {''}, 
                    intensity_in = {round(vertex_attr(g)$intensity_in, 4)},
                    intensity_out = {round(vertex_attr(g)$intensity_out, 4)},
                    '')
  
  e_label <- switch(edge_intensity, 
                    none = {''}, 
                    intensity = {round(edge_attr(g)$intensity, 4)},
                    '')
  
  geoplot_obj <- list(graph=g, distances_mtx = obj$distances_mtx)
  class(geoplot_obj) <- "netTools"
  
  GeoreferencedPlot(geoplot_obj, 
                    vertex_intensity = v_label, 
                    edge_intensity = e_label, 
                    xy_axes = xy_axes, 
                    enable_grid = enable_grid, 
                    ...)
}