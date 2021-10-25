calculateEventIntensities.intensitynetUnd = function(obj){
  g <- obj$graph
  intensities <- obj$intensities
  edge_counts <- c()
  counts <- c()
  
  pb = txtProgressBar(min = 0, max = gorder(g), initial = 0) 
  
  #TODO: Consider progressive calculation with the other methods
  
  # check if the intensities was previously calculated, if not, calculate them
  if(is.null(intensities)){
    cat("Calculating intensities...\n")
    for(node_id in V(g)){
      if(degree(g, node_id) > 0){
        
        ev_mat <- nodeIntensity(obj, node_id)$intensity
        
        #Adds result of Edgewise intenisty function to 'edge_counts'
        edge_counts[[node_id]] <- ev_mat   
        #Adds result of Nodewise mean intenisty function to 'counts'
        counts[[node_id]] <- Reduce('+', ev_mat)/degree(g, node_id)
        
        
        if(length(counts[[node_id]]) == 0){
          counts[[node_id]] <- 0
        }
      }
      
      # Counts for isolated nodes
      else{
        counts[[node_id]] <- 0
      }
      setTxtProgressBar(pb,node_id)
    }
    close(pb)
    
    intensities <- list(node_mean_intensity = counts, edge_intensity = edge_counts)
    
    g <- g %>% set_vertex_attr(name = "intensity", value = as.matrix(counts)) %>%
      set_edge_attr(name = "intensity", value = as.matrix(edge_counts))
    
    intnet <- list(graph = g, events = obj$events_mtx, graph_type = obj$graph_type, distances = obj$dist_mtx, intensities = intensities)
    attr(intnet, 'class') <- c("intensitynet", "intensitynetUnd")
    return(intnet)
  }
  else{
    print("The intensities of this net have already been calculated.")
    return(obj)
  }
}

edgeIntensity.intensitynet= function(obj,  node_id1, node_id2){
  g <- obj$graph
  distances_mtx <- obj$distances
  events_mtx <- obj$events
  
  # Note that the igraph library already handle the error when one of the node id's 
  # are not part of the graph. Also gives the proper information about it.
  edge_id <- get.edge.ids(g, c(node_id1, node_id2))
  
  # If the intensity of this edge was previously calculated, then return it
  if(edge_id != 0 & !is.null(edge_attr(g, "intensity", index=edge_id))){
    if(length(is.na(vertex_attr(g, "intensity", edge_id)))==0){
      return(edge_attr(g, 'intensity', index=edge_id))
    }
  }
  
  # Distance between the node and its neighbor
  res <- tryCatch(
    {
      abs(distances_mtx[node_id1, node_id2]) # Distance between the node and its neighbor 
    },
    # If the nodes are not part of the graph, give the proper information of the error.
    error=function(cond) {
      if(is.null(match(node_id1, V(g))) & is.null(match(node_id2, V(g)))){
        message("First and second vertices (node_id1, node_id2) doesn't exist in the graph.")
      }else if( is.null(match(node_id1, V(g))) ){
        message("First vertice ID (node_id1) doesn't exist in the graph.")
      }else if( is.null(match(node_id2, V(g))) ){
        message("Second vertice ID (v2) doesn't exist in the graph.")
      }else{
        neighbors_list <- neighbors(g, node_id1)
        if(! V(g)[node_id2] %in% neighbors_list){
          message("Second vertice (node_id2) it's not a neighbor of first vertice (node_id1)")
        }else{
          message(cond)
        }
      }
    }
  )    
  
  node1_coords <- list(xcoord = vertex_attr(g, "xcoord", node_id1),
                       ycoord = vertex_attr(g, "ycoord", node_id1))
  
  node2_coords <- list(xcoord = vertex_attr(g, "xcoord", node_id2),
                       ycoord = vertex_attr(g, "ycoord", node_id2))
  
  # Count events inside a window formed by the node and its neighbor
  
  # Defining event window
  x_min <- min(c(node1_coords$xcoord, node2_coords$xcoord))
  x_max <- max(c(node1_coords$xcoord, node2_coords$xcoord))
  
  y_min <- min(c(node1_coords$ycoord, node2_coords$ycoord))
  y_max <- max(c(node1_coords$ycoord, node2_coords$ycoord))
  
  indicator <- 0
  # Counting events
  for(row in 1:nrow(events_mtx)) {
    event_x <- events_mtx[row, 1]
    event_y <- events_mtx[row, 2]
    
    if(x_min <= event_x & x_max >= event_x & y_min <= event_y & y_max >= event_y){
      indicator <- indicator + 1
    } 
  }
  edge_intensity <- indicator/res
  
  edge_intensity
}

edgeIntensity.intensitynet= function(obj,  node_id1, node_id2, z=100){
  
  if(node_id1 == node_id2){
    stop("Both vertices cannot be the same.")
  }
  
  if(z <= 0){
    print("Warning: 'z' cannot be equal or less than 0, using default.")
    z <- 0
  }
  
  g <- obj$graph
  distances_mtx <- obj$distances
  events_mtx <- obj$events
  
  # Note that the igraph library already handle the error when one of the node id's 
  # are not part of the graph. Also gives the proper information about it.
  edge_id <- get.edge.ids(g, c(node_id1, node_id2))
  
  # If the intensity of this edge was previously calculated, then return it
  if(edge_id != 0 & !is.null(edge_attr(g, "intensity", index=edge_id))){
    if(length(is.na(vertex_attr(g, "intensity", edge_id)))==0){
      return(edge_attr(g, 'intensity', index=edge_id))
    }
  }
  
  # Distance between the node and its neighbor
  res <- tryCatch(
    {
      abs(distances_mtx[node_id1, node_id2]) # Distance between the node and its neighbor 
    },
    # If the nodes are not part of the graph, give the proper information of the error.
    error=function(cond) {
      neighbors_list <- neighbors(g, node_id1)
      if(! V(g)[node_id2] %in% neighbors_list){
        message("Second vertice (node_id2) it's not a neighbor of first vertice (node_id1)")
      }else{
        message(cond)
      }
    }
  )    
  
  node1_coords <- list(xcoord = vertex_attr(g, "xcoord", node_id1),
                       ycoord = vertex_attr(g, "ycoord", node_id1))
  
  node2_coords <- list(xcoord = vertex_attr(g, "xcoord", node_id2),
                       ycoord = vertex_attr(g, "ycoord", node_id2))
  
  # Count events inside a window formed by the node and its neighbor
  
  # Defining event window
  dx <- node1_coords$xcoord - node2_coords$xcoord
  dy <- node1_coords$ycoord - node2_coords$ycoord
  
  if(dx<dy){
    zx <- z/(sqrt(1+(dx/dy)^2))
    zy <- -(dx/dy)*zx
  }else{
    zy <- z/(sqrt(1+(dy/dx)^2))
    zx <- -(dy/dx)*zy
  }
  
  p1 <- c(node1_coords$xcoord - zx, node1_coords$ycoord - zy)
  p2 <- c(node1_coords$xcoord + zx, node1_coords$ycoord + zy)
  p3 <- c(node2_coords$xcoord - zx, node2_coords$ycoord - zy)
  p4 <- c(node2_coords$xcoord + zx, node2_coords$ycoord + zy)
  
  win_points <- list(p1, p2, p3, p4)
  anticlockwise <- orderPoints(x=c(p1[1], p2[1], p3[1], p4[1]), 
                               y=c(p1[2], p2[2], p3[2], p4[2]), clockwise = FALSE)
  
  win_points <- win_points[order(anticlockwise)]
  
  df <- data.frame(x=c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1]),
                   y=c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2]))
  obj_df <- list(df=df)
  class(obj_df) <- "netTools"
  
  if(clockwise(obj_df)){
    win <- owin(poly=list(x=rev(c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1])),
                          y=rev(c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2]))))
  }else{
    win <- owin(poly=list(x=c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1]),
                          y=c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2])))
  }
  
  
  indicator <- 0
  # Counting events
  for(row in 1:nrow(events_mtx)) {
    event_x <- events_mtx[row, 1]
    event_y <- events_mtx[row, 2]
    
    if(inside.owin(event_x, event_y, win)){
      indicator <- indicator + 1
    } 
  }
  edge_intensity <- indicator/res
  
  edge_intensity
}


GeoreferencedPlot.netTools = function(obj){
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  
  if(!is.null(distances_mtx)){
    norm_coords = layout.norm(matrix(cbind(vertex_attr(g)$xcoord, vertex_attr(g)$ycoord), ncol=2))
    plot(g,
         layout = norm_coords,
         vertex.label=NA,
         vertex.size=2,
         window=TRUE,
         axes=TRUE,
         edge.label = edge_attr(g)$intensity,
         edge.label.cex = 0.5)
  }
  else{
    plot(g, 
         vertex.label=NA, 
         vertex.size=2,
         vertex.size2=2)
  }
}
