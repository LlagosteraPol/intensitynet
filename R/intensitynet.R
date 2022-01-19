

#' Constructor of the class intensitynet. In order to create an intensitynet object, it is needed; an adjacency matrix, the
#' coordinates of the nodes and the coordinates of the events.
#'
#' @name intensitynet
#'
#' @param adjacency_mtx Network adjacency matrix
#' @param node_coords Nodes latitude and longitude matrix
#' @param event_coords Events latitude and longitude matrix
#' @param graph_type Network type: 'undirected' (default), 'directed' or 'mixed' 
#'
#' @return intensitynet class object containing: graph = <igraph>, events = <matrix>, graph_type = c('directed', 'undirected', 'mixed'), 
#' distances = <matrix>
#' 
#' @examples
#'
#' library(spatstat)
#' data(chicago)
#' chicago_df <- as.data.frame(chicago[["data"]]) # Get as dataframe the data from Chicago
#'
#' # Get the adjacency matrix. One way is to create an igraph object from the edge coordinates.
#' edges <- cbind(chicago[["domain"]][["from"]], chicago[["domain"]][["to"]])
#' chicago_net <- igraph::graph_from_edgelist(edges)
#'
#' # And then use the igraph function 'as_adjacency_matrix'
#' chicago_adj_mtx <- as.matrix(igraph::as_adjacency_matrix(chicago_net))
#' chicago_node_coords <- data.frame(xcoord = chicago[["domain"]][["vertices"]][["x"]], 
#'                                  ycoord = chicago[["domain"]][["vertices"]][["y"]])
#'
#' # Create a dataframe with the coordinates of the events 'assault'
#' chicago_assault <- chicago_df[chicago_df$marks == 'assault',]
#' assault_coordinates <- data.frame(xcoord = chicago_assault[,1],
#'                                   ycoord = chicago_assault[,2])
#'                                   
#' # Create the intensitynet object, in this case will be undirected 
#' intnet_chicago <- intensitynet(chicago_adj_mtx, 
#'                                node_coords = chicago_node_coords, 
#'                                event_coords = assault_coordinates)
#' 
#' @export
intensitynet <- function(adjacency_mtx, node_coords, event_coords, graph_type = 'undirected'){
  
  if (is.data.frame(adjacency_mtx)) {
    adjacency_mtx <- as.matrix(adjacency_mtx)
  }
  
  if (is.data.frame(node_coords)) {
    node_coords <- as.matrix(node_coords)
  }
  
  if (is.data.frame(event_coords)) {
    event_coords <- as.matrix(event_coords)
  }
  colnames(node_coords) <- c("xcoord", "ycoord")
  colnames(event_coords) <- c("xcoord", "ycoord")
  
  node_coords_obj <- list(node_coords = node_coords)
  class(node_coords_obj) <- "netTools"
  dist_mtx <- CalculateDistancesMtx(node_coords_obj)
  
  net_setup <- list(
    adjacency_mtx = adjacency_mtx,
    node_coords = node_coords,
    distances_mtx = dist_mtx,
    graph_type = graph_type
  )
  class(net_setup) <- "netTools"
  g <- InitGraph(net_setup)
  
  intnet <- list(graph = g, events = event_coords, graph_type = graph_type, distances_mtx = dist_mtx)
  attr(intnet, "class") <- "intensitynet"
  # Select the proper class
  switch(graph_type, 
         'undirected' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetUnd")},
         'directed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetDir")},
         'mixed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetMix")})
  
  intnet # return
}


# -------- Network functions ----------

#' Gives general node correlation of the network (choosing from normal correlation, covariance, 
#' moran-i or geary)
#' 
#' @name NodeGeneralCorrelation
#'
#' @param obj intensitynet object
#' @param dep_type 'correlation', 'covariance', moran', 'geary'. The type of 
#' dependence statistic to be computed.
#' @param lag_max Maximum geodesic lag at which to compute dependence
#' @param intensity vector containing the intensity values that the heatmaps
#' 
#' @return A vector containing the dependence statistics (ascending from order 0). 
#' 
#' @examples 
#' 
#' data("und_intnet_chicago")
#' g <- und_intnet_chicago$graph
#' gen_corr <- NodeGeneralCorrelation(und_intnet_chicago, dep_type = 'correlation', lag_max = 2, 
#'                                    intensity = igraph::vertex_attr(g)$intensity)
#' 
#' @export
NodeGeneralCorrelation <- function(obj, dep_type, lag_max, intensity){
  UseMethod("NodeGeneralCorrelation")
}


#' Gives node local moran-i or geary-c correlations
#' 
#' @name NodeLocalCorrelation
#' 
#' @source *Luc Anselin. A Local Indicator of Multivariate SpatialAssociation: Extending Geary's c, Geographical Analysis 2018; doi: https://doi.org/10.1111/gean.12164
#'
#' @param obj intensitynet object
#' @param dep_type 'moran', 'getis' or 'geary'. Type of local correlation to be computed (Moran-i, Getis-Gstar, Geary-c*),
#' default = 'moran'.
#' @param intensity vector containing the intensity values that which are used to calculate the correlation.
#' 
#' @return intensitynet class object which contains an igraph network with the selected correlation 
#' added into the vertices attributes
#' 
#' @examples 
#' 
#' data("und_intnet_chicago")
#' g <- und_intnet_chicago$graph
#' data_moran <- NodeLocalCorrelation(und_intnet_chicago, 
#'                                    dep_type = 'moran', 
#'                                    intensity = igraph::vertex_attr(g)$intensity)
#' moran_i <- data_moran$correlation
#' intnet <- data_moran$intnet
#' 
#' @export
NodeLocalCorrelation <- function(obj, dep_type = 'moran', intensity){
  UseMethod("NodeLocalCorrelation")
}


#' Plot the network and if specified, the correlation heatmap. Which could be:
#'
#' @name PlotHeatmap
#'
#' @param obj intensitynet object
#' @param heattype 'moran': Local Moran-i correlation (with 999 permutations), 'geary': Local Geary-c* 
#' correlation. The correlations will use the indicated intensity type.
#' The function also allow to only plot the intensity heatmap 'v_intensity' for vertices or 'e_intensity' for edges.
#' @param intensity_type name of the intenisty used to plot the heatmap. For undirected networks: 'intensity'. 
#' For directed networks: 'intensity_in' or 'intensity_out'. For mixed networks: 'intensity_in', 'intensity_out', 
#' 'intensity_und' or 'intensity_all'. If the intensity parameter is NULL, the function will use, if exist, 
#' the intensity (undirected) or intensity_in (directed) values from the network nodes.
#' @param net_vertices chosen vertices to plot the heatmap (or it related edges in case to plot the edge heatmap)
#' @param ... extra arguments for the class ggplot
#' 
#' @return The plot of the heatmap with class c("gg", "ggplot")
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' PlotHeatmap(und_intnet_chicago, heattype='moran')
#'
#' @export
PlotHeatmap <- function(obj, heattype = 'none', intensity_type = 'none', net_vertices = NULL, ...){
  UseMethod("PlotHeatmap")
}


#' Plot the net and the events in the neighborhood area of the given node
#' 
#' @name PlotNeighborhood
#' 
#' @param obj Intensitynet object
#' @param node_id Id of the node which the plot will be focused
#' @param ... Extra arguments for plotting
#' 
#' @return No return value, just plots the neighborhood and the events.
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' PlotNeighborhood(und_intnet_chicago, node_id = 'V300')
#' 
#' @export
PlotNeighborhood <- function(obj, node_id, ...){
  UseMethod("PlotNeighborhood")
}


#' Get the intensitynet object delimited by the given window
#' 
#' @name ApplyWindow
#' 
#' @param obj intensitynet object
#' @param x_coords vector containing the x coordinate limits of the window
#' @param y_coords vector containing the y coordinate limits of the window
#' 
#' @return intensitynet object delimited by the window (sub-part of the original)
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' sub_intnet_chicago <- ApplyWindow(und_intnet_chicago, 
#'                                   x_coords = c(300, 900), 
#'                                   y_coords = c(500, 1000))
#' 
#' @export
ApplyWindow <- function(obj, x_coords, y_coords){
  UseMethod("ApplyWindow")
}


ShortestNodeDistance <- function(obj, node_id1, node_id2){
  UseMethod("ShortestNodeDistance")
}


# -------- Intensity functions ----------

#' Calculates the intensity of the given path
#'
#' @name PathIntensity
#'
#' @param obj intensitynet object
#' @param path_nodes vector containing the node ID's of the path
#' 
#' @return intensity of the path
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' short_path <- ShortestPathIntensity(und_intnet_chicago, node_id1 = 'V1', node_id2 = 'V300')
#' PathIntensity(und_intnet_chicago, short_path$path)
#' 
#' @export
PathIntensity <- function(obj, path_nodes){
  UseMethod("PathIntensity")
}


#' Calculates the shortest path between two vertices and calculates its intensity
#'
#' @name ShortestPathIntensity
#'
#' @param obj intensitynet object
#' @param node_id1 starting node
#' @param node_id2 ending node
#' @param weighted TRUE or FALSE (default), tell if the distances must be taken into account 
#' 
#' @return intensity of the shortest path and the path vertices
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' ShortestPathIntensity(und_intnet_chicago, node_id1 = 'V1', node_id2 = 'V300')
#' 
#' @export
ShortestPathIntensity <- function(obj,  node_id1, node_id2, weighted = FALSE){
  UseMethod("ShortestPathIntensity")
}

#' Calculates edgewise and mean nodewise intensities for the given intensitynet object
#' 
#' @name CalculateEventIntensities
#' 
#' @param obj intensitynetDir object
#' 
#' @return proper intensitynet object (Undirected, Directed or Mixed) with a graph containing all the intensities as 
#' attributes of its nodes and edges
#' 
#' @examples 
#'
#' data("und_intnet_chicago")
#' intnet_chicago <- CalculateEventIntensities(und_intnet_chicago)
#' 
#' @export
CalculateEventIntensities <- function(obj){
  UseMethod("CalculateEventIntensities")
}


MeanNodeIntensity <- function(obj, node_id){
  UseMethod("MeanNodeIntensity")
}


EdgeIntensity <- function(obj, node_id1, node_id2, z){
  UseMethod("EdgeIntensity")
}


AllEdgeIntensities <- function(obj, z){
  UseMethod("AllEdgeIntensities")
}


SetNetworkAttribute <- function(obj, where, name, value){
  UseMethod("SetNetworkAttribute")
}


LaplacianGearyRepresentation <- function(obj,  intensity_type = 'none'){
  UseMethod("LaplacianGearyRepresentation")
}


#' If not calculated, calculates the intesnity of the edge with nodes; node_id1, node_id2. 
#' If the edge already contains an intensity, give it directly.
#'
#' @name EdgeIntensity.intensitynet
#' 
#' @param obj intensitynet object
#' @param node_id1 First node ID of the edge
#' @param node_id2 Second node ID of the edge
#' @param z Maximum distance between the event and the edge to consider the event part of the edge.
#' 
#' @return Intensity of the edge
#' 
EdgeIntensity.intensitynet <- function(obj,  node_id1, node_id2, z = 5){
  if(node_id1 == node_id2){
    stop("The two vertices cannot be the same.")
  }
  
  if(z <= 0){
    print("Warning: 'z' cannot be equal or less than 0, using default.")
    z <- 5
  }
  
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  event_coords <- obj$events
  
  # Note that the igraph library already handle the error when one of the node id's 
  # are not part of the graph. Also gives the proper information about it.
  edge_id <- igraph::get.edge.ids(g, c(node_id1, node_id2))
  
  # If the intensity of this edge was previously calculated, then return it
  if(edge_id != 0 & !is.null(igraph::edge_attr(g, "intensity", index = edge_id))){
    if(!is.na(igraph::edge_attr(g, "intensity", edge_id))[1]){
      return(igraph::edge_attr(g, 'intensity', index = edge_id))
    }
  }
  
  # Distance between the node and its neighbor
  edge_dist <- tryCatch(
    {
      abs(distances_mtx[node_id1, node_id2]) # Distance between the node and its neighbor 
    },
    # If the nodes are not part of the graph, give the proper information of the error.
    error=function(cond) {
      neighbors_list <- igraph::neighbors(g, node_id1)
      if(! igraph::V(g)[node_id2] %in% neighbors_list){
        message("Second vertice (node_id2) it's not a neighbor of first vertice (node_id1)")
      }else{
        message(cond)
      }
    }
  )    
  
  node1 <- c(igraph::vertex_attr(g, "xcoord", node_id1), igraph::vertex_attr(g, "ycoord", node_id1))
  node2 <- c(igraph::vertex_attr(g, "xcoord", node_id2), igraph::vertex_attr(g, "ycoord", node_id2))
  
  indicator <- 0
  # Counting events
  for(row in 1:nrow(event_coords)) {
    ep <- c(event_coords[row, 1], event_coords[row, 2])
    dist_obj <- list(p1 = node1, p2 = node2, ep = ep)
    class(dist_obj) <- 'netTools'
    d <- PointToSegment(dist_obj)
    
    # If the event is at a distance less or equal 'z' from the edge (segment) 
    # connecting both given points (the road), then is counted as an event of that road
    if(d <= z){
      indicator <- indicator + 1
    }
  }
  edge_intensity <- indicator / edge_dist
  
  edge_intensity
}


#' Calculate all the edge intensities of the graph. It's more accurate than using iteratively the 
#' function EdgeIntensity for all edges.
#' 
#' @name AllEdgeIntensities.intensitynet
#' 
#' @param obj intensitynet object
#' @param z Maximum distance between the event and the edge to consider the event part of the edge.
#' 
#' @return intensitynet class object where the graph contains all the edge intensities as an attribute
#' 
AllEdgeIntensities.intensitynet <- function(obj, z = 5){
  if(z <= 0){
    print("Warning: 'z' cannot be equal or less than 0, using default.")
    z <- 5
  }
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  event_coords <- obj$events
  edge_list <- igraph::ends(g, igraph::E(g), names=FALSE)
  
  if(length(event_coords) == 0){
    return(NA)
  }
  
  edge_events <- edge_list[,1]
  edge_events <- cbind(edge_events, edge_list[,2])
  edge_events <- cbind(edge_events, 0)
  edge_events <- cbind(edge_events, 0)
  colnames(edge_events) <- c('from', 'to', 'n_events', 'intensity')
  
  node_coords <- as.numeric(igraph::V(g))
  node_coords <- cbind(node_coords, igraph::vertex_attr(g, "xcoord"))
  node_coords <- cbind(node_coords, igraph::vertex_attr(g, "ycoord"))
  colnames(node_coords) <- c('node', 'xcoord', 'ycoord')
  
  #start_time <- Sys.time() # debug only
  pb = utils::txtProgressBar(min = 0, max = nrow(event_coords), initial = 0) 
  cat("Calculating edge intensities...\n")
  
  e_count <- 0
  for(row in 1:nrow(event_coords)){
    utils::setTxtProgressBar(pb, row)
    tmp_edge <- NULL
    shortest_d <- NULL
    for(edge_row in 1:nrow(edge_events)){
      # Check if the intensities are already calculated
      if(row == 1){
        if(!is.null(igraph::edge_attr(g, 'intensity', igraph::E(g)[edge_row]))){
          e_count <- e_count + 1
          next
        }
      }
      # Faster but only works if the node ID is the same as its index
      node1 <- node_coords[edge_events[edge_row, 'from'],][2:3]
      node2 <- node_coords[edge_events[edge_row, 'to'],][2:3]
      
      ep <- event_coords[row, ]
      dist_obj <- list(p1 = node1, p2 = node2, ep = ep)
      class(dist_obj) <- 'netTools'
      d <- PointToSegment(dist_obj)
      
      # If the event is at a distance less or equal 'z' from the edge (segment)
      # connecting both given points (the road), then is counted as an event of that road
      if(d <= z){
        if(d == 0){
          tmp_edge <- edge_row
          break
        }
        if (is.null(shortest_d) || d < shortest_d){
          tmp_edge <- edge_row
        }
      }
    }
    if (!is.null(tmp_edge)){
      edge_events[tmp_edge, 'n_events'] <- edge_events[tmp_edge, 'n_events'] + 1
    }
    # If the intensity of all edges is already calculated return the object
    if(row == 1){
      if(e_count == length(igraph::E(g))){
        print("Intensities were already calculated.")
        return(obj)
      } 
    }
  }
  close(pb)
  #cat(paste0("Time: ", Sys.time() - start_time, "\n")) # debug only
  
  for (edge_row in 1:nrow(edge_events)) {
    # Distance between the node and its neighbor
    edge_dist <- abs(distances_mtx[edge_events[edge_row, 'from'], edge_events[edge_row, 'to']])
    edge_events[edge_row, 'intensity'] <-  edge_events[edge_row, 'n_events'] / edge_dist
  }
  SetNetworkAttribute(obj = obj, 
                      where = 'edge', 
                      name = 'intensity', 
                      value = as.matrix(edge_events[, 'intensity']))
}


#' Calculates the intensity of the given path
#'
#' @name PathIntensity.intensitynet
#'
#' @param obj intensitynet object
#' @param path_nodes vector containing the node ID's of the path
#' 
#' @return intensity of the path
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' short_path <- ShortestPathIntensity(und_intnet_chicago, node_id1 = 'V1', node_id2 = 'V300')
#' PathIntensity(und_intnet_chicago, short_path$path)
#' 
#' @export
PathIntensity.intensitynet <- function(obj, path_nodes){
  edge_counts <- list()
  path_intensity <- 0
  
  prev <- NULL
  for(node_id in path_nodes){
    if(is.null(prev)){
      prev <- node_id
      next
    }
    
    path_intensity <- path_intensity + Reduce('+', EdgeIntensity(obj, prev, node_id))
    
    prev <- node_id
  }
  # Divide the intensity of the edges by their number (In a path -> N edges = N vertices - 1)
  path_intensity <- path_intensity / (length(path_nodes) - 1)
  
  path_intensity
}


#' Calculates the shortest path between two vertices and calculates its intensity
#'
#' @name ShortestPathIntensity.intensitynet
#'
#' @param obj intensitynet object
#' @param node_id1 starting node
#' @param node_id2 ending node
#' @param weighted TRUE or FALSE (default), tell if the distances must be taken into account 
#' 
#' @return intensity of the shortest path and the path vertices
#' 
#' @examples
#'
#' data("und_intnet_chicago")
#' ShortestPathIntensity(und_intnet_chicago, node_id1 = 'V1', node_id2 = 'V300')
#' 
#' @export
ShortestPathIntensity.intensitynet <- function(obj,  node_id1, node_id2, weighted = FALSE){
  g <- obj$graph
  
  if(weighted){
    path <- ShortestNodeDistance(obj, node_id1, node_id2)$path
  }else{
    path <- unlist(igraph::get.shortest.paths(g, node_id1, node_id2)$vpath)
  }
  
  return(list(intensity = PathIntensity(obj, path), path = path))
}


#' Gives general node correlation of the network (choosing from normal correlation, covariance, 
#' moran-i or geary)
#' 
#' @name NodeGeneralCorrelation.intensitynet
#'
#' @param obj intensitynet object
#' @param dep_type 'correlation', 'covariance', moran', 'geary'. The type of 
#' dependence statistic to be computed.
#' @param lag_max Maximum geodesic lag at which to compute dependence
#' @param intensity vector containing the intensity values that the heatmaps
#' 
#' @return A vector containing the dependence statistics (ascending from order 0). 
#' 
#' @examples 
#' 
#' data("und_intnet_chicago")
#' g <- und_intnet_chicago$graph
#' gen_corr <- NodeGeneralCorrelation(und_intnet_chicago, dep_type = 'correlation', lag_max = 2, 
#'                                    intensity = igraph::vertex_attr(g)$intensity)
#' 
#' @export
NodeGeneralCorrelation.intensitynet <- function(obj, dep_type, lag_max, intensity){
  g <- obj$graph
  g_sna <- intergraph::asNetwork(g)
  sna::nacf(g_sna, intensity, type = dep_type, mode = "graph", lag.max = lag_max)
}


#' Gives node local moran-i or geary-c correlations
#' 
#' @name NodeLocalCorrelation.intensitynet
#' 
#' @source *Luc Anselin. A Local Indicator of Multivariate SpatialAssociation: Extending Geary's c, Geographical Analysis 2018; doi: https://doi.org/10.1111/gean.12164
#'
#' @param obj intensitynet object
#' @param dep_type 'moran', 'getis' or 'geary'. Type of local correlation to be computed (Moran-i, Getis-Gstar, Geary-c*),
#' default = 'moran'.
#' @param intensity vector containing the intensity values that which are used to calculate the correlation.
#' 
#' @return intensitynet class object which contains an igraph network with the selected correlation 
#' added into the vertices attributes
#' 
#' @examples 
#' 
#' data("und_intnet_chicago")
#' g <- und_intnet_chicago$graph
#' data_moran <- NodeLocalCorrelation(und_intnet_chicago, 
#'                                    dep_type = 'moran', 
#'                                    intensity = igraph::vertex_attr(g)$intensity)
#' moran_i <- data_moran$correlation
#' intnet <- data_moran$intnet
#' 
#' @export
NodeLocalCorrelation.intensitynet <- function(obj, dep_type = 'moran', intensity){
  g <- obj$graph
  adj_mtx <- igraph::as_adj(graph = g)
  adj_listw <- spdep::mat2listw(adj_mtx)
  nb <- adj_listw$neighbours
  w_listw <- spdep::nb2listw(nb, style="W", zero.policy=TRUE) 
  
  if(dep_type == 'geary'){
    nb_b <- spdep::listw2mat(w_listw)
    b <- methods::as(nb_b, "CsparseMatrix")
    all(b == Matrix::t(b))
    
    # scale, with default settings, will calculate the mean and standard deviation of the entire vector,
    # then "scale" each element by those values by subtracting the mean and dividing by the sd
    val <- scale(intensity)[,1]
    
    # Calculates local geary-c.
    # Details in paper: A Local Indicator of Multivariate SpatialAssociation: Extending Geary’s c, from Luc Anselin
    n <- length(intensity)
    locgc <- numeric(n)
    for (i in c(1:n)) {
      locgc[i] <- sum(b[i,] * (val[i] - val)^2)
    }
    
    #--------------------------------------Comprovation:---------------------------------------
    # General Geary-c from local:
    # general <- sum(locgc)/(2*sum(nb_b))
    # 
    # g_sna <- intergraph::asNetwork(g)
    # general_ref <- sna::nacf(g_sna, intensity, type = 'geary', mode = "graph")[2]
    #------------------------------------------------------------------------------------------
    
    intnet <- SetNetworkAttribute(obj = obj, where = 'vertex', name = "geary", value = locgc)
    return(list(correlation = locgc, intnet = intnet))
    
  } else if (dep_type == 'moran'){
    locmoran <- spdep::localmoran(x = intensity, listw = w_listw, zero.policy=TRUE)
    intnet <- SetNetworkAttribute(obj = obj, where = 'vertex', name = 'moran', value = locmoran[, 'Ii'])
    return(list(correlation = locmoran, intnet = intnet))
    
  } else if (dep_type == 'getis'){
    b_listw <- spdep::nb2listw(nb, style="B", zero.policy=TRUE) 
    locgg <- spdep::localG(x = intensity, listw = b_listw)
    intnet <- SetNetworkAttribute(obj = obj, where = 'vertex', name = "getis", value = locgg)
    return(list(correlation = locgg, intnet = intnet))
  }else{
    stop("Error: wrong 'dep_type' argument. Allowed arguments are; 'moran', 'getis', 'geary'.")
  }
}


#' Plot the network and if specified, the correlation heatmap. Which could be:
#'
#' @name PlotHeatmap.intensitynet
#'
#' @param obj intensitynet object
#' @param heattype 'moran': Local Moran-i correlation (with 999 permutations), 'geary': Local Geary-c* 
#' correlation. The correlations will use the indicated intensity type.
#' The function also allow to only plot the intensity heatmap 'v_intensity' for vertices or 'e_intensity' for edges.
#' @param intensity_type name of the intenisty used to plot the heatmap. For undirected networks: 'intensity'. 
#' For directed networks: 'intensity_in' or 'intensity_out'. For mixed networks: 'intensity_in', 'intensity_out', 
#' 'intensity_und' or 'intensity_all'. If the intensity parameter is NULL, the function will use, if exist, 
#' the intensity (undirected) or intensity_in (directed) values from the network nodes.
#' @param net_vertices chosen vertices to plot the heatmap (or it related edges in case to plot the edge heatmap)
#' @param ... extra arguments for the class ggplot
#' 
#' @return The plot of the heatmap with class c("gg", "ggplot")
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' PlotHeatmap(und_intnet_chicago, heattype='moran')
#' 
#' @export
PlotHeatmap.intensitynet <- function(obj, heattype = 'none', intensity_type = 'none', net_vertices = NULL, ...){
  g <- obj$graph
  adj_mtx <- igraph::as_adj(graph = g)
  adj_listw <- spdep::mat2listw(adj_mtx)
  nb <- adj_listw$neighbours
  w_listw <- spdep::nb2listw(nb, style = "W",  zero.policy=TRUE)
  
  if(heattype != 'none' && heattype != 'moran' && heattype != 'geary' && 
     heattype != 'v_intensity' && heattype != 'e_intensity' && heattype != 'intensity'){
    warning('Parameter "heattype" should be "moran", "geary", "intensity", 
            "v_intensity", "e_intensity" or "none". Using default ("none").')
  }
  
  if (is.null(net_vertices)){
    net_vertices <- igraph::V(g)
  }
  
  # If the intensity is not provided, try to take it from the given network
  if( intensity_type == 'none' || is.null(igraph::vertex_attr(graph = g, name = intensity_type)) ){
    if( intensity_type != 'none' && is.null(igraph::vertex_attr(graph = g, name = intensity_type))){
      warning(paste0("The given intensity type '", intensity_type, "', doestn't match any of the network intensities."))
    }
    if(!is.null(igraph::vertex_attr(graph = g, name = 'intensity'))){
      intensity <- igraph::vertex_attr(graph = g, name = 'intensity')
      intensity_type <- 'intensity'
    }else if(!is.null(igraph::vertex_attr(graph = g, name = 'intensity_in'))){
      intensity <- igraph::vertex_attr(graph = g, name = 'intensity_in')
      intensity_type <- 'intensity_in'
    }else{
      intensity <- NA
    }
  }
  
  node_coords <- data.frame(xcoord = igraph::vertex_attr(graph = g, name = 'xcoord'), 
                            ycoord = igraph::vertex_attr(graph = g, name = 'ycoord'))
  rownames(node_coords) <- igraph::vertex_attr(graph = g, name = 'name')
  
  if(heattype == 'moran'){ # Local Moran-i
    locmoran <- spdep::localmoran_perm(x = intensity, 
                                       listw = w_listw, 
                                       zero.policy = TRUE, 
                                       na.action = stats::na.omit, 
                                       nsim = 999) # 999 permutations
    
    # Calculate deviations
    node_int_deviation <- intensity - mean(intensity)  
    locmoran_deviation <- locmoran[, 'Ii'] - mean(locmoran[, 'Ii'])
    
    # create a new variable identifying the moran plot quadrant for each observation, dismissing the non-significant ones
    sig_dstr <- NA
    significance <- 0.05
    
    # non-significant observations
    sig_dstr[(locmoran[, 5] > significance)] <- 2 # "insignificant"  
    # low-low quadrant
    sig_dstr[(node_int_deviation < 0 & locmoran_deviation < 0) & (locmoran[, 5] <= significance)] <- 3 # "low-low"
    # low-high quadrant
    sig_dstr[(node_int_deviation < 0 & locmoran_deviation > 0) & (locmoran[, 5] <= significance)] <- 4 # "low-high"
    # high-low quadrant
    sig_dstr[(node_int_deviation > 0 & locmoran_deviation < 0) & (locmoran[, 5] <= significance)] <- 5 # "high-low"
    # high-high quadrant
    sig_dstr[(node_int_deviation > 0 & locmoran_deviation > 0) & (locmoran[, 5] <= significance)] <- 6 # "high-high"
    
    sig_dstr[setdiff(as.numeric(igraph::V(g)), net_vertices)] <- 1 # Not contemplated vertices 
    
    data_df <- data.frame(xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord, 
                          value = sig_dstr)
    
  }else if(heattype == 'geary'){  # Local Geary-c
    nb_b <- spdep::listw2mat(w_listw)
    
    b <- methods::as(nb_b, "CsparseMatrix")
    all(b == Matrix::t(b))
    
    # scale, with default settings, will calculate the mean and standard deviation of the entire vector,
    # then "scale" each element by those values by subtracting the mean and dividing by the sd
    val <- scale(intensity)[,1]
    
    # Calculates local geary-c.
    # Details in paper: A Local Indicator of Multivariate SpatialAssociation: Extending Geary’s c, from Luc Anselin
    n <- length(intensity)
    locgc <- numeric(n)
    for (i in c(1:n)) {
      locgc[i] <- sum(b[i,] * (val[i] - val)^2)
    }
    
    
    #--------------------------------------Comprovation:---------------------------------------
    # General Geary-c from local:
    # general <- sum(locgc)/(2*sum(nb_b))
    # 
    # g_sna <- intergraph::asNetwork(g)
    # general_ref <- sna::nacf(g_sna, intensity, type = 'geary', mode = "graph")[2]
    #------------------------------------------------------------------------------------------
    
    
    # create a new variable identifying the moran plot quadrant for each observation, dismissing the non-significant ones
    sig_dstr <- NA
    
    sig_dstr[locgc < 1] <- 2 # positive spatial autocorrelation
    sig_dstr[locgc == 1 || is.na(locgc)] <- 3 # no spatial autocorrelation
    sig_dstr[locgc > 1] <- 4 # negative spatial autocorrelation
    
    sig_dstr[setdiff(as.numeric(igraph::V(g)), net_vertices)] <- 1 # Not contemplated vertices 
    
    data_df <- data.frame(xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord, 
                          value = sig_dstr)
    
  }else if(heattype == 'getis'){ # Local Getis-G*
    print("Needs implementation")
    # TODO: Implement Getis G.
    
    # locgg <- spdep::localG(x = intensity, listw = b_listw)
    # 
    # # create a new variable identifying the moran plot quadrant for each observation, dismissing the non-significant ones
    # sig_dstr <- NA
    # 
    # data_df <- data.frame(xcoord = node_coords$xcoord, 
    #                       ycoord = node_coords$ycoord, 
    #                       value = sig_dstr)
    
    
    
  }else if(heattype == 'v_intensity'){
    norm_int <- (intensity - min(intensity)) / (max(intensity) - min(intensity))
    data_df <- data.frame(xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord, 
                          value = norm_int)
    
  }else if(heattype == 'e_intensity'){
    data_df <- data.frame(xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord)
    
  }else if(heattype == 'intensity'){
    norm_int <- (intensity - min(intensity)) / (max(intensity) - min(intensity))
    data_df <- data.frame(xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord, 
                          value = norm_int)
    
  }else{
    data_df <- data.frame(xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord, 
                          value = NA)
  }
  geoplot_obj <- list(graph = g, data_df = data_df, net_vertices = net_vertices, mode = heattype)
  class(geoplot_obj) <- "netTools"
  
  GeoreferencedGgplot2(geoplot_obj, ...)
}


#' Plot the net and the events in the neighborhood area of the given node
#' 
#' @name PlotNeighborhood.intensitynet
#' 
#' @param obj Intensitynet object
#' @param node_id Id of the node which the plot will be focused
#' @param ... Extra arguments for plotting
#' 
#' @return No return value, just plots the neighborhood and the events.
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' PlotNeighborhood(und_intnet_chicago, node_id = 'V300')
#' 
#' @export
PlotNeighborhood.intensitynet <- function(obj, node_id, ...){
  g <- obj$graph
  events <- obj$events
  w_margin <- 50
  
  nei <- igraph::neighbors(g, node_id)
  
  v_coords <- cbind(igraph::vertex_attr(g, 'xcoord', nei), igraph::vertex_attr(g, 'ycoord', nei))
  v_coords <- rbind(v_coords, c(igraph::vertex_attr(g, 'xcoord', node_id), igraph::vertex_attr(g, 'ycoord', node_id)))
  colnames(v_coords) <- c('xcoord', 'ycoord')
  rownames(v_coords) <- c(names(nei), node_id)
  
  window_coords <- list(min_x = min(v_coords[, 'xcoord']), min_y = min(v_coords[, 'ycoord']),
                        max_x = max(v_coords[, 'xcoord']), max_y = max(v_coords[, 'ycoord']))
  
  e_coords <- NULL
  for(row in 1:nrow(events)){
    if(events[row, 1] >= window_coords$min_x - w_margin & 
       events[row, 1] <= window_coords$max_x + w_margin & 
       events[row, 2] >= window_coords$min_y - w_margin & 
       events[row, 2] <= window_coords$max_y + w_margin){
      e_coords <- rbind(e_coords, events[row,])
    }
  }
  if(!is.null(e_coords)) colnames(e_coords) <- c('xcoord', 'ycoord')
  
  # Plot vertices
  graphics::plot(v_coords, xlim = c(window_coords$min_x - w_margin , window_coords$max_x + w_margin), 
                 ylim = c(window_coords$min_y - w_margin , window_coords$max_y + w_margin), ...)
  graphics::text(x = v_coords[, 'xcoord'], y = v_coords[, 'ycoord'], c(names(nei), node_id), cex = 1, col = 'blue')
  
  # Draw edges
  for(row in 1:nrow(v_coords)-1){
    graphics::lines(x = c(v_coords[nrow(v_coords),]['xcoord'], v_coords[row,]['xcoord']),
                    y = c(v_coords[nrow(v_coords),]['ycoord'], v_coords[row,]['ycoord']),
                    type = "l", lty = 1)
  }
  
  # Plot events
  graphics::points(x = e_coords[,'xcoord'], y = e_coords[,'ycoord'], col = 'red')
}


#' Set attributes to the network edges or nodes
#'
#' @name SetNetworkAttribute.intensitynet
#'
#' @param obj intensitynet object
#' @param where 'vertex' or 'edge', where to set the attribute
#' @param name name of the attribute
#' @param value vector containing the data for the attribute
#' 
#' @return intensitynet object containing the network with the added attributes
#' 
SetNetworkAttribute.intensitynet <- function(obj, where, name, value){
  g <- obj$graph
  
  if(where == 'edge'){
    #g <- g %>% igraph::set_edge_attr(name = name, value = value)
    g <- igraph::set_edge_attr(g, name = name, value = value)
  }else{
    #g <- g %>% igraph::set_vertex_attr(name = name, value = value)
    g <- igraph::set_vertex_attr(g, name = name, value = value)
  } 
  
  intnet <- list(graph = g, 
                 events = obj$events, 
                 graph_type = obj$graph_type, 
                 distances_mtx = obj$distances_mtx)
  attr(intnet, 'class') <- class(obj)
  intnet
}


#' Get the intensitynet object delimited by the given window
#' 
#' @name ApplyWindow.intensitynet
#' 
#' @param obj intensitynet object
#' @param x_coords vector containing the x coordinate limits of the window
#' @param y_coords vector containing the y coordinate limits of the window
#' 
#' @return intensitynet object delimited by the window (sub-part of the original)
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' sub_intnet_chicago <- ApplyWindow(und_intnet_chicago, 
#'                                   x_coords = c(300, 900), 
#'                                   y_coords = c(500, 1000))
#' 
#' @export
ApplyWindow.intensitynet <- function(obj, x_coords, y_coords){
  g <- obj$graph
  events <- obj$events
  nodes <- igraph::V(g)
  
  window_coords <- list(min_x = min(x_coords), min_y = min(y_coords),
                        max_x = max(x_coords), max_y = max(y_coords))
  
  sub_v_coords <- c()
  for(row in 1:length(nodes)){
    v <- igraph::V(g)[row]
    vertex_info <- igraph::vertex_attr(graph = g, index = v)
    if(vertex_info$xcoord >= window_coords$min_x & 
       vertex_info$xcoord <= window_coords$max_x & 
       vertex_info$ycoord >= window_coords$min_y & 
       vertex_info$ycoord <= window_coords$max_y){
      sub_v_coords <- c(sub_v_coords, v)
    }
  }
  
  sub_e_coords <- NULL
  for(row in 1:nrow(events)){
    if(events[row, 1] >= window_coords$min_x & 
       events[row, 1] <= window_coords$max_x & 
       events[row, 2] >= window_coords$min_y & 
       events[row, 2] <= window_coords$max_y){
      sub_e_coords <- rbind(sub_e_coords, events[row,])
    }
  }
  if(!is.null(sub_e_coords)) colnames(sub_e_coords) <- c('xcoord', 'ycoord')
  
  # Get sub-graph
  sub_g <- igraph::induced_subgraph(graph = g, vids = sub_v_coords)
  
  # Keep only the biggest component
  # sub_g_components <- igraph::components(sub_g)$membership
  # greatest_sub <- names(which(table(sub_g_components) == max(table(sub_g_components))))
  # sub_g <- induced_subgraph(graph = g, vids = names(which(sub_g_components == greatest_sub))) 
  
  node_coords_obj <- list(node_coords = cbind(igraph::vertex_attr(sub_g, 'xcoord'), 
                                              igraph::vertex_attr(sub_g, 'ycoord')))
  class(node_coords_obj) <- "netTools"
  sub_dist_mtx <- CalculateDistancesMtx(node_coords_obj)
  
  intnet <- list(graph = sub_g, 
                 events = sub_e_coords, 
                 graph_type = obj$graph_type, 
                 distances_mtx = sub_dist_mtx)
  attr(intnet, 'class') <- class(obj)
  
  intnet
}


#' Calculates the shortest distance path between two nodes 
#'
#' @name ShortestNodeDistance.intensitynet
#'
#' @param obj intensitynet object
#' @param node_id1 id of the starting node
#' @param node_id2 id of the end node
#' 
#' @return distance of the path and the nodes of the path
#'
ShortestNodeDistance.intensitynet <- function(obj, node_id1, node_id2){
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  
  weighted_path <- unlist(igraph::get.shortest.paths(g, node_id1, node_id2)$vpath)
  if(!is.null(distances_mtx)){
    weight_sum <- sum(igraph::E(g, path = unlist(weighted_path))$weight)
  }
  else{
    weight_sum <- length(weighted_path)
  }
  list(weight = weight_sum, path = weighted_path)  
}