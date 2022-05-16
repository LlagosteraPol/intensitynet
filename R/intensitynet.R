#' Constructor of the class intensitynet. In order to create an intensitynet object, it is needed; an adjacency matrix, the
#' coordinates of the nodes and the coordinates of the events.
#'
#' @name intensitynet
#'
#' @param adjacency_mtx Network adjacency matrix
#' @param node_coords Nodes latitude and longitude matrix (coordinates)
#' @param event_data DataFrame with event latitude and longitude coordinates (mandatory columns) and optional attributes related to the events
#' @param graph_type Network type: 'undirected' (default), 'directed' or 'mixed' 
#' @param event_correction Value that determines how far can be an event to be considered part of a segment (default 5). 
#' This value highly depends on the given coordinate system
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
#' # Create the intensitynet object, in this case will be undirected 
#' intnet_chicago <- intensitynet(chicago_adj_mtx, 
#'                                node_coords = chicago_node_coords, 
#'                                event_data = chicago_df)
#' 
#' @export
intensitynet <- function(adjacency_mtx, node_coords, event_data, graph_type = 'undirected', event_correction = 5){
  
  if(event_correction < 0){
    message("Warning: event correction value cannot be less than 0, using default.")
    event_correction <- 5
  }
  
  if (is.data.frame(adjacency_mtx)) {
    adjacency_mtx <- as.matrix(adjacency_mtx)
  }
  
  if (is.data.frame(node_coords)) {
    node_coords <- as.matrix(node_coords)
  }
  colnames(node_coords) <- c("xcoord", "ycoord")
  
  
  if (is.matrix(event_data)) {
    event_data <- as.data.frame(event_data)
  }
  names(event_data)[1:2] <- c("xcoord", "ycoord")
  
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
  
  intnet <- list(graph = g, events = event_data, graph_type = graph_type, 
                 distances_mtx = dist_mtx, event_correction = event_correction)
  attr(intnet, "class") <- "intensitynet"
  # Select the proper class
  switch(graph_type, 
         'undirected' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetUnd")},
         'directed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetDir")},
         'mixed' = {attr(intnet, 'class') <- c(class(intnet), "intensitynetMix")})
  
  intnet # return
}


# -------- Network functions ----------

#' It allows to compute different dependence statistics on the network
#' for the given vector and for neighborhoods of distinct order. Such statistics are; correlation,
#' covariance, Moran’s I and Geary’s C. 
#' 
#' @name NodeGeneralCorrelation
#'
#' @param obj intensitynet object
#' @param dep_type 'correlation', 'covariance', moran', 'geary'. The type of 
#' dependence statistic to be computed.
#' @param lag_max Maximum geodesic lag at which to compute dependence
#' @param intensity Vector containing the values to calculate the specified dependency in the network. Usually the node mean intensities.
#' @param partial_neighborhood use partial neighborhood (TRUE) or cumulative (FALSE). TRUE by default
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
NodeGeneralCorrelation <- function(obj, dep_type, lag_max, intensity, partial_neighborhood = TRUE){
  UseMethod("NodeGeneralCorrelation")
}


#' Gives the node local Moran-I, Getis-Gstar or Geary-c correlations
#' 
#' @name NodeLocalCorrelation
#' 
#' @source *"A Local Indicator of Multivariate SpatialAssociation: Extending Geary's c, Geographical Analysis" Luc Anselin (2018) <doi:10.1111/gean.12164>
#'
#' @param obj intensitynet object
#' @param dep_type 'moran', 'getis' or 'geary'. Type of local correlation to be computed (Moran-i, Getis-Gstar, Geary-c),
#' default = 'moran'.
#' @param intensity vector containing the values to calculate the specified correlation for each node in the network.
#' 
#' @return a vector containing two values. The first value is a vector with the specified local correlations for each node. 
#' The second values is the  given intensitynet class object but with the correlations added to the node attributes of its network. 
#' 
#' @examples 
#' \dontrun{
#' data("und_intnet_chicago")
#' g <- und_intnet_chicago$graph
#' data_moran <- NodeLocalCorrelation(und_intnet_chicago, 
#'                                    dep_type = 'moran', 
#'                                    intensity = igraph::vertex_attr(g)$intensity)
#' moran_i <- data_moran$correlation
#' intnet <- data_moran$intnet
#' }
#' @export
NodeLocalCorrelation <- function(obj, dep_type = 'moran', intensity){
  UseMethod("NodeLocalCorrelation")
}


#' Plot the network correlations or intensities.
#'
#' @name PlotHeatmap
#'
#' @param obj intensitynet object
#' @param heat_type a string with the desired heatmap to be plotted, the options are; 
#' 'moran': Local Moran-i correlation (with 999 permutations), 
#' 'geary': Local Geary-c correlation. The correlations will use the indicated intensity type,
#' 'v_intensity': vertice mean intensity,
#' 'e_intensity': edge intensity,
#' mark name: name of the mark (string) to plot its edge proportion,
#' 'none': plain map.
#' @param intensity_type name of the vertex intensity used to plot the heatmap for moran, geary and v_intensity options (of the heat_type argument).
#' The options are; 
#' For undirected networks: 'intensity'. 
#' For directed networks: 'intensity_in' or 'intensity_out'. For mixed networks: 'intensity_in', 'intensity_out', 
#' 'intensity_und' or 'intensity_all'. If the intensity parameter is 'none', the function will use, if exist, 
#' the intensity (undirected) or intensity_in (directed) values from the network nodes. If the heat_type is 'e_intensity', this
#' parameter will be skiped and plot the edge intensities instead.
#' @param net_vertices chosen vertices to plot the heatmap (or its related edges in case to plot the edge heatmap)
#' @param net_edges chosen edges to plot the heatmap, can be either the edge id's or its node endpoints (e.j. c(1,2, 2,3, 7,8))
#' @param show_events option to show the events as orange squares, FALSE by default
#' @param alpha optional argument to set the transparency of the events (show_events = TRUE). The range is from 0.1 (transparent) to 1 (opaque). Default: alpha = 1
#' @param ... extra arguments for the class ggplot
#' 
#' @return The plot of the heatmap with class c("gg", "ggplot")
#' 
#' @examples
#' 
#' \dontrun{
#' data("und_intnet_chicago")
#' PlotHeatmap(und_intnet_chicago, heat_type='moran')
#' }
#' 
#' @export
PlotHeatmap <- function(obj, heat_type = 'none', intensity_type = 'none', net_vertices = NULL, net_edges = NULL, show_events = FALSE, alpha = 1, ...){
  UseMethod("PlotHeatmap")
}


#' Plot the net and the events in the neighborhood area of the given node
#' 
#' @name PlotNeighborhood
#' 
#' @param obj intensitynet object
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

#' Calculates the total weight of the given path
#'
#' @name PathTotalWeight
#'
#' @param obj intensitynet object
#' @param path_nodes vector containing the node ID's of the path
#' @param weight an string specfiying the type of weight to be computed. If no weight type is provided,
#' the function will calculate the toatl amount of edges. Default NA.
#' 
#' @return total weight of the path
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' PathTotalWeight(und_intnet_chicago, c('V115', 'V123', 'V125', 'V134'), weight = 'intensity')
#' 
#' @export
PathTotalWeight <- function(obj, path_nodes, weight = NA){
  UseMethod("PathTotalWeight")
}


#' Calculates the shortest path between two vertices (based on the minimum amount of edges) and 
#' calculates its toatl weight
#'
#' @name ShortestPath
#'
#' @param obj intensitynet object
#' @param node_id1 starting node
#' @param node_id2 ending node
#' @param weight an string, calculate the shortest path based on this type of weight. If no weight type is provided,
#' the function will calculate the shortest path based on the minimum amount of edges. Default NA.
#' @param mode Character 'in', 'out', 'all' (default). Gives whether the shortest paths to or from the given vertices 
#' should be calculated for directed graphs. If out then the shortest paths from the vertex, if in 
#' then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not 
#' directed paths are searched. This argument is ignored for undirected graphs.
#' 
#' @return total weight of the shortest path and the path vertices with class igraph.vs
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' ShortestPath(und_intnet_chicago, node_id1 = 'V1', node_id2 = 'V300', weight = 'intensity')
#' 
#' @export
ShortestPath <- function(obj,  node_id1, node_id2, weight = NA, mode = 'all'){
  UseMethod("ShortestPath")
}

#' Calculates edgewise and mean nodewise intensities for the given intensitynet object and, for each edge, the proportions of
#' all event covariates.
#' 
#' @name RelateEventsToNetwork
#' 
#' @param obj intensitynet object
#' 
#' @return proper intensitynet object (Undirected, Directed, or Mixed) with a graph containing the nodewise intensity in the node 
#' attributes and the edgewise intensities and event covariate proportions as edge attributes.
#' 
#' @examples 
#'
#' data("und_intnet_chicago")
#' intnet_chicago <- RelateEventsToNetwork(und_intnet_chicago)
#' 
#' @export
RelateEventsToNetwork <- function(obj){
  UseMethod("RelateEventsToNetwork")
}


MeanNodeIntensity <- function(obj, node_id){
  UseMethod("MeanNodeIntensity")
}


EdgeIntensity <- function(obj, node_id1, node_id2){
  UseMethod("EdgeIntensity")
}


EdgeIntensitiesAndProportions <- function(obj){
  UseMethod("EdgeIntensitiesAndProportions")
}


SetNetworkAttribute <- function(obj, where, name, value){
  UseMethod("SetNetworkAttribute")
}


#' If not calculated, calculates the intensity of the edge with nodes; node_id1, node_id2. 
#' If the edge already contains an intensity, give it directly.
#'
#' @name EdgeIntensity.intensitynet
#' 
#' @param obj intensitynet object
#' @param node_id1 First node ID of the edge
#' @param node_id2 Second node ID of the edge
#' 
#' @return Intensity of the edge
#' 
EdgeIntensity.intensitynet <- function(obj,  node_id1, node_id2){
  if(node_id1 == node_id2){
    stop("The two vertices cannot be the same.")
  }
  
  if(obj$event_correction < 0){
    message("Warning: event correction value cannot be less than 0, using default.")
    z <- 5
  }
  else{
    z <- obj$event_correction
  }
  
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  event_data <- obj$events
  
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
  for(row in 1:nrow(event_data)) {
    ep <- c(event_data[row, 1], event_data[row, 2])
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


#' Calculate all the edge intensities of the graph. It's more fast than using iteratively the 
#' function EdgeIntensity for all edges.
#' 
#' @name EdgeIntensitiesAndProportions.intensitynet
#' 
#' @param obj intensitynet object
#' 
#' @return intensitynet class object where the graph contains all the edge intensities as an attribute
#' 
EdgeIntensitiesAndProportions.intensitynet <- function(obj){
  if(obj$event_correction < 0){
    message("Warning: event correction value cannot be less than 0, using default.")
    z <- 5
  }
  else{
    z <- obj$event_correction
  }
  
  g <- obj$graph
  distances_mtx <- obj$distances_mtx
  event_data <- obj$events
  edge_list <- igraph::ends(g, igraph::E(g), names=FALSE)
  
  if(length(event_data) == 0){
    return(NA)
  }
  
  from_coords <- igraph::vertex_attr(graph = g, name = "xcoord", index = edge_list[,1])
  from_coords <- cbind(from_coords, igraph::vertex_attr(graph = g, name = "ycoord", index = edge_list[,1]))
  
  to_coords <- igraph::vertex_attr(graph = g, name = "xcoord", index = edge_list[,2])
  to_coords <- cbind(to_coords, igraph::vertex_attr(graph = g, name = "ycoord", index = edge_list[,2]))
  
  # Prepare structure to set information in the edges
  edge_events <- data.frame(from = edge_list[,1], 
                            to = edge_list[,2],
                            n_events = 0,
                            intensity = 0)
  
  # Set up the names for the covariates of the events
  if(ncol(event_data) > 2){
    for(i in 3:ncol(event_data)){
      if(is.numeric(event_data[,i])){
        edge_events[colnames(event_data[i])] <- 0 # Set event column name
      }else{
        edge_events[as.character(unique(event_data[[i]]))] <- 0 # Set event unique variables as names
      }
    }
  }
  
  message(paste0("\nCalculating edge intensities with event error distance of ", obj$event_correction ,"..."))
  pb = utils::txtProgressBar(min = 0, max = nrow(event_data), initial = 0, style=3) 
  
  e_count <- 0
  for(row in 1:nrow(event_data)){
    utils::setTxtProgressBar(pb, row)
    tmp_edge <- NULL
    shortest_d <- NULL
    
    ep <- event_data[row, ]
    dist_obj <- list(p1 = from_coords, p2 = to_coords, ep = ep[,1:2])
    class(dist_obj) <- 'netTools'
    event_seg_dist <- PointToSegment(dist_obj)
    
    closest_e <- which.min(event_seg_dist)
    
    # Check if the intensities are already calculated
    if(!is.null(igraph::edge_attr(g, 'intensity', igraph::E(g)[closest_e]))){
      e_count <- e_count + 1
      next
    }
    
    # Check if the closest edge is in the required boundary and, if so,
    # set up the information (except intensity) to the edge_event DataFrame
    if ( event_seg_dist[closest_e] <= z ){
      edge_events[closest_e, 'n_events'] <- edge_events[closest_e, 'n_events'] + 1
      
      if(ncol(ep) > 2){
        for(i_col in 3:ncol(ep)){
          if( is.numeric(ep[,i_col]) ){
            tmp_str <- colnames(ep[i_col])
            edge_events[closest_e, tmp_str] <- edge_events[closest_e, tmp_str] + ep[,i_col]
          }else{
            tmp_str <-  as.character(ep[,i_col])
            edge_events[closest_e, tmp_str] <- edge_events[closest_e, tmp_str] + 1
          }
        }
      }
    }
  }
  close(pb)
  # If the intensity of all edges is already calculated return the object
  if(e_count == nrow(event_data)){
    return(obj)
  } 
  
  #Calculate intensity and proportions
  for (edge_row in 1:nrow(edge_events)) {
    if(edge_events[edge_row, 'n_events'] > 0 ){
      # Distance between the node and its neighbor
      edge_dist <- abs(distances_mtx[edge_events[edge_row, 'from'], edge_events[edge_row, 'to']])
      edge_events[edge_row, 'intensity'] <-  edge_events[edge_row, 'n_events'] / edge_dist
      
      if(ncol(edge_events) > 4){
        for(i_col in 5:ncol(edge_events)){
          edge_events[edge_row, i_col] <- edge_events[edge_row, i_col] / edge_events[edge_row, 'n_events'] 
        }
      }
    }
  }
  
  # Save information from 'edge_events' to the edge attributes of the network
  for(i_col in 3:ncol(edge_events)){
    obj <- SetNetworkAttribute(obj = obj, 
                               where = 'edge', 
                               name = colnames(edge_events[i_col]), 
                               value = as.matrix(edge_events[, i_col]))
  }
  return(obj)
}


#' Calculates the total weight of the given path
#'
#' @name PathTotalWeight.intensitynet
#'
#' @param obj intensitynet object
#' @param path_nodes vector containing the node ID's of the path
#' @param weight an string specfiying the type of weight to be computed. If no weight type is provided,
#' the function will calculate the toatl amount of edges. Default NA.
#' 
#' @return total weight of the path
#' 
#' @examples
#' 
#' data("und_intnet_chicago")
#' PathTotalWeight(und_intnet_chicago, c('V115', 'V123', 'V125', 'V134'), weight = 'intensity')
#' 
#' @export
PathTotalWeight.intensitynet <- function(obj, path_nodes, weight = NA){
  g <- obj$graph
  
  if(!is.na(weight) && !(weight %in% igraph::edge_attr_names(g))){
    warning("The given weight doens't exist in the edge attributes, using default instead (NA)")
    weight = NA
  }
  
  path_edges <- igraph::E(g, path = c(1,2,5,7))
  
  if (is.na(weight)){
    total_weight <- length(path_edges)
  }else{
    total_weight <- sum(igraph::edge_attr(g, name = weight, index = path_edges))
  }
  
  total_weight
}


#' Calculates the shortest path between two vertices (based on the minimum amount of edges) and 
#' calculates its total weight
#'
#' @name ShortestPath.intensitynet
#'
#' @param obj intensitynet object
#' @param node_id1 starting node
#' @param node_id2 ending node
#' @param weight an string, calculate the shortest path based on this type of weight. If no weight type is provided,
#' the function will calculate the shortest path based on the minimum amount of edges. Default NA.
#' @param mode Character 'in', 'out', 'all' (default). Gives whether the shortest paths to or from the given vertices 
#' should be calculated for directed graphs. If out then the shortest paths from the vertex, if in 
#' then to it will be considered. If all, the default, then the corresponding undirected graph will be used, ie. not 
#' directed paths are searched. This argument is ignored for undirected graphs.
#' 
#' @return total weight of the shortest path and the path vertices with class igraph.vs
#' 
#' @examples
#'
#' data("und_intnet_chicago")
#' ShortestPath(und_intnet_chicago, node_id1 = 'V1', node_id2 = 'V300', weight = 'intensity')
#' 
#' @export
ShortestPath.intensitynet <- function(obj,  node_id1, node_id2, weight = NA, mode = 'all'){
  g <- obj$graph
  
  if(!is.na(weight) && !(weight %in% igraph::edge_attr_names(g))){
    warning("The given weight doens't exist in the edge attributes, using default instead (NA)")
    weight = NA
  }
  
  if (is.na(weight)){
    path_data <- igraph::shortest_paths(graph = g, from = node_id1, to = node_id2, mode = mode, weights = NA, output = 'both')
    total_weight <- length(path_data$vpath[[1]])
  }else{
    
    # Remove any decimal by multiplitying with the maximum amount of possible decimals, then add 1 to remove 0 weight edges,
    # this is done to make sure that the shortest path igraph function take into account all decimals
    weight_vector <- igraph::edge_attr(g)[[weight]] * 100000000 + rep(1, igraph::gsize(g))
    
    path_data <- igraph::shortest_paths(graph = g, from = node_id1, to = node_id2, mode = mode, weights = weight_vector, output = 'both')
    total_weight <- sum(igraph::edge_attr(graph = g, weight, index = path_data$vpath[[1]]))
  }
  return(list(total_weight = total_weight, path = path_data$vpath[[1]]))
}


#' It allows to compute different dependence statistics on the network
#' for the given vector and for neighborhoods of distinct order. Such statistics are; correlation,
#' covariance, Moran’s I and Geary’s C. 
#' 
#' @name NodeGeneralCorrelation.intensitynet
#'
#' @param obj intensitynet object
#' @param dep_type 'correlation', 'covariance', moran', 'geary'. The type of 
#' dependence statistic to be computed.
#' @param lag_max Maximum geodesic lag at which to compute dependence
#' @param intensity Vector containing the values to calculate the specified dependency in the network. Usually the node mean intensities.
#' @param partial_neighborhood use partial neighborhood (TRUE) or cumulative (FALSE). TRUE by default
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
NodeGeneralCorrelation.intensitynet <- function(obj, dep_type, lag_max, intensity, partial_neighborhood = TRUE){
  g <- obj$graph
  g_sna <- intergraph::asNetwork(g)
  
  if(obj$graph_type == 'undirected') m <- 'graph'
  else m <- 'digraph'
  
  sna::nacf(g_sna, intensity, type = dep_type, mode = m, lag.max = lag_max, partial.neighborhood = partial_neighborhood)
}


#' Gives the node local Moran-I, Getis-Gstar or Geary-c correlations
#' 
#' @name NodeLocalCorrelation.intensitynet
#' 
#' @source *Luc Anselin. A Local Indicator of Multivariate SpatialAssociation: Extending Geary's c, Geographical Analysis 2018; doi: https://doi.org/10.1111/gean.12164
#'
#' @param obj intensitynet object
#' @param dep_type 'moran', 'getis' or 'geary'. Type of local correlation to be computed (Moran-i, Getis-Gstar, Geary-c),
#' default = 'moran'.
#' @param intensity vector containing the values to calculate the specified correlation for each node in the network.
#' 
#' @return a vector containing two values. The first value is a vector with the specified local correlations for each node. 
#' The second values is the  given intensitynet class object but with the correlations added to the node attributes of its network. 
#' 
#' @examples 
#' \dontrun{
#' data("und_intnet_chicago")
#' g <- und_intnet_chicago$graph
#' data_moran <- NodeLocalCorrelation(und_intnet_chicago, 
#'                                    dep_type = 'moran', 
#'                                    intensity = igraph::vertex_attr(g)$intensity)
#' moran_i <- data_moran$correlation
#' intnet <- data_moran$intnet
#' }
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


#' Plot the network correlations or intensities.
#'
#' @name PlotHeatmap.intensitynet
#'
#' @param obj intensitynet object
#' @param heat_type a string with the desired heatmap to be plotted, the options are; 
#' 'moran': Local Moran-i correlation (with 999 permutations), 
#' 'geary': Local Geary-c correlation. The correlations will use the indicated intensity type,
#' 'v_intensity': vertice mean intensity,
#' 'e_intensity': edge intensity,
#' mark name: name of the mark (string) to plot its edge proportion,
#' 'none': plain map.
#' @param intensity_type name of the vertex intensity used to plot the heatmap for moran, geary and v_intensity options (of the heat_type argument).
#' The options are; 
#' For undirected networks: 'intensity'. 
#' For directed networks: 'intensity_in' or 'intensity_out'. For mixed networks: 'intensity_in', 'intensity_out', 
#' 'intensity_und' or 'intensity_all'. If the intensity parameter is 'none', the function will use, if exist, 
#' the intensity (undirected) or intensity_in (directed) values from the network nodes. If the heat_type is 'e_intensity', this
#' parameter will be skiped and plot the edge intensities instead.
#' @param net_vertices chosen vertices to plot the heatmap
#' @param net_edges chosen edges to plot the heatmap, can be either the edge id's or its node endpoints (e.j. c(1,2, 2,3, 7,8))
#' @param show_events option to show the events as orange squares, FALSE by default
#' @param alpha optional argument to set the transparency of the events (show_events = TRUE). The range is from 0.1 (transparent) to 1 (opaque). Default: alpha = 1
#' @param ... extra arguments for the class ggplot
#' 
#' @return The plot of the heatmap with class c("gg", "ggplot")
#' 
#' @examples
#' 
#' \dontrun{
#' data("und_intnet_chicago")
#' PlotHeatmap(und_intnet_chicago, heat_type='moran')
#' }
#' 
#' @export
PlotHeatmap.intensitynet <- function(obj, heat_type = 'none', intensity_type = 'none', net_vertices = NULL, net_edges = NULL, show_events = FALSE, alpha = 1, ...){
  g <- obj$graph
  adj_mtx <- igraph::as_adj(graph = g)
  adj_listw <- spdep::mat2listw(adj_mtx)
  nb <- adj_listw$neighbours
  w_listw <- spdep::nb2listw(nb, style = "W",  zero.policy=TRUE)
  
  if(heat_type != 'none' && heat_type != 'moran' && heat_type != 'geary' && 
     heat_type != 'v_intensity' && heat_type != 'e_intensity'){
    
    if( !(heat_type %in% igraph::edge_attr_names(g) ))
    {
      warning('Parameter "heat_type" should be; for correlations: "moran" or "geary", 
                                               for intensities: "v_intensity" or "e_intensity", 
                                               for marks: the name of the mark. 
                                               Using default ("none").')
      heat_type <- 'none'
    }
  }
  
  if (is.null(net_vertices) && is.null(net_edges)){
    net_vertices <- igraph::V(g)
    net_edges <- igraph::E(g)
  } else if (!is.null(net_vertices)){
    net_vertices <- igraph::V(g)[net_vertices] # Convert to class 'igraph.v'
  }else if(class(net_edges) != 'igraph.es'){
    net_edges <- igraph::E(g, P = net_edges) # Convert to class 'igraph.es'
  }
  
  # If the intensity is not provided, try to take it from the given network
  if( is.null(igraph::vertex_attr(graph = g, name = intensity_type)) ){
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
  
  if(heat_type == 'moran'){ # Local Moran-i
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
    
  }else if(heat_type == 'geary'){  # Local Geary-c
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
    
  }else if(heat_type == 'getis'){ # Local Getis-G*
    message("Needs implementation")
    # TODO: Implement Getis G.
    
    # locgg <- spdep::localG(x = intensity, listw = b_listw)
    # 
    # # create a new variable identifying the moran plot quadrant for each observation, dismissing the non-significant ones
    # sig_dstr <- NA
    # 
    # data_df <- data.frame(xcoord = node_coords$xcoord, 
    #                       ycoord = node_coords$ycoord, 
    #                       value = sig_dstr)
    
    
    
  }else if(heat_type == 'v_intensity' || heat_type == 'e_intensity'){
    norm_int <- (intensity - min(intensity)) / (max(intensity) - min(intensity))
    data_df <- data.frame(xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord, 
                          value = norm_int)
    
  }else{
    data_df <- data.frame(xcoord = node_coords$xcoord, 
                          ycoord = node_coords$ycoord, 
                          value = NA)
  }
  geoplot_obj <- list(intnet = obj, 
                      data_df = data_df, 
                      net_vertices = net_vertices, 
                      net_edges = net_edges,
                      mode = heat_type, 
                      show_events = show_events,
                      alpha = alpha)
  class(geoplot_obj) <- "netTools"
  
  return( GeoreferencedGgplot2(geoplot_obj, ...) )
}


#' Plot the net and the events in the neighborhood area of the given node
#' 
#' @name PlotNeighborhood.intensitynet
#' 
#' @param obj intensitynet object
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
                 distances_mtx = obj$distances_mtx,
                 event_correction = obj$event_correction)
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
                 distances_mtx = sub_dist_mtx,
                 event_correction = obj$event_correction)
  attr(intnet, 'class') <- class(obj)
  
  intnet
}


#' Calculates the shortest distance path between two nodes (based on the minimum amount of edges).
#' The function also returns the total weight of the path, if the weight is not available, returns
#' the number of edges.
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