source("netintensityV2.R")

#' Subclass of NetintensityV2 specific to work with Undirected General graphs
#' 
UndirectedGeneral <- setRefClass(
  Class = "UndirectedGeneral",
  contains = "NetintensityV2",
  fields = list(),
  methods = list(
    
    #' Test function to know which class is refering
    #'
    #' @name UndirectedGeneral_what
    #'
    what = function(){
      print("Undirected General Graph Subclass")
    },
    
    
    #' Give the intensity of the given path
    #' 
    #' @name UndirectedGeneral_pathIntensity
    #' 
    #' @param path_nodes List of node ID's that form the path.
    #' 
    #' @return path_intensity - intensity of the given path
    #' 
    pathIntensity = function(path_nodes){
      
      edge_counts <- list()
      path_intensity <- 0
      
      prev <- NULL
      for(node_id in path_nodes){
        if(is.null(prev)){
          prev <- node_id
          next
        }
        
        path_intensity <- path_intensity + Reduce('+', edgeIntensity(prev, node_id))
        
        prev <- node_id
      }
      # Divide the intensity of the edges by their number (In a path -> N edges = N vertices - 1)
      path_intensity <- path_intensity / (length(path_nodes) - 1)
      
      return(path_intensity)
    },
    
    
    #' Give the intensity of the shortest path between two nodes
    #' 
    #' @name UndirectedGeneral_shortestPathIntensity
    #' 
    #' @param node_id1 Starting node ID
    #' @param node_id2 Final node ID
    #' @param wheighted Indicates if takes in account the wheight (TRUE) of the edges or not (FALSE).
    #' 
    #' @return a list with the intensity of the shortest path (intensity) and the path vertices (path)
    #' 
    shortestPathIntensity = function(node_id1, node_id2, weighted = FALSE){
      if(weighted){
        path <- shortestDistance(node_id1, node_id2)$path
      }
      else{
        path <- unlist(get.shortest.paths(.g, node_id1, node_id2)$vpath)
      }
      
      return(list(intensity = pathIntensity(path), path = path))
    },
    
    
    #' Calculates the intensity of the given node and input into the node attribute of the graph.
    #' 
    #' @name UndirectedGeneral_nodeIntensity
    #' 
    #' @param node_id ID of the node
    #' 
    #' @return list with the total intensity of the node ('intensity') and 
    #' the intensity respect each neighbor ('detailed') 
    #' 
    nodeIntensity = function(node_id){
      
      if(degree(.g, node_id) > 0){
        neighbors_list <- neighbors(.g, node_id)
        
        ev_mat <- matrix(0, ncol = length(neighbors_list)) 
        colnames(ev_mat) <- as.vector(neighbors_list) 
        rownames(ev_mat) <- node_id
        
        for (neighbor_id in neighbors_list){
          ev_mat[as.character(node_id), as.character(neighbor_id)] <- edgeIntensity(node_id, neighbor_id)
        }
        
        total_intensity <- Reduce('+', ev_mat)
        
        .setNodeIntensity(node_id, total_intensity)
        
        return(list(intensity = total_intensity, detailed = ev_mat))
      }
    },
    
    
    #' Gives the mean node intensity of the given node ID
    #'
    #' @name UndirectedGeneral_meanNodeIntensity
    #' 
    #' @param node_id ID of the node
    #' 
    #' @return Mean intensity of the node
    #' 
    meanNodeIntensity = function(node_id){
      Reduce('+', nodeIntensity(node_id)$intensity)/degree(.g, node_id)
    },
    
    
    #' If not calculated, calculates the intesnity of the edge with nodes; node_id1, node_id2 and
    #' input into the edge attribute of the graph. If the edge already contains an intensity,
    #' gives it directly.
    #'
    #' @name UndirectedGeneral_edgeIntensity
    #' 
    #' @param node_id1 First node ID of the edge
    #' @param node_id2 Second node ID of the edge
    #' 
    #' @return edge_intensity - Intensity of the edge
    #'
    edgeIntensity = function(node_id1, node_id2){
      
      # Note that the igraph library already handle the error when one of the node id's 
      # are not part of the graph and gives the proper information about it.
      edge_id <- get.edge.ids(.g, c(node_id1,node_id2))
      # If the intensity of this edge was previously calculated, then return it
      if(edge_id != 0 & !is.null(vertex_attr(.g, "intensity", edge_id)) ){
        if(!is.na(vertex_attr(.g, "intensity", edge_id))){
          return(vertex_attr(.g, "intensity", edge_id))
        }
      }
      
      # Distance between the node and its neighbor
      res <- tryCatch(
        {
          abs(.distances_mtx[node_id1, node_id2]) # Distance between the node and its neighbor 
        },
        # If the nodes are not part of the graph, give the proper information of the error.
        error=function(cond) {
          if( is.null(match(node_id1, V(.g))) ){
            message("First vertice ID (node_id1) doesn't exist in the graph.")
          }
          else if( is.null(match(node_id2, V(.g))) ){
            message("Second vertice ID (v2) doesn't exist in the graph.")
          }
          else if(is.null(match(node_id1, V(.g))) & is.null(match(node_id2, V(.g)))){
            message("First and second vertices (node_id1, node_id2) doesn't exist in the graph.")
          }
          else{
            neighbors_list <- neighbors(.g, node_id1)
            
            if(! node_id2 %in% neighbors_list){
              message("Second vertice (node_id2) it's not a neighbor of first vertice (node_id1)")
            }
            else{
              message(cond)
            }
          }
        }
      )    
      
      ev_mat <- matrix(0, ncol = 1) 
      colnames(ev_mat) <- node_id2 
      rownames(ev_mat) <- node_id1
      
      node1_coords <- list(xcoord = vertex_attr(.g, "xcoord", node_id1),
                           ycoord = vertex_attr(.g, "ycoord", node_id1))
      
      node2_coords <- list(xcoord = vertex_attr(.g, "xcoord", node_id2),
                           ycoord = vertex_attr(.g, "ycoord", node_id2))
      
      # Count events inside a window formed by the node and its neighbor
      
      # Defining event window
      x_min <- min(c(node1_coords$xcoord, node2_coords$xcoord))
      x_max <- max(c(node1_coords$xcoord, node2_coords$xcoord))
      
      y_min <- min(c(node1_coords$ycoord, node2_coords$ycoord))
      y_max <- max(c(node1_coords$ycoord, node2_coords$ycoord))
      
      indicator <- 0
      # Counting events
      for(row in 1:nrow(.events_mtx)) {
        event_x <- .events_mtx[row, 1]
        event_y <- .events_mtx[row, 2]
        
        if(x_min <= event_x & x_max >= event_x & y_min <= event_y & y_max >= event_y){
          indicator <- indicator + 1
        } 
      }
      edge_intensity <- indicator/res
      
      .setEdgeIntensity(node_id1, node_id2, edge_intensity)
      
      return(edge_intensity)
    },
    
    
    #' Calculates edgewise and mean nodewise intensity function for Undirected networks
    #' 
    #' @name UndirectedGeneral_calculateIntensities
    #' 
    #' @return mean_intensity - vector containing the mean intensity per node
    #' @return edge_intensity - array containing edge intensity for adjacent edges
    #' 
    calculateAllIntensities = function(){
      
      edge_counts <- list()
      #counts <- list()
      counts <- c()
      
      #TODO: Consider progressive calculation with the other methods
      # check if the intensities was previously calculated, if not, calculate them
      if(length(.intensities) == 0){
        for(node_id in V(.g)){
          if(degree(.g, node_id) > 0){
            
            ev_mat <- nodeIntensity(node_id)$detailed
            
            #Adds result of Edgewise intenisty function to 'edge_counts'
            edge_counts[[node_id]] <- ev_mat   
            #Adds result of Nodewise mean intenisty function to 'counts'
            counts[[node_id]] <- Reduce('+', ev_mat)/degree(.g, node_id)
            
            
            if(length(counts[[node_id]]) == 0){
              counts[[node_id]] <- 0
            }
          }
          
          # Counts for isolated nodes
          else{
            counts[[node_id]] <- 0
          }
          cat("Calculation for node", node_id, "\n")
        }
        names(counts) = vertex_attr(.g)$name
        .intensities <<- list(mean_intensity = counts, edge_intensity = edge_counts)
      }
      return(.intensities)
    },
    
    
    #' Calculates edgewise and mean nodewise intensity function for Undirected networks
    #' 
    #' @name UndirectedGeneral_calculateIntensities_old
    #' 
    #' @return counts vector containing the mean.intensity per node
    #' @return edge.counts array containing edge.intensity for adjacent edges
    #' @return g underlying graph object # TODO: Delete, can be recovered using the getter method 
    #' @return deg Vector of degrees per node # TODO: Delete, can be recovered using 'degree(igraph)'
    #' 
    calculateIntensities_old = function(){
      
      edge_counts <- list()
      counts <- list()
      
      for(node_id in V(.g)){
        if(degree(.g, node_id) > 0){
          neighbors_list <- neighbors(.g, node_id)
          
          ev_mat <- matrix(0, ncol = length(neighbors_list)) 
          colnames(ev_mat) <- as.vector(neighbors_list) 
          rownames(ev_mat) <- node_id
          
          for (neighbor_id in neighbors_list){
            
            res <- abs(.distances_mtx[node_id, neighbor_id]) # Distance between the node and its neighbor 
            
            node_coords <- list(xcoord = vertex_attr(.g, "xcoord", node_id),
                                ycoord = vertex_attr(.g, "ycoord", node_id))
            
            nei_coords <- list(xcoord = vertex_attr(.g, "xcoord", neighbor_id),
                               ycoord = vertex_attr(.g, "ycoord", neighbor_id))
            
            # Count events inside a window formed by the node and its neighbor
            
            # Defining event window
            x_min <- min(c(node_coords$xcoord, nei_coords$xcoord))
            x_max <- max(c(node_coords$xcoord, nei_coords$xcoord))
            
            y_min <- min(c(node_coords$ycoord, nei_coords$ycoord))
            y_max <- max(c(node_coords$ycoord, nei_coords$ycoord))
            
            
            indicator <- 0
            # Counting events
            for(row in 1:nrow(.events_mtx)) {
              event_x <- .events_mtx[row, 1]
              event_y <- .events_mtx[row, 2]
              
              if(x_min <= event_x & x_max >= event_x & y_min <= event_y & y_max >= event_y){
                indicator <- indicator + 1
              } 
            }
            ev_mat[as.character(node_id), as.character(neighbor_id)] <- indicator/res
          }
          
          #Adds result of Edgewise intenisty function to 'edge_counts'
          edge_counts[[node_id]] <- ev_mat   
          #Adds result of Nodewise mean intenisty function to 'counts'
          counts[[node_id]] <- Reduce('+', ev_mat)/degree(.g, node_id)
          
          
          if(length(counts[[node_id]]) == 0){
            counts[[node_id]] <- 0
          }
        }
        
        # Counts for isolated nodes
        else{
          counts[[node_id]] <- 0
        }
        cat("Calculation for node", node_id, "\n")
      }
      
      # TODO: Remove 'graph' and 'n_neighbors'. Can be recovered by the functions of the class.
      list(mean_intensity = counts, edge_intensity = edge_counts, graph = .g, n_neighbors = degree(.g))
    }
  )
)