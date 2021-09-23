require(spatstat)
require(igraph)
require(intervals)
require(aoos)
require(roxygen2)

setOldClass(c("igraph"))



#-----------------------------OBJECT CLASSES------------------------------
#' Reference class that encapsulates different function to calculate
#' intensities
#' 
#' @field .g (Private) igraph object
#' @field .distances_mtx (Private) Matrix containing all the distances 
#' between each pair of nodes
#' 
NetintensityV2 <- setRefClass(
  Class = "NetintensityV2",
  fields = list(.g = "igraph", 
                .distances_mtx = "matrix",
                .events_mtx = "matrix",
                .intensities = "list"),
  methods = list(
    
    
    #' Create a weighted graph (igraph object) wich their vertices contains
    #' coordinates as attributes. The resulting graph is stored in the class.
    #' 
    #' @name NetintensityV2_initGraph
    #' @param adjacency_mtx An adjacency matrix encoding network structure
    #' @param graph_type (character) Graph mode used by igraph.
    #' @param x_coord_node X-coordinates of node in network
    #' @param y_coord_node Y-coordinates of node in network
    #' @param x.event X-coordinates of event
    #' @param y.event Y-coordinates of event
    #' 
    initGraph = function(adjacency_mtx, graph_type, x_coord_node, y_coord_node, x_event, y_event){
      
      .events_mtx <<- matrix(cbind(x_event, y_event), ncol=2)
      
      .setDistancesMtx(x_coord_node, y_coord_node)
      
      wheighted_adj_mtx <- adjacency_mtx * .distances_mtx
      
      .g <<- graph_from_adjacency_matrix(wheighted_adj_mtx, mode = graph_type, weighted = TRUE)
      
      .setGraphCoords(x_coord_node, y_coord_node)
    },
    
    shortestDistance = function(node_id1, node_id2){
      weighted_path <- unlist(get.shortest.paths(.g, node_id1, node_id2)$vpath)
      if(!is.null(.distances_mtx)){
        weight_sum <- sum(E(.g, path = unlist(weighted_path))$weight)
      }
      else{
        weight_sum <- length(weighted_path)
      }
      list(path = weighted_path, weight = weight_sum)  
    },
    
    
    #' Private function. Set coordinates as attributes to the vertices 
    #' of the graph stored in the class.
    #' 
    #' @name NetintensityV2_.setGraphCoords
    #' 
    #' @param x_coord_node X-coordinates of node in network
    #' @param y_coord_node Y-coordinates of node in network
    #' 
    .setGraphCoords = function(x_coord_node, y_coord_node){
      # TODO: change x and y to long (longitude) and lat (latitude), beware of the coordinate system type (e.g WGS84)
      .g <<- .g %>% set_vertex_attr(name = "xcoord", value = x_coord_node) %>% 
        set_vertex_attr(name = "ycoord", value = y_coord_node)
    },
    
    
    #' Private function. Set the distances matrix of all pair of nodes. 
    #' The matrix is stored in the class.
    #' 
    #' @name NetintensityV2_.setDistancesMtx
    #' 
    #' @param x_coord_node X-coordinates of node in network
    #' @param y_coord_node Y-coordinates of node in network
    #' 
    .setDistancesMtx = function(x_coord_node, y_coord_node){
      .distances_mtx <<- pairdist(ppp(x_coord_node,
                                      y_coord_node,
                                      xrange=c(min(as.numeric(x_coord_node)), max(as.numeric(x_coord_node))),
                                      yrange=c(min(as.numeric(y_coord_node)), max(as.numeric(y_coord_node)))))
    },
    
    #' Private function. Set the given value as an "intensity" attribute of an edge with nodes; node_id1, node_id2.
    #' 
    #' @name NetintensityV2_.setEdgeIntensity
    #' 
    #' @param node_id1 First node ID of the edge
    #' @param node_id2 Second node ID of the edge
    #' @param value Value to set as intensity
    #' 
    .setEdgeIntensity = function(node_id1, node_id2, value){
      edge_id <- get.edge.ids(.g, c(node_id1, node_id2))
      .g <<- .g %>% set_edge_attr(name = "intensity", index = edge_id, value = value)
    },
    
    #' Private function. Set the given value as an "intensity" attribute of the given node ID.
    #' 
    #' @name NetintensityV2_.setNodeIntensity
    #' 
    #' @param node_id ID of the node
    #' @param value Value to set as intensity
    #' 
    .setNodeIntensity = function(node_id, value){
      .g <<- .g %>% set_vertex_attr(name = "intensity", index = node_id, value = value)
    },
    
    
    #' Getter of the graph stored in the class.
    #'
    #' @name NetintensityV2_getGraph
    #' 
    #' @return .g - the graph (igraph object)
    #' 
    getGraph = function(){
      .g
    },
    
    
    #' Abstract function that must be implemented by subclasses.
    #' 
    #' @name NetintensityV2_calculateIntensities
    #'  
    calculateIntensities = function(){
      # implemented by subclasses
    },
    
    
    #' Plot the graph using the coordinates if provided. If not, perform a normal plot.
    #' 
    #' @name NetintensityV2_georeferencedPlot
    #'
    georeferencedPlot = function(){
      if(!is.null(.distances_mtx)){
        norm_coords = layout.norm(matrix(cbind(vertex_attr(.g)$xcoord, vertex_attr(.g)$ycoord), ncol=2))
        plot(.g, layout = norm_coords, vertex.label=NA, vertex.size=2, window=TRUE, axes=TRUE, edge.label = edge_attr(.g)$weight, edge.label.cex = 0.5)
      }
      else{
        plot(.g, vertex.label=NA, vertex.size=2,vertex.size2=2)
      }
    }
  )
)

# Subclass of NetintensityV2

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

#' Subclass of NetintensityV2 specific to work with Undirected Plannar graphs
#' 
UndirectedPlannar <- setRefClass(
  Class = "UndirectedPlannar",
  contains = "NetintensityV2",
  methods = list(
    
    #' Test function to know which class is refering
    #'
    #' @name UndirectedPlannar_what
    #'
    what = function(){
      print("Undirected Plannar Graph")
    },
    
    
    #' Calculates edgewise and mean nodewise intensity function for Undirected Plannar networks
    #' 
    #' @name UndirectedPlannar_calculateIntensities
    #' 
    calculateIntensities = function(){
      # must be implemented
    }
  )
)


#' Subclass of NetintensityV2 specific to work with Directed Plannar graphs
#' 
DirectedPlannar <- setRefClass(
  Class = "DirectedPlannar",
  contains = "NetintensityV2",
  methods = list(
    
    #' Test function to know which class is refering
    #'
    #' @name DirectedPlannar_what
    #'
    what = function(){
      print("Directed Plannar Graph")
    },
    
    #' Calculates edgewise and mean nodewise intensity function for Directed Plannar networks
    #' 
    #' @name DirectedPlannar_calculateIntensities
    #' 
    calculateIntensities = function(){
      # must be implemented
    }
  )
)

#-----------------------------FACTORY METHODS------------------------------


#' Reference class providing a Factory interface for NetintensityV2 class
#' 
NetintensityV2GraphFactory <- setRefClass(
  Class = "NetintensityV2GraphFactory",
  fields = list(),
  methods = list(
    
    
    #' Abstract function to implement by the class subclasses
    #'
    #' @name NetintensityV2GraphFactory_initSubtype
    #'
    #' @param adjacency_mtx An adjacency matrix encoding network structure
    #' @param graph_characteristic (character) Graph structure, "General", "Plannar"....
    #' @param x_coord_node X-coordinates of node in network
    #' @param y_coord_node Y-coordinates of node in network
    #' @param x.event X-coordinates of event
    #' @param y.event Y-coordinates of event
    #' 
    initSubtype = function(adjacency_mtx, graph_characteristic, x_coord_node, y_coord_node, x_event, y_event){
      # implemented by subclasses
    },
    
    
    #' Calculate intensities for the correct graph type and structure
    #'
    #' @name NetintensityV2GraphFactory_constructGraph
    #'
    #' @param adjacency_mtx An adjacency matrix encoding network structure
    #' @param graph_characteristic (character) Graph structure, "General", "Plannar"....
    #' @param x_coord_node X-coordinates of node in network
    #' @param y_coord_node Y-coordinates of node in network
    #' @param x.event X-coordinates of event
    #' @param y.event Y-coordinates of event
    #' 
    #' @return netint - NetintensityV2 proper subclass object 
    #' 
    constructGraph = function(adjacency_mtx, graph_characteristic, x_coord_node, y_coord_node, x_event, y_event)
    {
      netint <- initSubtype(adjacency_mtx, graph_characteristic, x_coord_node, y_coord_node, x_event, y_event)
      netint
    }
  )
)

UndirectedFactory <- setRefClass(
  Class = "UndirectedFactory",
  contains = "NetintensityV2GraphFactory",
  fields = list(),
  methods = list(
    
    
    #' Initializes an NetintensityV2 object based on Undirected graphs and its proper structure.
    #'
    #' @name UndirectedFactory_initSubtype
    #' 
    #' @param adjacency_mtx An adjacency matrix encoding network structure
    #' @param graph_characteristic (character) Graph structure, "General", "Plannar"....
    #' @param x_coord_node X-coordinates of node in network
    #' @param y_coord_node Y-coordinates of node in network
    #' @param x.event X-coordinates of event
    #' @param y.event Y-coordinates of event
    #' 
    #' @return undirPln - UndirectedPlannar object
    #'
    initSubtype = function(adjacency_mtx, graph_characteristic, x_coord_node, y_coord_node, x_event, y_event){
      
      if(tolower(graph_characteristic) == "general"){
        undirGnr <- UndirectedGeneral()
        undirGnr$initGraph(as.matrix(adjacency_mtx), "undirected", x_coord_node, y_coord_node, x_event, y_event)
        return(undirGnr)
      }
      
      else if(tolower(graph_characteristic) == "plannar"){
        undirPln <- UndirectedPlannar()
        undirPln$what()
        return(undirPln)
      }
    }
  )
)

DirectedFactory <- setRefClass(
  Class = "DirectedFactory",
  contains = "NetintensityV2GraphFactory",
  fields = list(),
  methods = list(
    
    
    #' Initializes an NetintensityV2 object based on Directed graphs and its proper structure.
    #'
    #' @name DirectedFactory_initSubtype
    #' 
    #' @param adjacency_mtx An adjacency matrix encoding network structure
    #' @param graph_characteristic (character) Graph structure, "General", "Plannar"....
    #' @param x_coord_node X-coordinates of node in network
    #' @param y_coord_node Y-coordinates of node in network
    #' @param x.event X-coordinates of event
    #' @param y.event Y-coordinates of event
    #' 
    #' @return dirPln - DirectedPlannar object
    #'
    initSubtype = function(adjacency_mtx, graph_characteristic, x_coord_node, y_coord_node, x_event, y_event){
      
      if(tolower(graph_characteristic) == "plannar"){
        dirPln <- directedPlannar()
        dirPln$what()
        return(dirPln)
      }
    }
  )
)

#' Interface for NetIntensity
#' 
Netintensity <- setRefClass(
  Class = "Netintensity",
  #contains = "NetintensityV2",
  fields = list(adjacency_mtx = "matrix", 
                x_coord_node = "numeric", 
                y_coord_node = "numeric", 
                x_event = "numeric", 
                y_event = "numeric",
                graph_type = "character",
                graph_characteristic = "character",
                .graph_structure = "NetintensityV2"),
  methods = list(
    
    initialize = function(adjacency_mtx,
                          x_coord_node, 
                          y_coord_node, 
                          x_event, 
                          y_event,
                          graph_type = "undirected", 
                          graph_characteristic = "general"){
      
      if(missing(x_coord_node) || missing(y_coord_node) || missing(x_event) || missing(y_event)) {
        stop("Values are missing")
      }
      
      factory = NULL
      if(tolower(graph_type) == "undirected"){
        factory <- UndirectedFactory()
      }
      
      else if(tolower(graph_type) == "directed"){
        factory <- DirectedFactory()
      }
      
      # else if(tolower(graph_type) == "mixed"){
      #   # This kind of class is not yet implemented
      #   factory <- MixedFactory()
      # }
      
      else{
        stop("Graph type is not correct")
      }
      
      callSuper(graph_type = graph_type, 
                graph_characteristic = graph_characteristic,
                x_coord_node = x_coord_node, 
                y_coord_node = y_coord_node, 
                x_event = x_event, 
                y_event = y_event,
                .graph_structure = factory$constructGraph(adjacency_mtx, 
                                                          graph_characteristic, 
                                                          x_coord_node, 
                                                          y_coord_node, 
                                                          x_event, 
                                                          y_event))
    },
    
    summary = function(){
      data <- .graph_structure$calculateAllIntensities()
      data <- append(data, list(range = range( data[["mean_intensity"]] ) ) )
      return(data)
    },
    
    plot = function(tp = "hist"){
      
      if(tolower(tp) == "hist"){
        hist(degree(.graph_structure$getGraph()), 
             main = "",
             xlab = "degree")
      }
      
      else if (tp == "graph"){
        .graph_structure$georeferencedPlot()
      }
    },
    
    getGraphStructure = function(){
      return(.graph_structure)
    }
    
  )
)

#----------------------------- RUN ------------------------------
# Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
# data(Castellon) 
# 
# # Node coordinates: Georeferenced coordinates from 'castellon' nodes
# data(nodes) 
# 
# data(crimes) # Event (crime coordinates)



#subset of events

#crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)


# undirected <- UndirectedFactory()
# undirected_general <- undirected$constructGraph(Castellon, "general", nodes$cx, nodes$cy, crim$X, crim$Y)
# 
# 
# start_time <- Sys.time()
# intensities <- undirected_general$calculateAllIntensities()
# end_time <- Sys.time()
# 
# time1 <- end_time - start_time
# g <- undirected_general$getGraph()
# 
# start_time <- Sys.time()
# intensities <- undirected_general$calculateIntensities_old()
# end_time <- Sys.time()
# 
# time2 <- end_time - start_time


#edge_intensity <- undirected_general$edgeIntensity(252, 248)
# node_intensity <- undirected_general$nodeIntensity(252)
# mean_node_intensity <- undirected_general$meanNodeIntensity(252)
# 
#path_intensity <- undirected_general$pathIntensity(c(252, 248))
# shortest_path_intensity <- undirected_general$shortestPathIntensity(252, 248, TRUE)
# 

net <- Netintensity(Castellon, nodes$cx, nodes$cy, crim$X, crim$Y)
data <- net$summary()
#net$plot()

