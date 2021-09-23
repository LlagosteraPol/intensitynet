require(spatstat)
require(igraph)
require(intervals)
require(aoos)
require(roxygen2)

setOldClass(c("igraph"))


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
      
      weighted_adj_mtx <- adjacency_mtx * .distances_mtx
      
      .g <<- graph_from_adjacency_matrix(weighted_adj_mtx, mode = graph_type, weighted = TRUE)
      
      .setGraphCoords(x_coord_node, y_coord_node)
    },
    
    #' Gives the shortest path between two nodes. If the graph is weighted, takes the weight into account.  
    #' 
    #' @name NetintensityV2_shortestDistance
    #' 
    #' @param node_id1 Starting node ID
    #' @param node_id2 Final node ID
    #' 
    #' @return a list with the vertices of the shortest path (path) and its total weight (weight)
    #' 
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
        plot(.g, layout = norm_coords, vertex.label=NA, vertex.size=2, window=TRUE, axes=TRUE, edge.label = edge_attr(.g)$intensity, edge.label.cex = 0.5)
      }
      else{
        plot(.g, vertex.label=NA, vertex.size=2,vertex.size2=2)
      }
    }
  )
)