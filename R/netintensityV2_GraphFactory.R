source("netintensityV2.R")

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