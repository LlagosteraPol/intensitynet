source("netintensityV2_GraphFactory.R")
source("netintensityV2_DirectedPlannar.R")

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