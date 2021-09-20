source("netintensityV2_GraphFactory.R")
source("netintensityV2_UndirectedGeneral.R")
source("netintensityV2_UndirectedPlannar.R")

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