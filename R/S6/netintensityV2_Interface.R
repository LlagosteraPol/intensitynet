source("netintensityV2.R")
source("netintensityV2_UndirectedFactory.R")
source("netintensityV2_DirectedFactory.R")


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
