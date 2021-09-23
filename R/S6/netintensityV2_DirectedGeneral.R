source("netintensityV2.R")

#' Subclass of NetintensityV2 specific to work with Directed Plannar graphs
#' 
DirectedGeneral <- setRefClass(
  Class = "DirectedGeneral",
  contains = "NetintensityV2",
  methods = list(
    
    #' Test function to know which class is refering
    #'
    #' @name DirectedPlannar_what
    #'
    what = function(){
      print("Directed General Graph")
    },
    
    #' Calculates edgewise and mean nodewise intensity function for Directed Plannar networks
    #' 
    #' @name DirectedGeneral_calculateIntensities
    #' 
    calculateIntensities = function(){
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
        
        
        
      }
    }
  )
)