source("netintensityV2.R")

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