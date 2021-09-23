source("netintensityV2.R")

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