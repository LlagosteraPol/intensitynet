rm(list = ls())

#Set working directory
setwd("R")

source("S3/main.R")

# ---------------------------------------------DATA LOADING----------------------------------------------------
# Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
load("../Data/Castellon.RData")

# Node coordinates: Georeferenced coordinates from 'castellon' nodes
load("../Data/nodes.RData")

# Event (crime coordinates)
load("../Data/crimes.RData")

#subset of events
crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)

# --------------------------------------- INIT NETINTENSITY CLASS----------------------------------------------
intnet_und <- intensitynet(Castellon, nodes, crim)

Castellon_obj <-list(mtx = Castellon)
class(Castellon_obj) <- "netTools"
dirCastellon <-  Undirected2RandomDirectedAdjMtx(Castellon_obj)
intnet_dir <- intensitynet(dirCastellon, nodes, crim, graph_type='directed')

intnet_mix <- intensitynet(dirCastellon, nodes, crim, graph_type='mixed')

# ---------------------------------- NET TOOLS CLASS: FUNCTION TESTING ----------------------------------------

# Choose the type of graph for testing
#intnet <- intnet_und
intnet <- intnet_dir
#intnet <- intnet_mix

class(intnet)

# DISTANCES MATRIX
dist_mtx_test <- intnet$distances

# SHORTEST DISTANCE
short_dist_obj <- list(graph=intnet$graph, node_id1 = 'V601', node_id2 = 'V701', distances_mtx = intnet$distances)
class(short_dist_obj) <- "netTools"
short_dist <- ShortestDistance(short_dist_obj)

# GEORREFERENCED PLOT
plot(intnet)

# NODE INTENSITY
node_int <- MeanNodeIntensity(intnet, 'V601') 

# EDGE INTENSITY
edge_int <- EdgeIntensity(intnet, 'V1', 'V9')

# PATH INTENSITY
int_path <- PathIntensity(intnet, short_dist$path)

# All intensities
intnet_all <- CalculateEventIntensities(intnet)
g <- intnet_all$graph
edge_attr_names(g)
vertex_attr_names(g)

if(intnet_all$graph_type == 'undirected'){
  vertex_attr(g, 'intensity', V(g)['V1']) 
  
  for(node_id in V(g)){
    if(V(g)[node_id]$intensity>0) cat(node_id,": ",V(g)[node_id]$intensity, "\n")
  }
} else{
  vertex_attr(g, 'intensity_in', V(g)['V1']) 
  vertex_attr(g, 'intensity_out', V(g)['V1']) 
  
  for(node_id in V(g)){
    if(V(g)[node_id]$intensity_in>0) cat(node_id,": ",V(g)[node_id]$intensity_in, "\n")
  }
  
  for(node_id in V(g)){
    if(V(g)[node_id]$intensity_out>0) cat(node_id,": ",V(g)[node_id]$intensity_out, "\n")
  }
  
  if(intnet_all$graph_type == 'mixed'){
    vertex_attr(g, 'intensity_und', V(g)['V1']) 
    vertex_attr(g, 'intensity_all', V(g)['V1']) 
    
    for(node_id in V(g)){
      if(V(g)[node_id]$intensity_und>0) cat(node_id,": ",V(g)[node_id]$intensity_und, "\n")
    }
    
    for(node_id in V(g)){
      if(V(g)[node_id]$intensity_all>0) cat(node_id,": ",V(g)[node_id]$intensity_all, "\n")
    }
  }
}

for(edge_id in E(g)){
  if(E(g)[edge_id]$intensity>0) print(E(g)[edge_id]$intensity)
}

#-----------------------------INTENSITYNET CLASS: FUNCTION TESTING---------------------------------


correlations <- EventCorrelation(intnet_all, 'correlation', 2)

locmoran <- NodeLocalCorrelation(intnet_all, 'moran')
locg <- NodeLocalCorrelation(intnet_all, 'g')




