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
intnet <- intensitynet(Castellon, nodes, crim)

# ---------------------------------- NET TOOLS CLASS: FUNCTION TESTING ----------------------------------------

# DISTANCES MATRIX
dist_mtx_test <- intnet$distances

# SHORTEST DISTANCE
short_dist_obj <- list(graph=intnet$graph, node1 = 'V1', node2 = 'V10', distances_mtx = intnet$distances)
class(short_dist_obj) <- "netTools"
short_dist <- shortestDistance(short_dist_obj)

# GEORREFERENCED PLOT
plot(intnet)

# ---------------------------- INTENSITYNET UNDIRECTED CLASS: FUNCTION TESTING --------------------------------
# NODE INTENSITY
node_int <- nodeIntensity(intnet, 'V1')

# EDGE INTENSITY
edge_int <- edgeIntensity(intnet, 'V1', 'V9')

# PATH INTENSITY
int_path <- pathIntensity(intnet, short_dist$path)

# All intensities
intnet_all <- calculateEventIntensities(intnet)
g <- intnet_all$graph
edge_attr_names(g)
vertex_attr_names(g)

for(node_id in V(g)){
  if(V(g)[node_id]$intensity>0) print(V(g)[node_id]$intensity)
}

for(edge_id in E(g)){
  if(E(g)[edge_id]$intensity>0) print(E(g)[edge_id]$intensity)
}










