
rm(list = ls())

#Set working directory
setwd("R")

source("netintensityV2_UndirectedFactory.R")

# Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
load("..//Castellon.RData")

# Node coordinates: Georeferenced coordinates from 'castellon' nodes
load("../Data/nodes.RData")

# Event (crime coordinates)
load("../Data/crimes.RData")

#subset of events
crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)

undirected <- UndirectedFactory()
undirected_general <- undirected$constructGraph(Castellon, "general", nodes$cx, nodes$cy, crim$X, crim$Y)

# Calculates the intensity of the given node ID
node_intensity <- undirected_general$nodeIntensity(252)
node_intensity

# Calculates the mean intensity of the given node ID
mean_node_intensity <- undirected_general$meanNodeIntensity(252)
mean_node_intensity

# Calculates the intensity of the given edge
edge_intensity <- undirected_general$edgeIntensity(252, 248)
edge_intensity

# Calculates the intensity of the given path
path_intensity <- undirected_general$pathIntensity(c(252, 248, 250, 246, 254, 249, 242, 235, 230))
path_intensity

# Calculates the intensity of the path with start at Node_ID: 253 and end at Node_ID: 230. Takes into 
# account the distances of the edges.
shortest_path_intensity <- undirected_general$shortestPathIntensity(252, 230, TRUE)
shortest_path_intensity

# Calculates all the intensities of the graph
start_time <- Sys.time()
intensities <- undirected_general$calculateAllIntensities()
end_time <- Sys.time()

time1 <- end_time - start_time
time1

# Checking the attributes of the nodes and edges of the graph
g <- undirected_general$getGraph()
vertex_attr(g)
edge_attr(g)
