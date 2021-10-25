rm(list = ls())

#Set working directory
setwd("R/S3/")

source("./main.R")

# ---------------------------------------------DATA LOADING----------------------------------------------------
# Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
load("../../Data/Castellon.RData")

# Node coordinates: Georeferenced coordinates from 'castellon' nodes
load("../../Data/nodes.RData")

# Event (crime coordinates)
load("../../Data/crimes.RData")

#subset of events
crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)
#crim <- crimes
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

  pdf("S3/Plots/area_with_grid_plot.pdf")
  plot(intnet_all, enable_grid = TRUE, axis=TRUE)
  dev.off()
  
  plot(intnet_all, node_label = 'intensity', edge_label='none', vertex.color='red')
  
  vertex_attr(g, 'intensity', V(g)['V1']) 
  
  for(node_id in V(g)){
    if(V(g)[node_id]$intensity>0) cat(node_id,": ",V(g)[node_id]$intensity, "\n")
  }
  
  correlations <- NodeGeneralCorrelation(intnet_all, dep_type = 'correlation', lag_max = 2, 
                                         intensity = vertex_attr(g)$intensity)
  
  data_moran <- NodeLocalCorrelation(intnet_all, dep_type = 'moran_i', intensity = vertex_attr(g)$intensity)
  moran_i <- data_moran$correlation
  intnet_all <- data_moran$intnet_all
  
  data_geary <- NodeLocalCorrelation(intnet_all, dep_type = 'geary_g', intensity = vertex_attr(g)$intensity)
  geary <- data_geary$correlation
  intnet_all <- data_geary$intnet_all
  
} else{
  vertex_attr(g, 'intensity_in', V(g)['V1']) 
  vertex_attr(g, 'intensity_out', V(g)['V1']) 
  
  pdf("Plots/area_with_grid_Dir_plot.pdf")
  plot(intnet_all, enable_grid = TRUE, vertex_intensity = 'intensity_out')
  dev.off()
  
  pdf("Plots/area_with_grid_Dir_gplot.pdf")
  gplot(intnet_all)
  dev.off()
  
  pdf("Plots/area_with_grid_Dir_moran_gplot.pdf")
  gplot(intnet_all, intensity = vertex_attr(intnet_all$graph)$intensity_in, heatmap = 'moran_i')
  dev.off()
  
  pdf("Plots/area_with_grid_Dir_g_gplot.pdf")
  gplot(intnet_all, intensity = vertex_attr(intnet_all$graph)$intensity_in, heatmap = 'geary_g')
  dev.off()
  
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
  
  correlations <- NodeGeneralCorrelation(intnet_all, dep_type = 'correlation', lag_max = 2, 
                                         intensity = vertex_attr(g)$intensity_in)
  
  data_moran <- NodeLocalCorrelation(intnet_all, dep_type = 'moran_i', intensity = vertex_attr(g)$intensity_in)
  moran_i <- data_moran$correlation
  intnet_all <- data_moran$intnet
  
  data_geary <- NodeLocalCorrelation(intnet_all, dep_type = 'geary_g', intensity = vertex_attr(g)$intensity_in)
  geary <- data_geary$correlation
  intnet_all <- data_geary$intnet
}

for(edge_id in E(g)){
  if(E(g)[edge_id]$intensity>0) print(E(g)[edge_id]$intensity)
}

#-----------------------------INTENSITYNET CLASS: FUNCTION TESTING---------------------------------




#--------------------------------------------PLOTS------------------------------------------------
pdf("S3/Plots/area_with_grid_Dir.pdf")
plot(intnet_all, enable_grid = TRUE, vertex_intensity = 'intensity_in')
dev.off()

pdf("Plots/area_with_grid.pdf")
plot(intnet_all, enable_grid = FALSE, axis=TRUE)
dev.off()



plot(intnet_all, node_label = 'intensity', edge_label='none', vertex.color='red')

gplot(intnet_all)
gplot(intnet_all, heatmap = 'locmoran')
gplot(intnet_all, heatmap = 'locg')



