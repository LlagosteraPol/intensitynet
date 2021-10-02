
# to throw error: stop(sprintf("The object must be of class 'intensitynet', current class is %s", class(tt)[1]))


rm(list = ls())

source("S3/main.R")


# Adjacency matrix (undirected): Segmenting locations of the traffic network treated as the vertex set of the network.
load("../Data/Castellon.RData")

# Node coordinates: Georeferenced coordinates from 'castellon' nodes
load("../Data/nodes.RData")

# Event (crime coordinates)
load("../Data/crimes.RData")


#subset of events
crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)


init_graph <- function(obj){
  weighted_mtx = obj$adjacency_mtx * obj$distances
  g <- graph_from_adjacency_matrix(weighted_mtx, mode = obj$graph_type, weighted=TRUE)
  
  net_coords <- list(graph = g, node_coords = obj$node_coords)
  class(net_coords) <- "netTools"
  g <- setNetCoords(net_coords)
  
  g # return
}

setNetCoords <- function(obj){
  UseMethod("setNetCoords")
}

#TODO: Declare non-visible at namespace
setNetCoords.netTools = function(obj){
  x_coord_node <- obj$node_coords[1]
  y_coord_node <- obj$node_coords[2]
  
  # TODO: change x and y to long (longitude) and lat (latitude), beware of the coordinate system type (e.g WGS84)
  g <- obj$g %>% set_vertex_attr(name = "xcoord", value = x_coord_node) %>% 
    set_vertex_attr(name = "ycoord", value = y_coord_node)
  g
}

# ----------------------TESTING 'shortestDistance'-------------------------

shortestDistance = function(g, node_id1, node_id2, distances_mtx){
  
  weighted_path <- unlist(get.shortest.paths(g, node_id1, node_id2)$vpath)
  if(!is.null(distances_mtx)){
    weight_sum <- sum(E(g, path = unlist(weighted_path))$weight)
  }
  else{
    weight_sum <- length(weighted_path)
  }
  list(path = weighted_path, weight = weight_sum)  
}

nn <- as.matrix(nodes)
node_coords_obj <- list(node_coords=nn)
class(node_coords_obj) <- "netTools"
dist_mtx <- calculateDistancesMtx(node_coords_obj)

net_setup <- list(adjacency_mtx = as.matrix(Castellon), 
                  node_coords = nn, 
                  events = as.matrix(crimes), 
                  distances = dist_mtx, 
                  graph_type = 'undirected')

g <- init_graph(net_setup)
intnet <- intensitynet(Castellon, nodes, crim)
shortest <- shortestDistance(intnet$graph, 'V1', 'V20', intnet$distances)

# ---------------Setting vertice attributes-----------------------
vertex_attr(g, 'intensity', index='V1') 
g <- g %>% set_vertex_attr("intensity", index='V1', value = 333)
vertex_attr(g, 'intensity','V2') 

if(!is.null(vertex_attr(g, "intensity", index='V2'))){
  if(!is.na(vertex_attr(g, "intensity", index='V2')))  vertex_attr(g, "intensity", 'V1')
}

g <- g %>% set_edge_attr(name = "test", index = E(g)[1,2,3], value = c(10,20,30))
class(edge_attr(g)[1])
edge_attr(g)[2]
edge_attr_names(g)
edge_attr(g, "test", E(gt)[1,2,3])


library(sna)
g <- intnet_all$graph
nacf(g, vertex_attr(g, 'intensity'))


# ----------------------------------Window between two points--------------------
node1 <- c(2, 5)
node2 <- c(7, 1)
z <- 0.5

dx <- node1[1] - node2[1]
dy <- node1[2] - node2[2]

if(dx<dy){
  zx <- z/(sqrt(1+(dx/dy)^2))
  zy <- -(dx/dy)*zx
}else{
  zy <- z/(sqrt(1+(dy/dx)^2))
  zx <- -(dy/dx)*zy
}

p1 <- c(node1[1] - zx, node1[2] - zy)
p2 <- c(node1[1] + zx, node1[2] + zy)
p3 <- c(node2[1] - zx, node2[2] - zy)
p4 <- c(node2[1] + zx, node2[2] + zy)

plot(1, type = "n",                         # Remove all elements of plot
     xlab = "", ylab = "",
     xlim = c(749000, 750000), ylim = c(4430500, 4431000))
lines(c(node1[1], node2[1]), c(node1[2], node2[2]), type = "l", lty = 1)
points(c(p1[1], p2[1], p3[1], p4[1]), c(p1[2], p2[2], p3[2], p4[2]))

lines(c(p1[1], p3[1]), c(p1[2], p3[2]), type = "l", lty = 1)
lines(c(p2[1], p4[1]), c(p2[2], p4[2]), type = "l", lty = 1)

lines(c(p1[1], p2[1]), c(p1[2], p2[2]), type = "l", lty = 1)
lines(c(p3[1], p4[1]), c(p3[2], p4[2]), type = "l", lty = 1)

library(contoureR)
win_points <- list(p1, p2, p3, p4)
anticlockwise <- orderPoints(x=c(p1[1], p2[1], p3[1], p4[1]), 
                             y=c(p1[2], p2[2], p3[2], p4[2]), clockwise = FALSE)

win_points <- win_points[order(anticlockwise)]

win <- tryCatch(
  {
    owin(poly=list(x=c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1]),
                   y=c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2])))
  },
  error=function(cond) {
    owin(poly=list(x=rev(c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1])),
                   y=rev(c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2]))))
  }
)


plot(win)


df <- data.frame(x=c(win_points[[1]][1], win_points[[2]][1], win_points[[3]][1], win_points[[4]][1]),
                 y=c(win_points[[1]][2], win_points[[2]][2], win_points[[3]][2], win_points[[4]][2]))
clockwise(df)
clockwise <- function(x) {
  
  x.coords <- c(x[[1]], x[[1]][1])
  y.coords <- c(x[[2]], x[[2]][1])
  
  double.area <- sum(sapply(2:length(x.coords), function(i) {
    (x.coords[i] - x.coords[i-1])*(y.coords[i] + y.coords[i-1])
  }))
  
  double.area > 0
} 

dist_obj <- list(p1= node1, p2= node2, ep=c(3,3))
class(dist_obj) <- 'netTools'
perpenidularDistance(dist_obj)

