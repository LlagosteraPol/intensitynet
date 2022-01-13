load("../Und_Intnet_Chicago_Data.RData")


test_that("Calculate edge and nodemeans intensities from an undirected network", {
  library(spatstat)
  data(chicago)
  
  chicago_df <- as.data.frame(chicago[["data"]])
  edges <- cbind(chicago[["domain"]][["from"]], chicago[["domain"]][["to"]])
  chicago_net <- igraph::graph_from_edgelist(edges)
  chicago_adj_mtx <- as.matrix(igraph::as_adjacency_matrix(chicago_net))
  chicago_node_coords <- data.frame(xcoord = chicago[["domain"]][["vertices"]][["x"]], 
                                    ycoord = chicago[["domain"]][["vertices"]][["y"]])
  chicago_assault <- chicago_df[chicago_df$marks == 'assault',]
  assault_coordinates <- data.frame(xcoord = chicago_assault[,1],
                                    ycoord = chicago_assault[,2])
  
  # Generate undirected intensitynet object
  intnet_chicago <- intensitynet(chicago_adj_mtx, 
                                 node_coords = chicago_node_coords, 
                                 event_coords = assault_coordinates)
  
  
  intnet_chicago <- CalculateEventIntensities(intnet_chicago)
  
  
  expect_s3_class(intnet_chicago, c("intensitynet", "intensitynetUnd"), exact = TRUE)
})


test_that("Network general covariance", {
  intnet <- und_intnet_chicago
  gen_cov <- NodeGeneralCorrelation(intnet, 
                                    dep_type = 'covariance', 
                                    lag_max = 2, 
                                    intensity = igraph::vertex_attr(intnet$graph)$intensity)

  expect_gte(length(gen_cov), 1)
})


test_that('Node local moran i', {
  intnet <- und_intnet_chicago
  g <- intnet$graph
  
  data_moran <- NodeLocalCorrelation(intnet, dep_type = 'moran', intensity = igraph::vertex_attr(g)$intensity)
  moran_i <- data_moran$correlation
  intnet <- data_moran$intnet
  
  expect_gte(length(moran_i), 1)
})


test_that('Node local geary c', {
  intnet <- und_intnet_chicago
  g <- intnet$graph
  
  data_geary <- NodeLocalCorrelation(intnet, dep_type = 'geary', intensity = igraph::vertex_attr(g)$intensity)
  geary <- data_geary$correlation
  intnet <- data_geary$intnet
  
  expect_gte(length(geary), 1)
})


test_that('Node local getis g', {
  intnet <- und_intnet_chicago
  g <- intnet$graph
  
  data_getis <- NodeLocalCorrelation(intnet, dep_type = 'getis', intensity = igraph::vertex_attr(g)$intensity)
  getis <- data_getis$correlation
  intnet <- data_getis$intnet
  
  expect_gte(length(getis), 1)
})


test_that('Path Intensity', {
  intnet <- und_intnet_chicago
  
  short_dist <- ShortestNodeDistance(intnet, node_id1 = 'V1', node_id2 = 'V150')
  int_path <- PathIntensity(intnet, short_dist$path)
  
  expect_gte(int_path, 0)
})