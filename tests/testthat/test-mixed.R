data("mix_intnet_chicago")



test_that("Calculate edge and nodemeans intensities from a mixed network", {
  library(spatstat)
  data(chicago)
  
  chicago_df <- as.data.frame(chicago[["data"]])
  edges <- cbind(chicago[["domain"]][["from"]], chicago[["domain"]][["to"]])
  chicago_net <- igraph::graph_from_edgelist(edges)
  chicago_adj_mtx <- as.matrix(igraph::as_adjacency_matrix(chicago_net))
  chicago_node_coords <- data.frame(xcoord = chicago[["domain"]][["vertices"]][["x"]], 
                                    ycoord = chicago[["domain"]][["vertices"]][["y"]])
  
  # Generate undirected intensitynet object
  intnet_chicago <- intensitynet(chicago_adj_mtx, 
                                 node_coords = chicago_node_coords, 
                                 event_data = chicago_df,
                                 graph_type='mixed')
  
  
  intnet_chicago <- RelateEventsToNetwork(intnet_chicago)
  
  
  expect_s3_class(intnet_chicago, c("intensitynet", "intensitynetMix"), exact = TRUE)
})


test_that("Network general covariance", {
  intnet <- mix_intnet_chicago
  gen_cov <- NodeGeneralCorrelation(intnet, 
                                    dep_type = 'covariance', 
                                    lag_max = 2, 
                                    intensity = igraph::vertex_attr(intnet$graph)$intensity_in)
  
  expect_gte(length(gen_cov), 1)
})


test_that('Node local moran i', {
  intnet <- mix_intnet_chicago
  
  data_moran <- NodeLocalCorrelation(intnet, 
                                     dep_type = 'moran', 
                                     intensity = igraph::vertex_attr(intnet$graph)$intensity_in)
  moran_i <- data_moran$correlation
  intnet <- data_moran$intnet
  
  expect_gte(length(moran_i), 1)
})


test_that('Node local geary c', {
  intnet <- mix_intnet_chicago
  
  data_geary <- NodeLocalCorrelation(intnet, 
                                     dep_type = 'geary', 
                                     intensity = igraph::vertex_attr(intnet$graph)$intensity_in)
  geary <- data_geary$correlation
  intnet <- data_geary$intnet
  
  expect_gte(length(data_geary), 1)
})

test_that('Node local getis g', {
  intnet <- mix_intnet_chicago
  g <- intnet$graph
  
  data_getis <- NodeLocalCorrelation(intnet, dep_type = 'getis', intensity = igraph::vertex_attr(g)$intensity_in)
  getis <- data_getis$correlation
  intnet <- data_getis$intnet
  
  expect_gte(length(getis), 1)
})


test_that('Path Intensity', {
  intnet <- mix_intnet_chicago
  
  short_dist <- ShortestPathIntensity(intnet, node_id1 = 'V1', node_id2 = 'V300')
  int_path <- short_dist$intensity
  
  expect_gte(int_path, 0)
})