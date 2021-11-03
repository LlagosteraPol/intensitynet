load("../Castellon.RData")
load("../nodes.RData")
load("../crimes.RData")

load("../undirected_intensitynet.RData")

test_that("Calculate edge and nodemeans intensities from an undirected network", {
  crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)
  intnet <- intensitynet(Castellon, nodes, crim) # Generate undirected intensitynet object
  intnet_all <- CalculateEventIntensities(intnet)
  
  expect_s3_class(intnet_all, c("intensitynet", "intensitynetUnd"), exact = TRUE)
})

test_that("Network general covariance", {
  gen_cov <- NodeGeneralCorrelation(intnet_und, 
                                    dep_type = 'covariance', 
                                    lag_max = 2, 
                                    intensity = vertex_attr(intnet_und$graph)$intensity)

  expect_gte(length(gen_cov), 1)
})

test_that('Node local moran i', {
  intnet <- intnet_und
  g <- intnet$graph
  
  data_moran <- NodeLocalCorrelation(intnet, dep_type = 'moran_i', intensity = vertex_attr(g)$intensity)
  moran_i <- data_moran$correlation
  intnet <- data_moran$intnet
  
  expect_gte(length(moran_i), 1)
})

test_that('Node local geary g', {
  intnet <- intnet_und
  g <- intnet$graph
  
  data_geary <- NodeLocalCorrelation(intnet, dep_type = 'geary_g', intensity = vertex_attr(g)$intensity)
  geary <- data_geary$correlation
  intnet <- data_geary$intnet
  
  expect_gte(length(data_geary), 1)
})

test_that('Path Intensity', {
  short_dist_obj <- list(graph=intnet_und$graph, node_id1 = 'V601', node_id2 = 'V701', distances_mtx = intnet_und$distances)
  class(short_dist_obj) <- "netTools"
  short_dist <- ShortestDistance(short_dist_obj)
  
  int_path <- PathIntensity(intnet_und, short_dist$path)
  
  expect_gte(int_path, 0)
})