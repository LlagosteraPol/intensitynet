load("../Castellon.RData")
load("../nodes.RData")
load("../crimes.RData")

load("../directed_intensitynet.RData")

test_that("Calculate edge and nodemeans intensities from a directed network", {
  crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)
  intnet <- intensitynet(dir_castellon, nodes, crim, graph_type='directed') # Generate directed intensitynet object
  intnet_all <- CalculateEventIntensities(intnet)
  
  expect_s3_class(intnet_all, c("intensitynet", "intensitynetDir"), exact = TRUE)
})

test_that("Network general covariance", {
  gen_cov <- NodeGeneralCorrelation(intnet_dir, 
                                    dep_type = 'covariance', 
                                    lag_max = 2, 
                                    intensity = vertex_attr(intnet_dir$graph)$intensity_in)
  
  expect_gte(length(gen_cov), 1)
})

test_that('Node local moran i', {
  intnet <- intnet_dir
  g <- intnet$graph
  
  data_moran <- NodeLocalCorrelation(intnet, dep_type = 'moran_i', intensity = vertex_attr(g)$intensity_in)
  moran_i <- data_moran$correlation
  intnet <- data_moran$intnet
  
  expect_gte(length(moran_i), 1)
})

test_that('Node local geary g', {
  intnet <- intnet_dir
  g <- intnet$graph
  
  data_geary <- NodeLocalCorrelation(intnet, dep_type = 'geary_g', intensity = vertex_attr(g)$intensity_in)
  geary <- data_geary$correlation
  intnet <- data_geary$intnet
  
  expect_gte(length(data_geary), 1)
})

test_that('Path Intensity', {
  short_dist_obj <- list(graph=intnet_dir$graph, node_id1 = 'V601', node_id2 = 'V701', distances_mtx = intnet_dir$distances)
  class(short_dist_obj) <- "netTools"
  short_dist <- ShortestDistance(short_dist_obj)
  
  int_path <- PathIntensity(intnet_dir, short_dist$path)
  
  expect_gte(int_path, 0)
})