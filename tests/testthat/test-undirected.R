load("../Castellon.RData")
load("../nodes.RData")
load("../crimes.RData")

crim <- crimes[11:111,] # From crimes, take 11 to 111 (both included)

# Generate undirected intensitynet object
intnet <- intensitynet(Castellon, nodes, crim) 


test_that("Calculate edge and nodemeans intensities from the undirected network", {
  intnet_all <- CalculateEventIntensities(intnet)
  
  expect_s3_class(intnet_all, c("intensitynet", "intensitynetUnd"), exact = TRUE)
})

test_that("Get intensities", {
  vertex_intensity <- vertex_attr(intnet_all$graph, 'intensity') 
  
  expect_gt(vertex_intensity, 0)
})