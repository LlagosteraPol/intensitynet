load("../Und_Intnet_Chicago_Data.RData")


test_that('NodeLocalCorrelation function, bad argument handling', {
  intnet <- und_intnet_chicago
  g <- intnet$graph
  
  expect_error(NodeLocalCorrelation(intnet, dep_type = 'bad_arg', intensity = igraph::vertex_attr(g)$intensity),
              "Error: wrong 'dep_type' argument. Allowed arguments are; 'moran', 'getis', 'geary'.",
              fixed=TRUE)
})


test_that('ApplyWindow function',{
  intnet <- und_intnet_chicago
  sub_intnet_chicago <- ApplyWindow(intnet, 
                                    x_coords = c(300, 900), 
                                    y_coords = c(500, 1000))
  
  expect_s3_class(sub_intnet_chicago, c("intensitynet", "intensitynetUnd"), exact = TRUE)
})


