data("und_intnet_chicago")

test_that('ApplyWindow function',{
  intnet <- und_intnet_chicago
  sub_intnet_chicago <- ApplyWindow(intnet, 
                                    x_coords = c(300, 900), 
                                    y_coords = c(500, 1000))
  
  expect_s3_class(sub_intnet_chicago, c("intensitynetUnd", "intensitynet"), exact = TRUE)
})


