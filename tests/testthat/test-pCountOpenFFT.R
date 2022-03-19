test_that("pCountOpenFFT works, 1 site, 3 times", {
  # setup parallel processing for test
  cl = parallel::makeCluster(1)
  parallel::setDefaultCluster(cl=cl)
  doParallel::registerDoParallel(cl)
  
  nit = matrix(c(1,1,1), nrow=1)
  expect_error(pCountOpenFFT(nit = nit, K = 2), NA)
  
  # stop parallel processing
  parallel::stopCluster(cl)
})

test_that("pCountOpenFFT works, 2 sites, 4 times", {
  # setup parallel processing for test
  cl = parallel::makeCluster(1)
  parallel::setDefaultCluster(cl=cl)
  doParallel::registerDoParallel(cl)
  
  nit = matrix(c(1,1,1,1,1,1,1,1), nrow = 2)
  #expect_error(pCountOpenFFT(nit = nit, K = 12, outfile = "test_file.csv"), NA)
  expect_error(pCountOpenFFT(nit = nit, K = 2), NA)
  
  # stop parallel processing
  parallel::stopCluster(cl)
})