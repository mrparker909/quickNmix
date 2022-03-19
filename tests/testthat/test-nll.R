test_that("nll works, 1 site", {
  # setup parallel processing for test
  cl = parallel::makeCluster(2)
  parallel::setDefaultCluster(cl=cl)
  doParallel::registerDoParallel(cl)
  
  nit = matrix(anmu[1,], nrow=1)
  expect_error(nll(par = c(1,1,0,0), nit = nit, K = 10))
  expect_error(nll(par = c(1,1,0,0), nit = nit, K = 100), NA)
  
  # stop parallel processing
  parallel::stopCluster(cl)
})

test_that("nll works, 2 site", {
  # setup parallel processing for test
  cl = parallel::makeCluster(2)
  parallel::setDefaultCluster(cl=cl)
  doParallel::registerDoParallel(cl)
  
  nit = anmu[1:2,]
  expect_error(nll(par = c(1,1,0,0), nit = nit, K = 100))
  expect_error(nll(par = c(1,1,0,0), nit = nit, K = 140), NA)
  
  # stop parallel processing
  parallel::stopCluster(cl)
})