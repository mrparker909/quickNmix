test_that("fitNmix works, 1 site, 3 times", {
  # setup parallel processing for test
  cl = parallel::makeCluster(1)
  parallel::setDefaultCluster(cl=cl)
  doParallel::registerDoParallel(cl)
  
  nit = matrix(c(1,1,1), nrow = 1)
  
  expect_error(fitNmix(nit = nit, K = 5), NA)
  
  # stop parallel processing
  parallel::stopCluster(cl)
})

test_that("fitNmix works, 2 sites, 3 times", {
  # setup parallel processing for test
  cl = parallel::makeCluster(1)
  parallel::setDefaultCluster(cl=cl)
  doParallel::registerDoParallel(cl)
  
  nit = matrix(c(1,1,1,1,1,1), nrow = 2)
  expect_error(fitNmix(nit = nit, K = 5), NA)
  #expect_error(fitNmix(nit = nit, K = 20, outfile = "test_file.csv"), NA)
  #expect_error(fitNmix(nit = nit, K = 20, l_s_c=list(c(0,1))), NA)
  #expect_error(fitNmix(nit = nit, K = 20, l_s_c=list(c(0,1)), g_t_c=list(c(0,1,0,1))), NA)
  #expect_error(fitNmix(nit = nit, K = 20, p_s_c=list(c(0,1)), o_t_c=list(c(0,1,0,1))), NA)
  
  # stop parallel processing
  parallel::stopCluster(cl)
})