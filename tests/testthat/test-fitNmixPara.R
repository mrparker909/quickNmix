
test_that("fitNmixPara works, 2 sites, 4 times", {
  # setup parallel processing for test
  cl = parallel::makeCluster(2)
  parallel::setDefaultCluster(cl=cl)
  doParallel::registerDoParallel(cl)
  
  nit = matrix(c(10,8,11,13,9,7,10,12), nrow=2)
  expect_error(fitNmixPara(cluster=cl, nit = nit, K = 20, NA))
  #expect_error(fitNmixPara(cluster=cl, nit = nit, K = 20, outfile = "test_file.csv"), NA)
  #expect_error(fitNmixPara(cluster=cl, nit = nit, K = 20, l_s_c=list(c(0,1))), NA)
  #expect_error(fitNmixPara(cluster=cl, nit = nit, K = 20, l_s_c=list(c(0,1)), g_t_c=list(c(0,1,0,1))), NA)
  #expect_error(fitNmixPara(cluster=cl, nit = nit, K = 20, p_s_c=list(c(0,1)), o_t_c=list(c(0,1,0,1))), NA)
  
  # stop parallel processing
  parallel::stopCluster(cl)
})
