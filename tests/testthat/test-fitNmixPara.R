test_that("fitNmixPara works, 1 site, 3 times", {
  nit = matrix(anmu[1,1:3], nrow=1)
  cl = parallel::makeCluster(parallel::detectCores()-1)
  parallel::setDefaultCluster(cl=cl)
  expect_error(fitNmixPara(cluster=cl, nit = nit, K = 50), NA)
  parallel::stopCluster(cl)
})

test_that("fitNmixPara works, 2 sites, 4 times", {
  nit = matrix(c(10,8,11,13,9,7,10,12), nrow=2)
  cl = parallel::makeCluster(parallel::detectCores()-1)
  parallel::setDefaultCluster(cl=cl)
  expect_error(fitNmixPara(cluster=cl, nit = nit, K = 20, outfile = "test_file.csv"), NA)
  expect_error(fitNmixPara(cluster=cl, nit = nit, K = 20, l_s_c=list(c(0,1))), NA)
  expect_error(fitNmixPara(cluster=cl, nit = nit, K = 20, l_s_c=list(c(0,1)), g_t_c=list(c(0,1,0,1))), NA)
  expect_error(fitNmixPara(cluster=cl, nit = nit, K = 20, p_s_c=list(c(0,1)), o_t_c=list(c(0,1,0,1))), NA)
  parallel::stopCluster(cl)
})