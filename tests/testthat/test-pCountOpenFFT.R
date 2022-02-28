test_that("pCountOpenFFT works, 1 site, 3 times", {
  nit = matrix(anmu[1,1:3], nrow=1)
  expect_error(pCountOpenFFT(nit = nit, K = 50), NA)
})

test_that("pCountOpenFFT works, 2 sites, 4 times", {
  nit = matrix(c(10,8,11,13,9,7,10,12), nrow=2)
  expect_error(pCountOpenFFT(nit = nit, K = 20, outfile = "test_file.csv"), NA)
})