test_that("nll works, 1 site", {
  nit = matrix(anmu[1,], nrow=1)
  expect_error(nll(par = c(1,1,0,0), nit = nit, K = 10))
  expect_error(nll(par = c(1,1,0,0), nit = nit, K = 100), NA)
})

test_that("nll works, 2 site", {
  nit = anmu[1:2,]
  expect_error(nll(par = c(1,1,0,0), nit = nit, K = 100))
  expect_error(nll(par = c(1,1,0,0), nit = nit, K = 140), NA)
})