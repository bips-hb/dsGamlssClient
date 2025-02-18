# Test-smk-getPooledMean ################################################

#
# Set up
#

connect.studies.dataset.gamlss(list("e3_bw", "e3_gac_None"))

#
# Tests
#

test_that("output_getPooledMean_D$e3_bw", {
  res <- getPooledMean(ds.test_env$connections, "D$e3_bw")
  expect_length(class(res), 1)
  expect_true(is.numeric(res))
  expect_equal(res, mean(ds.test_env$local.values$e3_bw), tolerance=1e-07)
})

test_that("output_getPooledMean_D_red$e3_bw", {
  res <- getPooledMean(ds.test_env$connections, "D_red$e3_bw")
  expect_length(class(res), 1)
  expect_true(is.numeric(res))
  expect_equal(res, mean(ds.test_env$gamlss_red$e3_bw), tolerance=1e-07)
})

test_that("output_getPooledMean_D$e3_gac_None", {
  res <- getPooledMean(ds.test_env$connections, "D$e3_gac_None")
  expect_length(class(res), 1)
  expect_true(is.numeric(res))
  expect_equal(res, mean(ds.test_env$local.values$e3_gac_None), tolerance=1e-07)
})

test_that("output_getPooledMean_D_red$e3_gac_None", {
  res <- getPooledMean(ds.test_env$connections, "D_red$e3_gac_None")
  expect_length(class(res), 1)
  expect_true(is.numeric(res))
  expect_equal(res, mean(ds.test_env$gamlss_red$e3_gac_None), tolerance=1e-07)
})

#
# Done
#

disconnect.studies.dataset.gamlss()