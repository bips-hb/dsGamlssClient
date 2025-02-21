# Test-smk-getPooledVar ################################################

#
# Set up
#

connect.studies.dataset.gamlss(list("e3_bw", "e3_gac_None"))

#
# Tests
#

test_that("output_getPooledVar_D$e3_bw", {
  res <- getPooledVar(ds.test_env$connections, "D$e3_bw")
  weighted.var1 <- nrow(ds.test_env$local.values.1) * var(ds.test_env$local.values.1$e3_bw)
  weighted.var2 <- nrow(ds.test_env$local.values.2) * var(ds.test_env$local.values.2$e3_bw)
  weighted.var3 <- nrow(ds.test_env$local.values.3) * var(ds.test_env$local.values.3$e3_bw)
  pooled.var <- (weighted.var1 + weighted.var2 + weighted.var3) / nrow(ds.test_env$local.values)
  expect_length(class(res), 1)
  expect_true(is.numeric(res))
  expect_equal(res, pooled.var, tolerance = 1e-07)
})

test_that("output_getPooledVar_D_red$e3_bw", {
  load(testthat::test_path("data_files", "GAMLSS", "gamlss_red.rda"))
  res <- getPooledVar(ds.test_env$connections, "D_red$e3_bw")
  weighted.var1 <- nrow(gamlss_red) * var(gamlss_red$e3_bw)
  pooled.var <- weighted.var1 / nrow(gamlss_red)
  expect_length(class(res), 1)
  expect_true(is.numeric(res))
  expect_equal(res, pooled.var, tolerance = 1e-07)
})

test_that("output_getPooledVar_D$e3_gac_None", {
  res <- getPooledVar(ds.test_env$connections, "D$e3_gac_None")
  weighted.var1 <- nrow(ds.test_env$local.values.1) * var(ds.test_env$local.values.1$e3_gac_None)
  weighted.var2 <- nrow(ds.test_env$local.values.2) * var(ds.test_env$local.values.2$e3_gac_None)
  weighted.var3 <- nrow(ds.test_env$local.values.3) * var(ds.test_env$local.values.3$e3_gac_None)
  pooled.var <- (weighted.var1 + weighted.var2 + weighted.var3) / nrow(ds.test_env$local.values)
  expect_length(class(res), 1)
  expect_true(is.numeric(res))
  expect_equal(res, pooled.var, tolerance = 1e-07)
})

test_that("output_getPooledVar_D_red$e3_gac_None", {
  load(testthat::test_path("data_files", "GAMLSS", "gamlss_red.rda"))
  res <- getPooledVar(ds.test_env$connections, "D_red$e3_gac_None")
  weighted.var1 <- nrow(gamlss_red) * var(gamlss_red$e3_gac_None)
  pooled.var <- weighted.var1 / nrow(gamlss_red)
  expect_length(class(res), 1)
  expect_true(is.numeric(res))
  expect_equal(res, pooled.var, tolerance = 1e-07)
})

#
# Done
#

disconnect.studies.dataset.gamlss()
