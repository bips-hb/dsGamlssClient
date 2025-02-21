# Test-smk-isDefined ################################################

#
# Set up
#

connect.studies.dataset.gamlss(list("e3_bw", "e3_gac_None"))

#
# Tests
#

test_that("output_isDefined_dataframeD", {
  res <- isDefined(ds.test_env$connections, "D")
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 3)
  expect_true(res$server1)
  expect_true(res$server2)
  expect_true(res$server3)
})

test_that("output_isDefined_dataframeE", {
  expect_error(isDefined(ds.test_env$connections, "E"), "The input object E is not defined in server1, server2, server3!", fixed = TRUE)
})

test_that("output_isDefined_dataframeD_columne3_bw", {
  res <- isDefined(ds.test_env$connections, "D$e3_bw")
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 3)
  expect_true(res$server1)
  expect_true(res$server2)
  expect_true(res$server3)
})

test_that("output_isDefined_dataframeD_columnA", {
  expect_error(isDefined(ds.test_env$connections, "D$A"), "The input object D$A is not defined in server1, server2, server3!", fixed = TRUE)
})

test_that("output_isDefined_dataframeE_columnA", {
  expect_error(isDefined(ds.test_env$connections, "E$A"), "There are some DataSHIELD errors:", fixed = TRUE)
})

# error.message = FALSE

test_that("output_isDefined_dataframeD_errormessageFALSE", {
  res <- isDefined(ds.test_env$connections, "D", error.message = FALSE)
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 3)
  expect_true(res$server1)
  expect_true(res$server2)
  expect_true(res$server3)
})

test_that("output_isDefined_dataframeE_errormessageFALSE", {
  res <- isDefined(ds.test_env$connections, "E", error.message = FALSE)
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 3)
  expect_false(res$server1)
  expect_false(res$server2)
  expect_false(res$server3)
})

test_that("output_isDefined_dataframeD_columne3_bw", {
  res <- isDefined(ds.test_env$connections, "D$e3_bw", error.message = FALSE)
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 3)
  expect_true(res$server1)
  expect_true(res$server2)
  expect_true(res$server3)
})

test_that("output_isDefined_dataframeD_columnA_errormessageFALSE", {
  res <- isDefined(ds.test_env$connections, "D$A", error.message = FALSE)
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 3)
  expect_false(res$server1)
  expect_false(res$server2)
  expect_false(res$server3)
})

test_that("output_isDefined_dataframeE_columnA_errormessageFALSE", {
  expect_error(isDefined(ds.test_env$connections, "E$A", error.message = FALSE), "There are some DataSHIELD errors:", fixed = TRUE)
})

#
# Done
#

disconnect.studies.dataset.gamlss()
