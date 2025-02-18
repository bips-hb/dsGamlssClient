# Test-smk-extract ################################################

#
# Set up
#

#
# Tests
#

test_that("output_single_string", {
  res <- extract("D$e3_bw")
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 2)
  expect_equal(names(res), c("holders", "elements"))
  expect_equal(res$holders, "D")
  expect_equal(res$elements, "e3_bw")
})

test_that("output_single_string_noholder", {
  res <- extract("e3_bw")
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 2)
  expect_equal(names(res), c("holders", "elements"))
  expect_true(is.na(res$holders))
  expect_equal(res$elements, "e3_bw")
})

test_that("output_extract_vector", {
  res <- extract(c("D$e3_bw", "E$e3_gac_None"))
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 2)
  expect_equal(names(res), c("holders", "elements"))
  expect_equal(res$holders, c("D", "E"))
  expect_equal(res$elements, c("e3_bw", "e3_gac_None"))
})

test_that("output_extract_vector_noholder", {
  res <- extract(c("D$e3_bw", "e3_gac_None"))
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 2)
  expect_equal(names(res), c("holders", "elements"))
  expect_equal(res$holders, c("D", NA))
  expect_equal(res$elements, c("e3_bw", "e3_gac_None"))
})

test_that("output_extract_list", {
  res <- extract(list("D$e3_bw", "E$e3_gac_None"))
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 2)
  expect_equal(names(res), c("holders", "elements"))
  expect_equal(res$holders, c("D", "E"))
  expect_equal(res$elements, c("e3_bw", "e3_gac_None"))
})

test_that("output_extract_list_noholder", {
  res <- extract(list("D$e3_bw", "e3_gac_None"))
  expect_length(class(res), 1)
  expect_true(all("list" %in% class(res)))
  expect_length(res, 2)
  expect_equal(names(res), c("holders", "elements"))
  expect_equal(res$holders, c("D", NA))
  expect_equal(res$elements, c("e3_bw", "e3_gac_None"))
})

#
# Done
#