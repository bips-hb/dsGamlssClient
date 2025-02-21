# Test-smk-gamlssChecks ################################################

#
# Set up
#

connect.studies.dataset.gamlss(list(
  "e3_bw", "e3_gac_None", "hs_zbmi_who", "hs_child_age_None",
  "h_mbmi_None", "hs_correct_raven", "hs_wgtgain_None"
))

#
# Tests
#

test_that("output_gamlssChecks_withprefix", {
  res <- gamlssChecks(
    formula = stats::as.formula(D$e3_bw ~ pb(D$e3_gac_None)),
    sigma.formula = stats::as.formula(~1),
    nu.formula = stats::as.formula(~1),
    tau.formula = stats::as.formula(~1),
    data = "D",
    datasources = ds.test_env$connections
  )
  expect_true(is.null(res))
})

test_that("output_gamlssChecks_withprefix_dataNULL", {
  res <- gamlssChecks(
    formula = stats::as.formula(D$e3_bw ~ pb(D$e3_gac_None)),
    sigma.formula = stats::as.formula(~1),
    nu.formula = stats::as.formula(~1),
    tau.formula = stats::as.formula(~1),
    datasources = ds.test_env$connections
  )
  expect_true(is.null(res))
})

test_that("output_gamlssChecks_noprefix", {
  res <- gamlssChecks(
    formula = stats::as.formula(e3_bw ~ pb(e3_gac_None)),
    sigma.formula = stats::as.formula(~1),
    nu.formula = stats::as.formula(~1),
    tau.formula = stats::as.formula(~1),
    data = "D",
    datasources = ds.test_env$connections
  )
  expect_true(is.null(res))
})

test_that("output_gamlssChecks_noprefix_dataNULL", {
  expect_error(
    gamlssChecks(
      formula = stats::as.formula(e3_bw ~ pb(e3_gac_None)),
      sigma.formula = stats::as.formula(~1),
      nu.formula = stats::as.formula(~1),
      tau.formula = stats::as.formula(~1),
      data = NULL,
      datasources = ds.test_env$connections
    ),
    "No data.frame for the column e3_bw given. Specify it explicitly as dataname$e3_bw or provide a valid data argument.",
    fixed = TRUE
  )
})

test_that("output_gamlssChecks_NA", {
  expect_error(
    gamlssChecks(
      formula = stats::as.formula(e3_bw ~ pb(na_var)),
      sigma.formula = stats::as.formula(~1),
      nu.formula = stats::as.formula(~1),
      tau.formula = stats::as.formula(~1),
      data = "D_na",
      datasources = ds.test_env$connections
    ),
    "The variable na_var in server2 is missing at complete (all values are 'NA').",
    fixed = TRUE
  )
})

#
# Done
#

disconnect.studies.dataset.gamlss()
