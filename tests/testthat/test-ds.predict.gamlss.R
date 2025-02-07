# Test-arg-ds.predict.gamlss ################################################

#
# Set up
#

connect.studies.dataset.gamlss(list("e3_bw", "e3_gac_None", "hs_zbmi_who", "hs_child_age_None",
                                    "h_mbmi_None", "hs_correct_raven", "hs_wgtgain_None"))

#
# Tests
#

test_that("input_errors", {
  model_e3_bw.DS <- ds.gamlss(formula = D$e3_bw ~ D$e3_gac_None,
                              data = 'D', family = 'NO()')
  newdata <- data.frame(xvar = seq(28, 44, by=0.1))
  names(newdata) <- "e3_gac_None"
  expect_error(ds.predict.gamlss(object="notagamlssobject", newdata=NULL), "Please provide a valid ds.gamlss.object to derive the predictions.", fixed=TRUE)
  expect_error(ds.predict.gamlss(object=model_e3_bw.DS, newdata="nodata"), "Please provide a valid data.frame as newdata that can be used to derive the predictions from the ds.gamlss.object.", fixed=TRUE)
  expect_error(ds.predict.gamlss(object=model_e3_bw.DS, newdata=newdata, what="nu"), "The parameter what=nu is not among the distribution parameters from the family of the ds.gamlss.object mu, sigma.", fixed=TRUE)
  expect_error(ds.predict.gamlss(object=model_e3_bw.DS, newdata=newdata, what="mu", type="notype"), "The specified type notype should either be 'link' or 'response'.", fixed=TRUE)
})

#
# Done
#

disconnect.studies.dataset.gamlss()


# Test-smk-ds.gamlss ################################################

#
# Set up
#

connect.studies.dataset.gamlss(list("e3_bw", "e3_gac_None", "hs_zbmi_who", "hs_child_age_None",
                                    "h_mbmi_None", "hs_correct_raven", "hs_wgtgain_None"))

#
# Tests
#

test_that("predictions_parametric_gamlss_normal_dist", {
  model_e3_bw.DS <- ds.gamlss(formula = D$e3_bw ~ D$e3_gac_None,
                              data = 'D', family = 'NO()')
  model_e3_bw <- gamlss::gamlss(formula = e3_bw ~ e3_gac_None,
                                data = ds.test_env$local.values)
  newdata <- data.frame(xvar = seq(28, 44, by=0.1))
  names(newdata) <- "e3_gac_None"
  response <- gamlss::predictAll(model_e3_bw, newdata=newdata, type="response", output="data.frame")
  mu.response.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="mu", type="response")
  sigma.response.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="sigma", type="response")
  link <- gamlss::predictAll(model_e3_bw, newdata=newdata, type="link", output="data.frame")
  mu.link.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="mu", type="link")
  sigma.link.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="sigma", type="link")
  expect_length(mu.response.ds, nrow(newdata))
  expect_length(sigma.response.ds, nrow(newdata))
  expect_length(mu.link.ds, nrow(newdata))
  expect_length(sigma.link.ds, nrow(newdata))
  expect_equal(as.vector(mu.response.ds), response$mu, tolerance=1e-07)
  expect_equal(as.vector(sigma.response.ds), response$sigma, tolerance=1e-07)
  expect_equal(as.vector(mu.link.ds), link$mu, tolerance=1e-07)
  expect_equal(as.vector(sigma.link.ds), link$sigma, tolerance=1e-07)
})

test_that("predictions_pb_gamlss_bcpe_dist", {
  model_e3_bw.DS <- ds.gamlss(formula = D$e3_bw ~ pb(D$e3_gac_None),
                              sigma.formula = D$e3_bw ~ pb(D$e3_gac_None),
                              nu.formula = D$e3_bw ~ pb(D$e3_gac_None),
                              tau.formula = D$e3_bw ~ pb(D$e3_gac_None),
                              min.values = min(ds.test_env$local.values$e3_gac_None), 
                              max.values = max(ds.test_env$local.values$e3_gac_None), 
                              min.max.names = "e3_gac_None",
                              data = 'D', family = 'BCPE()')
  pb <- getFromNamespace("pb", "gamlss")
  gamlss.control <- getFromNamespace("gamlss.control", "gamlss")
  glim.control <- getFromNamespace("glim.control", "gamlss")
  model_e3_bw <- gamlss::gamlss(formula = e3_bw ~ pb(e3_gac_None),
                                sigma.formula = e3_bw ~ pb(e3_gac_None),
                                nu.formula = e3_bw ~ pb(e3_gac_None),
                                tau.formula = e3_bw ~ pb(e3_gac_None),
                                data = ds.test_env$local.values,
                                family = gamlss.dist::BCPE())
  newdata <- data.frame(xvar = seq(28, 44, by=0.1))
  names(newdata) <- "e3_gac_None"
  response <- gamlss::predictAll(model_e3_bw, newdata=newdata, type="response", output="data.frame")
  mu.response.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="mu", type="response")
  sigma.response.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="sigma", type="response")
  nu.response.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="nu", type="response")
  tau.response.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="tau", type="response")
  link <- gamlss::predictAll(model_e3_bw, newdata=newdata, type="link", output="data.frame")
  mu.link.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="mu", type="link")
  sigma.link.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="sigma", type="link")
  nu.link.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="nu", type="link")
  tau.link.ds <- ds.predict.gamlss(model_e3_bw.DS, newdata, what="tau", type="link")
  expect_length(mu.response.ds, nrow(newdata))
  expect_length(sigma.response.ds, nrow(newdata))
  expect_length(nu.response.ds, nrow(newdata))
  expect_length(tau.response.ds, nrow(newdata))
  expect_length(mu.link.ds, nrow(newdata))
  expect_length(sigma.link.ds, nrow(newdata))
  expect_length(nu.link.ds, nrow(newdata))
  expect_length(tau.link.ds, nrow(newdata))
  expect_equal(as.vector(mu.response.ds), response$mu, tolerance=1e-03)
  expect_equal(as.vector(sigma.response.ds), response$sigma, tolerance=1e-03)
  expect_equal(as.vector(nu.response.ds), response$nu, tolerance=1e-03)
  expect_equal(as.vector(tau.response.ds), response$tau, tolerance=1e-03)
  expect_equal(as.vector(mu.link.ds), link$mu, tolerance=1e-03)
  expect_equal(as.vector(sigma.link.ds), link$sigma, tolerance=1e-03)
  expect_equal(as.vector(nu.link.ds), link$nu, tolerance=1e-03)
  expect_equal(as.vector(tau.link.ds), link$tau, tolerance=1e-03)

})

#
# Done
#

disconnect.studies.dataset.gamlss()
