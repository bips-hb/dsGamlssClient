# Test-arg-ds.gamlss ################################################

#
# Set up
#

connect.studies.dataset.gamlss(list("e3_bw", "e3_gac_None"))

#
# Tests
#

test_that("input_errors", {
  expect_error(ds.gamlss(), "Please provide a valid formula!", fixed=TRUE)
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, sigma.formula="notaformula"), "sigma.formula should be a formula object")
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, family="BE()"), "Argument 'family' must be either 'NO()', 'NO2()', 'BCCG()' or 'BCPE()'", fixed=TRUE)
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, data = "nodata"), "The input object nodata is not defined in server1, server2, server3!")
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, min.values="notnumeric"), "min.values should be numeric")
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, max.values="notnumeric"), "max.values should be numeric")
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, min.max.names=14), "min.max.names should be character")
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, min.values=c(12,20), max.values=c(30,35), min.max.names=c("smoothingvar1")), "The length 1 of min.max.names does not match the length 2 of min.values and the length 2 of max.values.")
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, checks="notlogical"), "checks should be logical TRUE or FALSE")
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, mu.fix="notlogical"), "mu.fix should be logical TRUE or FALSE")
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, control=c(1,2)), "control should be a numeric vector of length 7")
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, i.control=c(1,2)), "i.control should be a numeric vector of length 4")
  expect_error(ds.gamlss(formula=e3_bw ~ e3_gac_None, autostep="notlogical"), "autostep should be logical TRUE or FALSE")
})

test_that("formula_variants", {
  
})

#
# Done
#

disconnect.studies.dataset.gamlss()


# Test-discctrl-ds.gamlss ################################################

#
# Set up
#

connect.studies.dataset.gamlss(list("e3_bw", "e3_gac_None"))

#
# Tests
#

test_that("discctrl_gamlss", {
  expect_message(ds.gamlss(formula = e3_bw ~ pb(e3_gac_None), 
                         sigma.formula = e3_bw ~ pb(e3_gac_None),
                         nu.formula = e3_bw ~ pb(e3_gac_None),
                         tau.formula = e3_bw ~ pb(e3_gac_None),
                         data = "D_red", 
                         family = "BCPE()"), 
               "MODEL FITTING TERMINATED AT FIRST ITERATION:")
  
  #check output with suppress messages
  model_e3_bw.DS <- suppressMessages(ds.gamlss(formula = e3_bw ~ pb(e3_gac_None), 
                                               sigma.formula = e3_bw ~ pb(e3_gac_None),
                                               nu.formula = e3_bw ~ pb(e3_gac_None),
                                               tau.formula = e3_bw ~ pb(e3_gac_None),
                                               data = "D_red", 
                                               family = "BCPE()"))
  expect_length(model_e3_bw.DS, 7)
  expect_equal(model_e3_bw.DS$y.vector.error, 
               matrix(c(0,0,0), dimnames=list(c("server1", "server2", "server3"), c("Y VECTOR"))))
  expect_equal(model_e3_bw.DS$mu.x.matrix.error, 
               matrix(rep(0, times=6), ncol=2, dimnames = list(c("server1", "server2", "server3"))))
  expect_equal(model_e3_bw.DS$sigma.x.matrix.error, 
               matrix(rep(0, times=6), ncol=2, dimnames = list(c("server1", "server2", "server3"))))
  expect_equal(model_e3_bw.DS$nu.x.matrix.error, 
               matrix(rep(0, times=6), ncol=2, dimnames = list(c("server1", "server2", "server3"))))
  expect_equal(model_e3_bw.DS$tau.x.matrix.error, 
               matrix(rep(0, times=6), ncol=2, dimnames = list(c("server1", "server2", "server3"))))
  expect_equal(model_e3_bw.DS$gamlss.overparameterized, 
               matrix(c(1,1,1), dimnames=list(c("server1", "server2", "server3"), c("MODEL OVERPARAMETERIZED"))))
  expect_equal(model_e3_bw.DS$errorMessage, 
               matrix(rep("Study data or applied model invalid for this source", times=3), 
                      dimnames=list(c("server1", "server2", "server3"), c("ERROR MESSAGES"))))
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

test_that("output_parametric_gamlss_normal_dist", {
  model_e3_bw.DS <- ds.gamlss(formula = e3_bw ~ e3_gac_None,
                              data = 'D', family = 'NO()')
  model_e3_bw <- gamlss::gamlss(formula = e3_bw ~ e3_gac_None,
                                data = ds.test_env$local.values)
  expect_length(model_e3_bw.DS, length(model_e3_bw)+1)
  expect_equal(model_e3_bw.DS$family, model_e3_bw$family)
  expect_equal(model_e3_bw.DS$parameters, model_e3_bw$parameters)
  expect_equal(model_e3_bw.DS$call, quote(ds.gamlss(formula = e3_bw ~ e3_gac_None, family = "NO()", data = "D")))
  expect_equal(model_e3_bw.DS$y, "The response variable is not disclosed!")
  expect_equal(model_e3_bw.DS$control, model_e3_bw$control)
  expect_equal(model_e3_bw.DS$weights, model_e3_bw$weights)
  expect_equal(model_e3_bw.DS$G.deviance, model_e3_bw$G.deviance, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$N, model_e3_bw$N)
  expect_equal(model_e3_bw.DS$rqres, model_e3_bw$rqres)
  expect_equal(model_e3_bw.DS$iter, model_e3_bw$iter)
  expect_equal(model_e3_bw.DS$type, model_e3_bw$type)
  expect_equal(model_e3_bw.DS$method, model_e3_bw$method)
  expect_equal(model_e3_bw.DS$contrasts, model_e3_bw$contrasts)
  expect_equal(model_e3_bw.DS$converged, model_e3_bw$converged)
  expect_equal(model_e3_bw.DS$residuals, "The residuals of the model are not disclosed!")
  expect_equal(model_e3_bw.DS$noObs, model_e3_bw$noObs)
  expect_equal(model_e3_bw.DS$mu.fv, "The fitted values of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.lp, "The linear predictors of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.wv, "The working variable of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.wt, "The working weights of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.link, model_e3_bw$mu.link)
  expect_equal(model_e3_bw.DS$mu.x, "The design matrix of the mu model is not disclosed!")
  expect_equal(model_e3_bw.DS$mu.qr, "The QR decomposition of the mu model is not disclosed!")
  expect_equal(model_e3_bw.DS$mu.coefficients, model_e3_bw$mu.coefficients, tolerance=1e-07, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$mu.offset, model_e3_bw$mu.offset)
  expect_equal(model_e3_bw.DS$mu.xlevels, model_e3_bw$mu.xlevels)
  expect_equal(model_e3_bw.DS$mu.formula, e3_bw ~ e3_gac_None)
  expect_equal(model_e3_bw.DS$mu.df, model_e3_bw$mu.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$mu.nl.df, model_e3_bw$mu.nl.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$mu.pen, model_e3_bw$mu.pen, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$df.fit, model_e3_bw$df.fit, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$pen, model_e3_bw$pen, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$df.residual, model_e3_bw$df.residual, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$sigma.fv, "The fitted values of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.lp, "The linear predictors of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.wv, "The working variable of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.wt, "The working weights of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.link, model_e3_bw$sigma.link)
  expect_equal(model_e3_bw.DS$sigma.x, "The design matrix of the sigma model is not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.qr, "The QR decomposition of the sigma model is not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.coefficients, model_e3_bw$sigma.coefficients, tolerance=1e-07, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$sigma.offset, model_e3_bw$sigma.offset)
  expect_equal(model_e3_bw.DS$sigma.xlevels, model_e3_bw$sigma.xlevels)
  expect_equal(model_e3_bw.DS$sigma.df, model_e3_bw$sigma.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$sigma.nl.df, model_e3_bw$sigma.nl.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$sigma.pen, model_e3_bw$sigma.pen, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$P.deviance, model_e3_bw$P.deviance, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$aic, model_e3_bw$aic, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$sbc, model_e3_bw$sbc, tolerance=1e-07)
  expect_equal(DSI::datashield.symbols(conns=ds.test_env$connections)[[1]], c("D", "D_red"))
})

test_that("output_parametric_gamlss_bcpe_dist", {
  model_e3_bw.DS <- ds.gamlss(formula = e3_bw ~ e3_gac_None,
                              sigma.formula = e3_bw ~ e3_gac_None,
                              nu.formula = e3_bw ~ e3_gac_None,
                              tau.formula = e3_bw ~ e3_gac_None,
                              data = 'D', family = 'BCPE()')
  model_e3_bw <- gamlss::gamlss(formula = e3_bw ~ e3_gac_None,
                                sigma.formula = e3_bw ~ e3_gac_None,
                                nu.formula = e3_bw ~ e3_gac_None,
                                tau.formula = e3_bw ~ e3_gac_None,
                                data = ds.test_env$local.values,
                                family = gamlss.dist::BCPE())
  expect_length(model_e3_bw.DS, length(model_e3_bw)+1)
  expect_equal(model_e3_bw.DS$family, model_e3_bw$family)
  expect_equal(model_e3_bw.DS$parameters, model_e3_bw$parameters)
  expect_equal(model_e3_bw.DS$call, quote(ds.gamlss(formula = e3_bw ~ e3_gac_None, 
                                                    sigma.formula = e3_bw ~ e3_gac_None,
                                                    nu.formula = e3_bw ~ e3_gac_None,
                                                    tau.formula = e3_bw ~ e3_gac_None,
                                                    family = "BCPE()", data = "D")))
  expect_equal(model_e3_bw.DS$y, "The response variable is not disclosed!")
  expect_equal(model_e3_bw.DS$control, model_e3_bw$control)
  expect_equal(model_e3_bw.DS$weights, model_e3_bw$weights)
  expect_equal(model_e3_bw.DS$G.deviance, model_e3_bw$G.deviance, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$N, model_e3_bw$N)
  expect_equal(model_e3_bw.DS$rqres, model_e3_bw$rqres)
  expect_equal(model_e3_bw.DS$iter, model_e3_bw$iter)
  expect_equal(model_e3_bw.DS$type, model_e3_bw$type)
  expect_equal(model_e3_bw.DS$method, model_e3_bw$method)
  expect_equal(model_e3_bw.DS$contrasts, model_e3_bw$contrasts)
  expect_equal(model_e3_bw.DS$converged, model_e3_bw$converged)
  expect_equal(model_e3_bw.DS$residuals, "The residuals of the model are not disclosed!")
  expect_equal(model_e3_bw.DS$noObs, model_e3_bw$noObs)
  expect_equal(model_e3_bw.DS$mu.fv, "The fitted values of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.lp, "The linear predictors of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.wv, "The working variable of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.wt, "The working weights of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.link, model_e3_bw$mu.link)
  expect_equal(model_e3_bw.DS$mu.x, "The design matrix of the mu model is not disclosed!")
  expect_equal(model_e3_bw.DS$mu.qr, "The QR decomposition of the mu model is not disclosed!")
  expect_equal(model_e3_bw.DS$mu.coefficients, model_e3_bw$mu.coefficients, tolerance=1e-07, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$mu.offset, model_e3_bw$mu.offset)
  expect_equal(model_e3_bw.DS$mu.xlevels, model_e3_bw$mu.xlevels)
  expect_equal(model_e3_bw.DS$mu.formula, e3_bw ~ e3_gac_None)
  expect_equal(model_e3_bw.DS$mu.df, model_e3_bw$mu.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$mu.nl.df, model_e3_bw$mu.nl.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$mu.pen, model_e3_bw$mu.pen, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$df.fit, model_e3_bw$df.fit, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$pen, model_e3_bw$pen, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$df.residual, model_e3_bw$df.residual, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$sigma.fv, "The fitted values of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.lp, "The linear predictors of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.wv, "The working variable of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.wt, "The working weights of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.link, model_e3_bw$sigma.link)
  expect_equal(model_e3_bw.DS$sigma.x, "The design matrix of the sigma model is not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.qr, "The QR decomposition of the sigma model is not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.coefficients, model_e3_bw$sigma.coefficients, tolerance=1e-07, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$sigma.offset, model_e3_bw$sigma.offset)
  expect_equal(model_e3_bw.DS$sigma.xlevels, model_e3_bw$sigma.xlevels)
  expect_equal(model_e3_bw.DS$sigma.formula, e3_bw ~ e3_gac_None)
  expect_equal(model_e3_bw.DS$sigma.df, model_e3_bw$sigma.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$sigma.nl.df, model_e3_bw$sigma.nl.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$sigma.pen, model_e3_bw$sigma.pen, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$nu.fv, "The fitted values of the nu model are not disclosed!")
  expect_equal(model_e3_bw.DS$nu.lp, "The linear predictors of the nu model are not disclosed!")
  expect_equal(model_e3_bw.DS$nu.wv, "The working variable of the nu model are not disclosed!")
  expect_equal(model_e3_bw.DS$nu.wt, "The working weights of the nu model are not disclosed!")
  expect_equal(model_e3_bw.DS$nu.link, model_e3_bw$nu.link)
  expect_equal(model_e3_bw.DS$nu.x, "The design matrix of the nu model is not disclosed!")
  expect_equal(model_e3_bw.DS$nu.qr, "The QR decomposition of the nu model is not disclosed!")
  expect_equal(model_e3_bw.DS$nu.coefficients, model_e3_bw$nu.coefficients, tolerance=1e-07, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$nu.offset, model_e3_bw$nu.offset)
  expect_equal(model_e3_bw.DS$nu.xlevels, model_e3_bw$nu.xlevels)
  expect_equal(model_e3_bw.DS$nu.formula, e3_bw ~ e3_gac_None)
  expect_equal(model_e3_bw.DS$nu.df, model_e3_bw$nu.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$nu.nl.df, model_e3_bw$nu.nl.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$nu.pen, model_e3_bw$nu.pen, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$tau.fv, "The fitted values of the tau model are not disclosed!")
  expect_equal(model_e3_bw.DS$tau.lp, "The linear predictors of the tau model are not disclosed!")
  expect_equal(model_e3_bw.DS$tau.wv, "The working variable of the tau model are not disclosed!")
  expect_equal(model_e3_bw.DS$tau.wt, "The working weights of the tau model are not disclosed!")
  expect_equal(model_e3_bw.DS$tau.link, model_e3_bw$tau.link)
  expect_equal(model_e3_bw.DS$tau.x, "The design matrix of the tau model is not disclosed!")
  expect_equal(model_e3_bw.DS$tau.qr, "The QR decomposition of the tau model is not disclosed!")
  expect_equal(model_e3_bw.DS$tau.coefficients, model_e3_bw$tau.coefficients, tolerance=1e-07, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$tau.offset, model_e3_bw$tau.offset)
  expect_equal(model_e3_bw.DS$tau.xlevels, model_e3_bw$tau.xlevels)
  expect_equal(model_e3_bw.DS$tau.formula, e3_bw ~ e3_gac_None)
  expect_equal(model_e3_bw.DS$tau.df, model_e3_bw$tau.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$tau.nl.df, model_e3_bw$tau.nl.df, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$tau.pen, model_e3_bw$tau.pen, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$P.deviance, model_e3_bw$P.deviance, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$aic, model_e3_bw$aic, tolerance=1e-07)
  expect_equal(model_e3_bw.DS$sbc, model_e3_bw$sbc, tolerance=1e-07)
  expect_equal(DSI::datashield.symbols(conns=ds.test_env$connections)[[1]], c("D", "D_red"))
})

test_that("output_pb_gamlss_bcpe_dist", {
  model_e3_bw.DS <- ds.gamlss(formula = e3_bw ~ pb(e3_gac_None),
                              sigma.formula = e3_bw ~ pb(e3_gac_None),
                              nu.formula = e3_bw ~ pb(e3_gac_None),
                              tau.formula = e3_bw ~ pb(e3_gac_None),
                              min.values = min(ds.test_env$local.values$e3_gac_None), 
                              max.values = max(ds.test_env$local.values$e3_gac_None), 
                              min.max.names = "e3_gac_None",
                              data = 'D', family = 'BCPE()')
  pb <- utils::getFromNamespace("pb", "gamlss")
  gamlss.control <- utils::getFromNamespace("gamlss.control", "gamlss")
  glim.control <- utils::getFromNamespace("glim.control", "gamlss")
  model_e3_bw <- gamlss::gamlss(formula = e3_bw ~ pb(e3_gac_None),
                                sigma.formula = e3_bw ~ pb(e3_gac_None),
                                nu.formula = e3_bw ~ pb(e3_gac_None),
                                tau.formula = e3_bw ~ pb(e3_gac_None),
                                data = ds.test_env$local.values,
                                family = gamlss.dist::BCPE())
  expect_length(model_e3_bw.DS, length(model_e3_bw)+1)
  expect_equal(model_e3_bw.DS$family, model_e3_bw$family)
  expect_equal(model_e3_bw.DS$parameters, model_e3_bw$parameters)
  expect_equal(model_e3_bw.DS$y, "The response variable is not disclosed!")
  expect_equal(model_e3_bw.DS$control, model_e3_bw$control)
  expect_equal(model_e3_bw.DS$weights, model_e3_bw$weights)
  expect_equal(model_e3_bw.DS$G.deviance, model_e3_bw$G.deviance, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$N, model_e3_bw$N)
  expect_equal(model_e3_bw.DS$rqres, model_e3_bw$rqres)
  expect_equal(model_e3_bw.DS$iter, model_e3_bw$iter)
  expect_equal(model_e3_bw.DS$type, model_e3_bw$type)
  expect_equal(model_e3_bw.DS$method, model_e3_bw$method)
  expect_equal(model_e3_bw.DS$contrasts, model_e3_bw$contrasts)
  expect_equal(model_e3_bw.DS$converged, model_e3_bw$converged)
  expect_equal(model_e3_bw.DS$residuals, "The residuals of the model are not disclosed!")
  expect_equal(model_e3_bw.DS$noObs, model_e3_bw$noObs)
  expect_equal(model_e3_bw.DS$mu.fv, "The fitted values of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.lp, "The linear predictors of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.wv, "The working variable of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.wt, "The working weights of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.link, model_e3_bw$mu.link)
  expect_equal(model_e3_bw.DS$mu.x, "The design matrix of the mu model is not disclosed!")
  expect_equal(model_e3_bw.DS$mu.qr, "The QR decomposition of the mu model is not disclosed!")
  expect_equal(model_e3_bw.DS$mu.coefficients, model_e3_bw$mu.coefficients, tolerance=1e-03, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$mu.offset, model_e3_bw$mu.offset)
  expect_equal(model_e3_bw.DS$mu.xlevels, model_e3_bw$mu.xlevels)
  expect_equal(model_e3_bw.DS$mu.formula, e3_bw ~ pb(e3_gac_None))
  expect_equal(model_e3_bw.DS$mu.df, model_e3_bw$mu.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.nl.df, model_e3_bw$mu.nl.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.s, "The smoothing fitted values of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.var, "The variances for the smoothing fitted values of the mu model are not disclosed!")
  expect_length(model_e3_bw.DS$mu.coefSmo[[1]], length(model_e3_bw$mu.coefSmo[[1]]))
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$coef, model_e3_bw$mu.coefSmo[[1]]$coef, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$fv, "The smoothing fitted values of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$lambda, model_e3_bw$mu.coefSmo[[1]]$lambda, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$edf, model_e3_bw$mu.coefSmo[[1]]$edf, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$sigb, model_e3_bw$mu.coefSmo[[1]]$sigb, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$sige, model_e3_bw$mu.coefSmo[[1]]$sige, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$method, model_e3_bw$mu.coefSmo[[1]]$method)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$knots, model_e3_bw$mu.coefSmo[[1]]$knots, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$fun, "The function for the knots of the mu model is not disclosed!")
  expect_equal(model_e3_bw.DS$mu.pen, model_e3_bw$mu.pen, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$df.fit, model_e3_bw$df.fit, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$pen, model_e3_bw$pen, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$df.residual, model_e3_bw$df.residual, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.fv, "The fitted values of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.lp, "The linear predictors of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.wv, "The working variable of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.wt, "The working weights of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.link, model_e3_bw$sigma.link)
  expect_equal(model_e3_bw.DS$sigma.x, "The design matrix of the sigma model is not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.qr, "The QR decomposition of the sigma model is not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.coefficients, model_e3_bw$sigma.coefficients, tolerance=1e-03, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$sigma.offset, model_e3_bw$sigma.offset)
  expect_equal(model_e3_bw.DS$sigma.xlevels, model_e3_bw$sigma.xlevels)
  expect_equal(model_e3_bw.DS$sigma.formula, e3_bw ~ pb(e3_gac_None))
  expect_equal(model_e3_bw.DS$sigma.df, model_e3_bw$sigma.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.nl.df, model_e3_bw$sigma.nl.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.s, "The smoothing fitted values of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.var, "The variances for the smoothing fitted values of the sigma model are not disclosed!")
  expect_length(model_e3_bw.DS$sigma.coefSmo[[1]], length(model_e3_bw$sigma.coefSmo[[1]]))
  expect_equal(model_e3_bw.DS$sigma.coefSmo[[1]]$coef, model_e3_bw$sigma.coefSmo[[1]]$coef, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.coefSmo[[1]]$fv, "The smoothing fitted values of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.coefSmo[[1]]$lambda, model_e3_bw$sigma.coefSmo[[1]]$lambda, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.coefSmo[[1]]$edf, model_e3_bw$sigma.coefSmo[[1]]$edf, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.coefSmo[[1]]$sigb, model_e3_bw$sigma.coefSmo[[1]]$sigb, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.coefSmo[[1]]$sige, model_e3_bw$sigma.coefSmo[[1]]$sige, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.coefSmo[[1]]$method, model_e3_bw$sigma.coefSmo[[1]]$method)
  expect_equal(model_e3_bw.DS$sigma.coefSmo[[1]]$knots, model_e3_bw$sigma.coefSmo[[1]]$knots, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.coefSmo[[1]]$fun, "The function for the knots of the sigma model is not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.pen, model_e3_bw$sigma.pen, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$nu.fv, "The fitted values of the nu model are not disclosed!")
  expect_equal(model_e3_bw.DS$nu.lp, "The linear predictors of the nu model are not disclosed!")
  expect_equal(model_e3_bw.DS$nu.wv, "The working variable of the nu model are not disclosed!")
  expect_equal(model_e3_bw.DS$nu.wt, "The working weights of the nu model are not disclosed!")
  expect_equal(model_e3_bw.DS$nu.link, model_e3_bw$nu.link)
  expect_equal(model_e3_bw.DS$nu.x, "The design matrix of the nu model is not disclosed!")
  expect_equal(model_e3_bw.DS$nu.qr, "The QR decomposition of the nu model is not disclosed!")
  expect_equal(model_e3_bw.DS$nu.coefficients, model_e3_bw$nu.coefficients, tolerance=1e-03, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$nu.offset, model_e3_bw$nu.offset)
  expect_equal(model_e3_bw.DS$nu.xlevels, model_e3_bw$nu.xlevels)
  expect_equal(model_e3_bw.DS$nu.formula, e3_bw ~ pb(e3_gac_None))
  expect_equal(model_e3_bw.DS$nu.df, model_e3_bw$nu.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$nu.nl.df, model_e3_bw$nu.nl.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$nu.s, "The smoothing fitted values of the nu model are not disclosed!")
  expect_equal(model_e3_bw.DS$nu.var, "The variances for the smoothing fitted values of the nu model are not disclosed!")
  expect_length(model_e3_bw.DS$nu.coefSmo[[1]], length(model_e3_bw$nu.coefSmo[[1]]))
  expect_equal(model_e3_bw.DS$nu.coefSmo[[1]]$coef, model_e3_bw$nu.coefSmo[[1]]$coef, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$nu.coefSmo[[1]]$fv, "The smoothing fitted values of the nu model are not disclosed!")
  expect_equal(model_e3_bw.DS$nu.coefSmo[[1]]$lambda, model_e3_bw$nu.coefSmo[[1]]$lambda, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$nu.coefSmo[[1]]$edf, model_e3_bw$nu.coefSmo[[1]]$edf, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$nu.coefSmo[[1]]$sigb, model_e3_bw$nu.coefSmo[[1]]$sigb, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$nu.coefSmo[[1]]$sige, model_e3_bw$nu.coefSmo[[1]]$sige, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$nu.coefSmo[[1]]$method, model_e3_bw$nu.coefSmo[[1]]$method)
  expect_equal(model_e3_bw.DS$nu.coefSmo[[1]]$knots, model_e3_bw$nu.coefSmo[[1]]$knots, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$nu.coefSmo[[1]]$fun, "The function for the knots of the nu model is not disclosed!")
  expect_equal(model_e3_bw.DS$nu.pen, model_e3_bw$nu.pen, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$tau.fv, "The fitted values of the tau model are not disclosed!")
  expect_equal(model_e3_bw.DS$tau.lp, "The linear predictors of the tau model are not disclosed!")
  expect_equal(model_e3_bw.DS$tau.wv, "The working variable of the tau model are not disclosed!")
  expect_equal(model_e3_bw.DS$tau.wt, "The working weights of the tau model are not disclosed!")
  expect_equal(model_e3_bw.DS$tau.link, model_e3_bw$tau.link)
  expect_equal(model_e3_bw.DS$tau.x, "The design matrix of the tau model is not disclosed!")
  expect_equal(model_e3_bw.DS$tau.qr, "The QR decomposition of the tau model is not disclosed!")
  expect_equal(model_e3_bw.DS$tau.coefficients, model_e3_bw$tau.coefficients, tolerance=1e-03, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$tau.offset, model_e3_bw$tau.offset)
  expect_equal(model_e3_bw.DS$tau.xlevels, model_e3_bw$tau.xlevels)
  expect_equal(model_e3_bw.DS$tau.formula, e3_bw ~ pb(e3_gac_None))
  expect_equal(model_e3_bw.DS$tau.df, model_e3_bw$tau.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$tau.nl.df, model_e3_bw$tau.nl.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$tau.s, "The smoothing fitted values of the tau model are not disclosed!")
  expect_equal(model_e3_bw.DS$tau.var, "The variances for the smoothing fitted values of the tau model are not disclosed!")
  expect_length(model_e3_bw.DS$tau.coefSmo[[1]], length(model_e3_bw$tau.coefSmo[[1]]))
  expect_equal(model_e3_bw.DS$tau.coefSmo[[1]]$coef, model_e3_bw$tau.coefSmo[[1]]$coef, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$tau.coefSmo[[1]]$fv, "The smoothing fitted values of the tau model are not disclosed!")
  expect_equal(model_e3_bw.DS$tau.coefSmo[[1]]$lambda, model_e3_bw$tau.coefSmo[[1]]$lambda, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$tau.coefSmo[[1]]$edf, model_e3_bw$tau.coefSmo[[1]]$edf, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$tau.coefSmo[[1]]$sigb, model_e3_bw$tau.coefSmo[[1]]$sigb, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$tau.coefSmo[[1]]$sige, model_e3_bw$tau.coefSmo[[1]]$sige, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$tau.coefSmo[[1]]$method, model_e3_bw$tau.coefSmo[[1]]$method)
  expect_equal(model_e3_bw.DS$tau.coefSmo[[1]]$knots, model_e3_bw$tau.coefSmo[[1]]$knots, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$tau.coefSmo[[1]]$fun, "The function for the knots of the tau model is not disclosed!")
  expect_equal(model_e3_bw.DS$tau.pen, model_e3_bw$tau.pen, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$P.deviance, model_e3_bw$P.deviance, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$aic, model_e3_bw$aic, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sbc, model_e3_bw$sbc, tolerance=1e-03)
  expect_equal(DSI::datashield.symbols(conns=ds.test_env$connections)[[1]], c("D", "D_red"))
})

test_that("output_pb_gamlss_formula_variant", {
  model_e3_bw.DS <- ds.gamlss(formula = D$e3_bw ~ pb(D$e3_gac_None),
                              sigma.formula = D$e3_bw ~ D$e3_gac_None,
                              min.values = min(ds.test_env$local.values$e3_gac_None), 
                              max.values = max(ds.test_env$local.values$e3_gac_None), 
                              min.max.names = "e3_gac_None",
                              family = 'NO()')
  pb <- utils::getFromNamespace("pb", "gamlss")
  gamlss.control <- utils::getFromNamespace("gamlss.control", "gamlss")
  glim.control <- utils::getFromNamespace("glim.control", "gamlss")
  model_e3_bw <- gamlss::gamlss(formula = e3_bw ~ pb(e3_gac_None),
                                sigma.formula = e3_bw ~ e3_gac_None,
                                data = ds.test_env$local.values,
                                family = gamlss.dist::NO())
  expect_length(model_e3_bw.DS, length(model_e3_bw)+1)
  expect_equal(model_e3_bw.DS$family, model_e3_bw$family)
  expect_equal(model_e3_bw.DS$parameters, model_e3_bw$parameters)
  expect_equal(model_e3_bw.DS$y, "The response variable is not disclosed!")
  expect_equal(model_e3_bw.DS$control, model_e3_bw$control)
  expect_equal(model_e3_bw.DS$weights, model_e3_bw$weights)
  expect_equal(model_e3_bw.DS$G.deviance, model_e3_bw$G.deviance, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$N, model_e3_bw$N)
  expect_equal(model_e3_bw.DS$rqres, model_e3_bw$rqres)
  expect_equal(model_e3_bw.DS$iter, model_e3_bw$iter)
  expect_equal(model_e3_bw.DS$type, model_e3_bw$type)
  expect_equal(model_e3_bw.DS$method, model_e3_bw$method)
  expect_equal(model_e3_bw.DS$contrasts, model_e3_bw$contrasts)
  expect_equal(model_e3_bw.DS$converged, model_e3_bw$converged)
  expect_equal(model_e3_bw.DS$residuals, "The residuals of the model are not disclosed!")
  expect_equal(model_e3_bw.DS$noObs, model_e3_bw$noObs)
  expect_equal(model_e3_bw.DS$mu.fv, "The fitted values of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.lp, "The linear predictors of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.wv, "The working variable of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.wt, "The working weights of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.link, model_e3_bw$mu.link)
  expect_equal(model_e3_bw.DS$mu.x, "The design matrix of the mu model is not disclosed!")
  expect_equal(model_e3_bw.DS$mu.qr, "The QR decomposition of the mu model is not disclosed!")
  expect_equal(model_e3_bw.DS$mu.coefficients, model_e3_bw$mu.coefficients, tolerance=1e-03, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$mu.offset, model_e3_bw$mu.offset)
  expect_equal(model_e3_bw.DS$mu.xlevels, model_e3_bw$mu.xlevels)
  expect_equal(model_e3_bw.DS$mu.formula, D$e3_bw ~ pb(D$e3_gac_None))
  expect_equal(model_e3_bw.DS$mu.df, model_e3_bw$mu.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.nl.df, model_e3_bw$mu.nl.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.s, "The smoothing fitted values of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.var, "The variances for the smoothing fitted values of the mu model are not disclosed!")
  expect_length(model_e3_bw.DS$mu.coefSmo[[1]], length(model_e3_bw$mu.coefSmo[[1]]))
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$coef, model_e3_bw$mu.coefSmo[[1]]$coef, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$fv, "The smoothing fitted values of the mu model are not disclosed!")
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$lambda, model_e3_bw$mu.coefSmo[[1]]$lambda, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$edf, model_e3_bw$mu.coefSmo[[1]]$edf, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$sigb, model_e3_bw$mu.coefSmo[[1]]$sigb, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$sige, model_e3_bw$mu.coefSmo[[1]]$sige, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$method, model_e3_bw$mu.coefSmo[[1]]$method)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$knots, model_e3_bw$mu.coefSmo[[1]]$knots, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$mu.coefSmo[[1]]$fun, "The function for the knots of the mu model is not disclosed!")
  expect_equal(model_e3_bw.DS$mu.pen, model_e3_bw$mu.pen, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$df.fit, model_e3_bw$df.fit, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$pen, model_e3_bw$pen, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$df.residual, model_e3_bw$df.residual, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.fv, "The fitted values of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.lp, "The linear predictors of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.wv, "The working variable of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.wt, "The working weights of the sigma model are not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.link, model_e3_bw$sigma.link)
  expect_equal(model_e3_bw.DS$sigma.x, "The design matrix of the sigma model is not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.qr, "The QR decomposition of the sigma model is not disclosed!")
  expect_equal(model_e3_bw.DS$sigma.coefficients, model_e3_bw$sigma.coefficients, tolerance=1e-03, ignore_attr=TRUE)
  expect_equal(model_e3_bw.DS$sigma.offset, model_e3_bw$sigma.offset)
  expect_equal(model_e3_bw.DS$sigma.xlevels, model_e3_bw$sigma.xlevels)
  expect_equal(model_e3_bw.DS$sigma.formula, D$e3_bw ~ D$e3_gac_None)
  expect_equal(model_e3_bw.DS$sigma.df, model_e3_bw$sigma.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sigma.nl.df, model_e3_bw$sigma.nl.df, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$P.deviance, model_e3_bw$P.deviance, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$aic, model_e3_bw$aic, tolerance=1e-03)
  expect_equal(model_e3_bw.DS$sbc, model_e3_bw$sbc, tolerance=1e-03)
  expect_equal(DSI::datashield.symbols(conns=ds.test_env$connections)[[1]], c("D", "D_red"))
})

#
# Done
#

disconnect.studies.dataset.gamlss()
