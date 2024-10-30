#'
#' @title Generalized Additive Models for Location Scale and Shape
#' @description Fits a generalized additive model for location, scale and shape (GAMLSS)
#' on data from a single source or multiple sources on the server side. 
#' @details Fits a generalized additive model for location, scale and shape (GAMLSS)
#' on data from a single source or multiple sources on the server side. In the latter
#' case, the data are co-analysed (when using \code{ds.gamlss})  by using an approach 
#' that is mathematically equivalent to placing all individual-level data from all sources
#' in one central warehouse and analysing those data using the conventional 
#' \code{gamlss()} function in R. For additional details see the help header of gamlss
#' functions in native R gamlss package.
#' 
#' Server functions called: \code{gamlssDS1}, 
#'                          \code{gamlssDS2},
#'                          \code{gamlssDS3},
#'                          \code{gamlssDS4},
#'                          \code{gamlssDS5},
#'                          \code{gamlssDS6},
#'                          \code{gamlssDS7}
#' 
#' @param formula a formula object, with the response on the left of an ~ operator, 
#' and the terms, separated by + operators, on the right. Nonparametric smoothing
#' terms are indicated by pb() for penalised beta splines, cs for smoothing splines, 
#' lo for loess smooth terms and random or ra for random terms, 
#' e.g. 'y~cs(x,df=5)+x1+x2*x3'. 
#' @param sigma.formula a formula object for fitting a model to the sigma parameter,
#' as in the formula above, e.g. sigma.formula='~cs(x,df=5)'.
#' @param nu.formula a formula object for fitting a model to the nu parameter, 
#' e.g. nu.formula='~x'.
#' @param tau.formula a formula object for fitting a model to the tau parameter, 
#' e.g. tau.formula='~cs(x,df=2)'.
#' @param family a gamlss.family object, which is used to define the distribution 
#' and the link functions of the various parameters. The distribution families 
#' supported by gamlss() can be found in gamlss.family. Functions such as 'BI()' 
#' (binomial) produce a family object. Also can be given without the parentheses
#' i.e. 'BI'. Family functions can take arguments, as in 'BI(mu.link=probit)'.
#' @param data a data frame containing the variables occurring in the formula. 
#' If this is missing, the variables should be on the parent environment.
#' @param min.values numeric vector. Allows to specify minimum values for the covariates
#' which are used to determine the knots for pb(). If NULL an anonymized (noisy)
#' minimum is used instead. Default NULL.
#' @param max.values numeric vector. Allows to specify maximum values for the covariates
#' which are used to determine the knots for pb(). If NULL an anonymized (noisy)
#' maximum is used instead. Default NULL.
#' @param min.max.names a character vector giving the names for the minimum and
#' maximum values. Only required if \code{min.values} or \code{max.values} are given. 
#' Default NULL.
#' @param checks logical. If TRUE \code{ds.gamlss} checks the structural integrity 
#' of the model. Default FALSE.
#' @param method a character indicating the algorithm for GAMLSS. Currently only the
#' Rigby and Stasinopolous algorithm (\code{method='RS'}) is implemented. 
#' @param mu.coef.start optional vector of regression coefficients to compute improved 
#' start values for mu. Default NULL.
#' @param sigma.coef.start optional vector of regression coefficients to compute improved
#' start values for sigma. Default NULL.
#' @param nu.coef.start optional vector of regression coefficients to compute improved
#' start values for nu. Default NULL.
#' @param tau.coef.start optional vector of regression coefficients to compute improved
#' start values for tau. Default NULL.
#' @param mu.coef.start.names vector with names for the regression coefficients in 
#' \code{mu.coef.start}. These names are needed to obtain the design matrix to compute the
#' improved start values for mu. If values are given in \code{mu.coef.start} but 
#' \code{mu.coef.start.names} is NULL then the same formula as in \code{formula} is used to
#' obtain the design matrix. Default NULL.
#' @param sigma.coef.start.names vector with names for the regression coefficients in 
#' \code{sigma.coef.start}. These names are needed to obtain the design matrix to compute the
#' improved start values for sigma. If values are given in \code{sigma.coef.start} but 
#' \code{sigma.coef.start.names} is NULL then the same formula as in \code{sigma.formula} 
#' is used to obtain the design matrix. Default NULL.
#' @param nu.coef.start.names vector with names for the regression coefficients in 
#' \code{nu.coef.start}. These names are needed to obtain the design matrix to compute the
#' improved start values for nu. If values are given in \code{nu.coef.start} but 
#' \code{nu.coef.start.names} is NULL then the same formula as in \code{nu.formula} 
#' is used to obtain the design matrix. Default NULL.
#' @param tau.coef.start.names vector with names for the regression coefficients in 
#' \code{tau.coef.start}. These names are needed to obtain the design matrix to compute the
#' improved start values for tau. If values are given in \code{tau.coef.start} but 
#' \code{tau.coef.start.names} is NULL then the same formula as in \code{tau.formula} 
#' is used to obtain the design matrix. Default NULL.
#' @param mu.fix logical, indicating whether the mu parameter should be kept fixed
#' in the fitting processes.
#' @param sigma.fix logical, indicating whether the sigma parameter should be kept
#' fixed in the fitting processes.
#' @param nu.fix logical, indicating whether the nu parameter should be kept fixed 
#' in the fitting processes.
#' @param tau.fix logical, indicating whether the tau parameter should be kept fixed
#' in the fitting processes.
#' @param control this sets the control parameters of the outer iterations algorithm 
#' using the gamlss.control function. This is a vector of 7 numeric values: (i) c.crit 
#' (the convergence criterion for the algorithm), (ii) n.cyc (the number of cycles of 
#' the algorithm), (iii) mu.step (the step length for the parameter mu), (iv) sigma.step 
#' (the step length for the parameter sigma), (v) nu.step (the step length for the
#' parameter nu), (vi) tau.step (the step length for the parameter tau), (vii) gd.tol
#' (global deviance tolerance level). The default values for these 7 parameters are 
#' set to c(0.001, 20, 1, 1, 1, 1, Inf).
#' @param i.control this sets the control parameters of the inner iterations of the 
#' RS algorithm using the glim.control function. This is a vector of 4 numeric values: 
#' (i) cc (the convergence criterion for the algorithm), (ii) cyc (the number of 
#' cycles of the algorithm), (iii) bf.cyc (the number of cycles of the backfitting 
#' algorithm), (iv) bf.tol (the convergence criterion (tolerance level) for the 
#' backfitting algorithm). The default values for these 4 parameters are set to 
#' c(0.001, 50, 30, 0.001).
#' @param autostep logical, indicating whether the steps should be halved automatically 
#' if the new global deviance is greater than the old one. The default is TRUE.
#' @param k  the number of the nearest neighbours for which the mean is calculated
#' to obtain the anonymized quantile residuals. The default is 3.
#' @param datasources  a list of \code{\link[DSI]{DSConnection-class}} 
#' objects obtained after login. If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link[DSI]{datashield.connections_default}}.
#' @return a gamlss object with all components as in the native R gamlss function. 
#' Individual-level information like the components y (the response response) and 
#' residuals (the normalised quantile residuals of the model) are not disclosed to 
#' the client-side.
#' @author Annika Swenne
#' @export

ds.gamlss <- function(formula = NULL, sigma.formula = ~1, nu.formula = ~1, tau.formula = ~1,
                      family = 'NO()', data = NULL, min.values = NULL, max.values = NULL, 
                      min.max.names = NULL, checks = FALSE, method = 'RS', 
                      mu.coef.start = NULL, sigma.coef.start = NULL, 
                      nu.coef.start = NULL, tau.coef.start = NULL, 
                      mu.coef.start.names = NULL, sigma.coef.start.names = NULL, 
                      nu.coef.start.names = NULL, tau.coef.start.names = NULL, 
                      mu.fix = FALSE, sigma.fix = FALSE, nu.fix = FALSE, tau.fix = FALSE, 
                      control = c(0.001, 20, 1, 1, 1, 1, Inf),
                      i.control = c(0.001, 50, 30, 0.001), 
                      autostep = TRUE, k = 3, datasources = NULL){
  
  #**************************************************************************
  # I) Preparation ----
  # Are the required parameters appropriately specified?
  # Create objects that can be transferred to the server
  # Initialization by calling gamlssDS1
  #**************************************************************************
  
  # help function to integrate returned output
  .select <- function(l, field){
    lapply(l, function(obj) {obj[[field]]})
  }
  
  # help function to extract control parameters for pb smoother
  getpbcontrol <- function(args, pattern="pb.control"){
    if (length(grep(pattern="pb.control", args, fixed=TRUE))>0) {
      # control parameters specified
      pb.control <- eval(parse(text=args[grep(pattern="pb.control", x=args, fixed=TRUE)]))
    } else {
      # no control parameters specified - use default
      pb.control <- eval(parse(text="pb.control()"))
    }
  }
  
  # look for DS connections
  if(is.null(datasources)){
    datasources <- DSI::datashield.connections_find()
  }
  
  # ensure datasources is a list of DSConnection-class
  if(!(is.list(datasources) && all(unlist(lapply(datasources, function(d) {methods::is(d,"DSConnection")}))))){
    stop("The 'datasources' were expected to be a list of DSConnection-class objects", call.=FALSE)
  }
  
  # verify that 'formula' was set
  if(is.null(formula)){
    stop(" Please provide a valid formula!", call.=FALSE)
  }
  
  formula <- stats::as.formula(formula)
  sigma.formula <- stats::as.formula(sigma.formula)
  nu.formula <- stats::as.formula(nu.formula)
  tau.formula <- stats::as.formula(tau.formula)
  
  formulatext <- Reduce(paste, deparse(formula))
  sigma.formulatext <- Reduce(paste, deparse(sigma.formula))
  nu.formulatext <- Reduce(paste, deparse(nu.formula))
  tau.formulatext <- Reduce(paste, deparse(tau.formula))
  
  ## check that 'family' was set
  if(is.null(family)){
    stop(" Please provide a valid 'family' argument!", call.=FALSE)
  }
  
  ## if the argument 'data' is set, check that the data frame is defined (i.e. exists) on the server site
  if(!(is.null(data))){
    defined <- isDefined(datasources, data)
  }
  
  ## beginning of optional checks - the process stops if any of these checks fails #
  if(checks){
    message(" -- Verifying the variables in the model")
    # call the function that checks the variables in the formula are defined (exist) on the server site and are not missing at complete
    gamlssChecks(formula=formula, sigma.formula=sigma.formula, nu.formula=nu.formula, tau.formula=tau.formula, 
                 data=data, datasources=datasources)
  }else{
    #message("WARNING:'checks' is set to FALSE; variables in the model are not checked and error messages may not be intelligible!")
  }
  
  ## Outer iteration count (before assignment of beta.vect):
  # Iterations need to be counted. Start off with the count at 0
  # and increment by 1 at each new iteration
  outer.iteration.count <- 0
  
  ## number of 'valid' studies (those that passed the checks)
  numstudies <- length(datasources)
  
  ## Create strings for transfer to the server
  # special symbols (brackets, etc. cannot be transmitted to the server)
  # --> convert them to plain text and reconvert them back on the server side
  formula.trans <- gsub("(", "left_parenthesis", formulatext, fixed = TRUE)
  formula.trans <- gsub(")", "right_parenthesis", formula.trans, fixed = TRUE)
  formula.trans <- gsub("~", "tilde_symbol", formula.trans, fixed = TRUE)
  formula.trans <- gsub("=", "equal_symbol", formula.trans, fixed = TRUE)
  formula.trans <- gsub(",", "comma_symbol", formula.trans, fixed = TRUE)
  formula.trans <- gsub("*", "asterisk_symbol", formula.trans, fixed = TRUE)
  formula.trans <- gsub("^", "caret_symbol", formula.trans, fixed = TRUE)
  formula.trans <- gsub(" ", "", formula.trans, fixed = TRUE)
  
  sigma.formula.trans <- gsub("(", "left_parenthesis", sigma.formulatext, fixed = TRUE)
  sigma.formula.trans <- gsub(")", "right_parenthesis", sigma.formula.trans, fixed = TRUE)
  sigma.formula.trans <- gsub("~", "tilde_symbol", sigma.formula.trans, fixed = TRUE)
  sigma.formula.trans <- gsub("=", "equal_symbol", sigma.formula.trans, fixed = TRUE)
  sigma.formula.trans <- gsub(",", "comma_symbol", sigma.formula.trans, fixed = TRUE)
  sigma.formula.trans <- gsub("*", "asterisk_symbol", sigma.formula.trans, fixed = TRUE)
  sigma.formula.trans <- gsub("^", "caret_symbol", sigma.formula.trans, fixed = TRUE)
  sigma.formula.trans <- gsub(" ", "", sigma.formula.trans, fixed = TRUE)
  
  nu.formula.trans <- gsub("(", "left_parenthesis", nu.formulatext, fixed = TRUE)
  nu.formula.trans <- gsub(")", "right_parenthesis", nu.formula.trans, fixed = TRUE)
  nu.formula.trans <- gsub("~", "tilde_symbol", nu.formula.trans, fixed = TRUE)
  nu.formula.trans <- gsub("=", "equal_symbol", nu.formula.trans, fixed = TRUE)
  nu.formula.trans <- gsub(",", "comma_symbol", nu.formula.trans, fixed = TRUE)
  nu.formula.trans <- gsub("*", "asterisk_symbol", nu.formula.trans, fixed = TRUE)
  nu.formula.trans <- gsub("^", "caret_symbol", nu.formula.trans, fixed = TRUE)
  nu.formula.trans <- gsub(" ", "", nu.formula.trans, fixed = TRUE)
  
  tau.formula.trans <- gsub("(", "left_parenthesis", tau.formulatext, fixed = TRUE)
  tau.formula.trans <- gsub(")", "right_parenthesis", tau.formula.trans, fixed = TRUE)
  tau.formula.trans <- gsub("~", "tilde_symbol", tau.formula.trans, fixed = TRUE)
  tau.formula.trans <- gsub("=", "equal_symbol", tau.formula.trans, fixed = TRUE)
  tau.formula.trans <- gsub(",", "comma_symbol", tau.formula.trans, fixed = TRUE)
  tau.formula.trans <- gsub("*", "asterisk_symbol", tau.formula.trans, fixed = TRUE)
  tau.formula.trans <- gsub("^", "caret_symbol", tau.formula.trans, fixed = TRUE)
  tau.formula.trans <- gsub(" ", "", tau.formula.trans, fixed = TRUE)
  
  family.trans <- gsub("(", "left_parenthesis", family, fixed = TRUE)
  family.trans <- gsub(")", "right_parenthesis", family.trans, fixed = TRUE)
  family.trans <- gsub("=", "equal_symbol", family.trans, fixed = TRUE)
  family.trans <- gsub(",", "comma_symbol", family.trans, fixed = TRUE)
  family.trans <- gsub(" ", "", family.trans, fixed = TRUE)
  familytext <- family
  family <- gamlss.dist::as.family(eval(parse(text=family), envir=environment()))
  
  if(is.null(mu.coef.start.names)){
    mu.coef.start.names.trans <- NULL
  } else {
    mu.coef.start.names.trans <- paste0(mu.coef.start.names, collapse=",")
    mu.coef.start.names.trans <- gsub("(", "left_parenthesis", mu.coef.start.names.trans, fixed = TRUE)
    mu.coef.start.names.trans <- gsub(")", "right_parenthesis", mu.coef.start.names.trans, fixed = TRUE)
    mu.coef.start.names.trans <- gsub("*", "asterisk_symbol", mu.coef.start.names.trans, fixed = TRUE)
    mu.coef.start.names.trans <- gsub("^", "caret_symbol", mu.coef.start.names.trans, fixed = TRUE)
    mu.coef.start.names.trans <- gsub(" ", "", mu.coef.start.names.trans, fixed = TRUE)
  }
  
  if(is.null(sigma.coef.start.names)){
    sigma.coef.start.names.trans <- NULL
  } else {
    sigma.coef.start.names.trans <- paste0(sigma.coef.start.names, collapse=",")
    sigma.coef.start.names.trans <- gsub("(", "left_parenthesis", sigma.coef.start.names.trans, fixed = TRUE)
    sigma.coef.start.names.trans <- gsub(")", "right_parenthesis", sigma.coef.start.names.trans, fixed = TRUE)
    sigma.coef.start.names.trans <- gsub("*", "asterisk_symbol", sigma.coef.start.names.trans, fixed = TRUE)
    sigma.coef.start.names.trans <- gsub("^", "caret_symbol", sigma.coef.start.names.trans, fixed = TRUE)
    sigma.coef.start.names.trans <- gsub(" ", "", sigma.coef.start.names.trans, fixed = TRUE)
  }

  if(is.null(nu.coef.start.names)){
    nu.coef.start.names.trans <- NULL
  } else {
    nu.coef.start.names.trans <- paste0(nu.coef.start.names, collapse=",")
    nu.coef.start.names.trans <- gsub("(", "left_parenthesis", nu.coef.start.names.trans, fixed = TRUE)
    nu.coef.start.names.trans <- gsub(")", "right_parenthesis", nu.coef.start.names.trans, fixed = TRUE)
    nu.coef.start.names.trans <- gsub("*", "asterisk_symbol", nu.coef.start.names.trans, fixed = TRUE)
    nu.coef.start.names.trans <- gsub("^", "caret_symbol", nu.coef.start.names.trans, fixed = TRUE)
    nu.coef.start.names.trans <- gsub(" ", "", nu.coef.start.names.trans, fixed = TRUE)
  }

  if(is.null(tau.coef.start.names)){
    tau.coef.start.names.trans <- NULL
  } else {
    tau.coef.start.names.trans <- paste0(tau.coef.start.names, collapse=",")
    tau.coef.start.names.trans <- gsub("(", "left_parenthesis", tau.coef.start.names.trans, fixed = TRUE)
    tau.coef.start.names.trans <- gsub(")", "right_parenthesis", tau.coef.start.names.trans, fixed = TRUE)
    tau.coef.start.names.trans <- gsub("*", "asterisk_symbol", tau.coef.start.names.trans, fixed = TRUE)
    tau.coef.start.names.trans <- gsub("^", "caret_symbol", tau.coef.start.names.trans, fixed = TRUE)
    tau.coef.start.names.trans <- gsub(" ", "", tau.coef.start.names.trans, fixed = TRUE)
  }

  # transform the control parameters into characters
  control.trans <- paste0(as.character(control), collapse=",")
  i.control.trans <- paste0(as.character(i.control), collapse=",")
  
  ## Vectors of beta values
  # arbitrary length for start betas at this stage but in legal transmission format ("0,0,0,0,0")
  # To DO: check whether this is really necessary
  mu.beta.vect <- rep(0,5)
  sigma.beta.vect <- rep(0,5)
  nu.beta.vect <- rep(0,5)
  tau.beta.vect <- rep(0,5)
  mu.beta.vect.trans <- paste0(as.character(mu.beta.vect), collapse=",")
  sigma.beta.vect.trans <- paste0(as.character(sigma.beta.vect), collapse=",")
  nu.beta.vect.trans <- paste0(as.character(nu.beta.vect), collapse=",")
  tau.beta.vect.trans <- paste0(as.character(tau.beta.vect), collapse=",")
  
  ## Get global mean & sd for outcome
  # (necessary for initialization of distribution parameters on the server side
  # for certain families)
  outcome <- strsplit(formulatext, "~", fixed=TRUE)[[1]][1]
  # global mean (required by most distributions)
  global.mean <- getPooledMean(datasources, outcome)
  if (familytext=="NO()"){
    # global sd
    # attention: leads slightly different results than sd() on pooled data
    global.sd <- sqrt(getPooledVar(datasources, outcome))
  } else {
    global.sd <- NULL
  }
  
  ## Identify the correct dimension for start betas via calling first component of gamlssDS
  ## Initialize the distribution parameter estimates on the server side
  cally1 <- call('gamlssDS1', formula=formula.trans, sigma.formula=sigma.formula.trans, 
                 nu.formula=nu.formula.trans, tau.formula=tau.formula.trans,
                 family=family.trans, data=data, mu.coef.start=mu.coef.start, 
                 sigma.coef.start=sigma.coef.start, nu.coef.start=nu.coef.start,
                 tau.coef.start=tau.coef.start, mu.coef.start.names=mu.coef.start.names.trans, 
                 sigma.coef.start.names=sigma.coef.start.names.trans, 
                 nu.coef.start.names=nu.coef.start.names.trans, 
                 tau.coef.start.names=tau.coef.start.names.trans, 
                 mu.fix=mu.fix, sigma.fix=sigma.fix, 
                 nu.fix=nu.fix, tau.fix=tau.fix, global.mean=global.mean, global.sd=global.sd,
                 control=control.trans, i.control=i.control.trans)
  
  study.summary.0 <- DSI::datashield.aggregate(datasources, cally1)
  
  #**************************************************************************
  # II) Disclosure Risk ----
  # summarizes disclosure risk from the individual servers
  # returns error messages & stops execution of function if necessary
  #**************************************************************************
  
  at.least.one.study.data.error <- 0
  
  for(hh in 1:numstudies) {
    if(study.summary.0[[hh]]$errorMessage!="No errors"){
      at.least.one.study.data.error <- 1
    }
  }
  
  mu.num.par <- NULL
  sigma.num.par <- NULL
  nu.num.par <- NULL
  tau.num.par <- NULL
  mu.num.par.gamma <- NULL
  sigma.num.par.gamma <- NULL
  nu.num.par.gamma <- NULL
  tau.num.par.gamma <- NULL
  
  mu.coef.names <- NULL
  sigma.coef.names <- NULL
  nu.coef.names <- NULL
  tau.coef.names <- NULL
  
  
  if(at.least.one.study.data.error==0){
    mod.gamlss.ds <- study.summary.0[[1]]$mod.gamlss.ds
    mu.num.par <- study.summary.0[[1]]$dim.mu.x[[2]]
    sigma.num.par <- study.summary.0[[1]]$dim.sigma.x[[2]]
    nu.num.par <- study.summary.0[[1]]$dim.nu.x[[2]]
    tau.num.par <- study.summary.0[[1]]$dim.tau.x[[2]]
    mu.coef.names <- names(mod.gamlss.ds$mu.coefficients)
    sigma.coef.names <- names(mod.gamlss.ds$sigma.coefficients)
    nu.coef.names <- names(mod.gamlss.ds$nu.coefficients)
    tau.coef.names <- names(mod.gamlss.ds$tau.coefficients)
    smoother.xmin <- study.summary.0[[1]]$smoother.xmin
    smoother.xmax <- study.summary.0[[1]]$smoother.xmax
    G.dev <- study.summary.0[[1]]$G.dev
    ## get the length of the gamma vectors
    # mu
    if (!is.null(mod.gamlss.ds$mu.coefSmo)){
      for (s in 1:length(mod.gamlss.ds$mu.coefSmo)){
        mu.num.par.gamma <- c(mu.num.par.gamma, dim(mod.gamlss.ds$mu.coefSmo[[s]]$coef)[1])
      }
    }
    # sigma
    if (!is.null(mod.gamlss.ds$sigma.coefSmo)){
      for (s in 1:length(mod.gamlss.ds$sigma.coefSmo)){
        sigma.num.par.gamma <- c(sigma.num.par.gamma, dim(mod.gamlss.ds$sigma.coefSmo[[s]]$coef)[1])
      }
    }
    # nu
    if (!is.null(mod.gamlss.ds$nu.coefSmo)){
      for (s in 1:length(mod.gamlss.ds$nu.coefSmo)){
        nu.num.par.gamma <- c(nu.num.par.gamma, dim(mod.gamlss.ds$nu.coefSmo[[s]]$coef)[1])
      }
    }
    # tau
    if (!is.null(mod.gamlss.ds$tau.coefSmo)){
      for (s in 1:length(mod.gamlss.ds$tau.coefSmo)){
        tau.num.par.gamma <- c(tau.num.par.gamma, dim(mod.gamlss.ds$tau.coefSmo[[s]]$coef)[1])
      }
    }
  }
  y.invalid <- NULL
  mu.par.invalid <- NULL
  sigma.par.invalid <- NULL
  nu.par.invalid <- NULL
  tau.par.invalid <- NULL
  gamlss.saturation.invalid <- NULL
  errorMessage <- NULL
  weights <- NULL
  N <- NULL
  noObs <- NULL
  mu.offset <- NULL
  sigma.offset <- NULL
  nu.offset <- NULL
  tau.offset <- NULL
  
  for(ss in 1:numstudies){
    # overall minimum and maximum (anonymized) for the smoother variables
    # Note that the minimum and maximum are needed to use the same knots on all servers
    smoother.xmin <- pmin(smoother.xmin, study.summary.0[[ss]]$smoother.xmin)  # pairwise minimum (for each variable)
    smoother.xmax <- pmax(smoother.xmax, study.summary.0[[ss]]$smoother.xmax)  # pairwise minimum (for each variable)
    
    # weights
    weights <- c(weights, study.summary.0[[ss]]$mod.gamlss.ds$weights)
    
    # deviance
    G.dev <- c(G.dev, study.summary.0[[ss]]$G.dev)
    
    # sample size
    N <- sum(N, study.summary.0[[ss]]$mod.gamlss.ds$N)
    noObs <- sum(noObs, study.summary.0[[ss]]$mod.gamlss.ds$noObs)
    
    # offset
    mu.offset <- c(mu.offset, study.summary.0[[ss]]$mod.gamlss.ds$mu.offset)
    sigma.offset <- c(sigma.offset, study.summary.0[[ss]]$mod.gamlss.ds$sigma.offset)
    nu.offset <- c(nu.offset, study.summary.0[[ss]]$mod.gamlss.ds$nu.offset)
    tau.offset <- c(tau.offset, study.summary.0[[ss]]$mod.gamlss.ds$tau.offset)
    
    # disclosure risk
    y.invalid <- c(y.invalid, study.summary.0[[ss]]$y.invalid)
    mu.par.invalid <- rbind(mu.par.invalid, study.summary.0[[ss]]$mu.par.invalid)
    sigma.par.invalid <- rbind(sigma.par.invalid, study.summary.0[[ss]]$sigma.par.invalid)
    nu.par.invalid <- rbind(nu.par.invalid, study.summary.0[[ss]]$nu.par.invalid)
    tau.par.invalid <- rbind(tau.par.invalid, study.summary.0[[ss]]$tau.par.invalid)
    gamlss.saturation.invalid <- c(gamlss.saturation.invalid, study.summary.0[[ss]]$gamlss.saturation.invalid)
    errorMessage <- c(errorMessage, study.summary.0[[ss]]$errorMessage)
  }
  
  # assign to gamlss model
  mod.gamlss.ds$weights <- weights
  if("mu" %in% names(family$parameters)){
    mod.gamlss.ds$mu.formula <- formula
    mod.gamlss.ds$mu.offset <- mu.offset
  }
  if("sigma" %in% names(family$parameters)){
    mod.gamlss.ds$sigma.formula <- sigma.formula
    mod.gamlss.ds$sigma.offset <- sigma.offset
  }
  if("nu" %in% names(family$parameters)){
    mod.gamlss.ds$nu.formula <- nu.formula
    mod.gamlss.ds$nu.offset <- nu.offset
  }
  if("tau" %in% names(family$parameters)){
    mod.gamlss.ds$tau.formula <- tau.formula
    mod.gamlss.ds$tau.offset <- tau.offset
  }
  # the dataname is added for later prediction
  mod.gamlss.ds$dataname <- data
  
  y.invalid <- as.matrix(y.invalid)
  sum.y.invalid <- sum(y.invalid)
  dimnames(y.invalid) <- list(names(datasources), "Y VECTOR")
  
  if(!is.null(mu.par.invalid)){
    mu.par.invalid <- as.matrix(mu.par.invalid)
    dimnames(mu.par.invalid) <- list(names(datasources), mu.coef.names)
  }
  sum.mu.par.invalid <- sum(mu.par.invalid)
  
  if(!is.null(mu.par.invalid)){
    sigma.par.invalid <- as.matrix(sigma.par.invalid)
    dimnames(sigma.par.invalid) <- list(names(datasources), sigma.coef.names)
  }
  sum.sigma.par.invalid <- sum(sigma.par.invalid)
  
  if(!is.null(nu.par.invalid)){
    nu.par.invalid <- as.matrix(nu.par.invalid)
    dimnames(nu.par.invalid) <- list(names(datasources), nu.coef.names)
  }
  sum.nu.par.invalid <- sum(nu.par.invalid)
  
  if(!is.null(tau.par.invalid)){
    tau.par.invalid <- as.matrix(tau.par.invalid)
    dimnames(tau.par.invalid) <- list(names(datasources), tau.coef.names)
  }
  sum.tau.par.invalid <- sum(tau.par.invalid)
  
  gamlss.saturation.invalid <- as.matrix(gamlss.saturation.invalid)
  sum.gamlss.saturation.invalid <- sum(gamlss.saturation.invalid)
  dimnames(gamlss.saturation.invalid) <- list(names(datasources), "MODEL OVERPARAMETERIZED")
  
  errorMessage <- as.matrix(errorMessage)
  dimnames(errorMessage) <- list(names(datasources), "ERROR MESSAGES")
  
  output.blocked.information.1 <- "MODEL FITTING TERMINATED AT FIRST ITERATION:"
  output.blocked.information.2 <- "Any values of 1 in the following tables denote potential disclosure risks"
  output.blocked.information.3 <- "please use the argument <datasources> to include only valid studies."
  output.blocked.information.4 <- "Errors by study are as follows:"
  
  if(sum.y.invalid>0||sum.mu.par.invalid>0||sum.sigma.par.invalid>0||sum.nu.par.invalid>0||sum.tau.par.invalid>0||
     sum.gamlss.saturation.invalid>0||at.least.one.study.data.error==1){
    message("\n\nMODEL FITTING TERMINATED AT FIRST ITERATION:\n",
            "Any values of 1 in the following tables denote potential disclosure risks\n",
            "please use the argument <datasources> to include only valid studies.\n",
            "Errors by study are as follows:\n")
    print(as.matrix(y.invalid))
    print(as.matrix(mu.par.invalid))
    print(as.matrix(sigma.par.invalid))
    print(as.matrix(nu.par.invalid))
    print(as.matrix(tau.par.invalid))
    print(as.matrix(gamlss.saturation.invalid))
    print(as.matrix(errorMessage))
    
    return(list(
      output.blocked.information.1,
      output.blocked.information.2,
      output.blocked.information.3,
      output.blocked.information.4,
      y.vector.error = y.invalid,
      mu.x.matrix.error = mu.par.invalid,
      sigma.x.matrix.error = sigma.par.invalid,
      nu.x.matrix.error = nu.par.invalid,
      tau.x.matrix.error = tau.par.invalid,
      gamlss.overparameterized = gamlss.saturation.invalid,
      errorMessage = errorMessage
    ))
    stop("DATA ERROR")
  }
  
  #**************************************************************************
  # III) Initialization ----
  #**************************************************************************
  
  ## Names of the parameters
  parameters <- mod.gamlss.ds$parameters
  smoother.names <- study.summary.0[[1]]$smoother.names
  
  ## Vectors of beta values with correct dimensions
  if(!is.null(mu.num.par)){
    mu.beta.vect <- rep(0, mu.num.par)
  }else{
    mu.beta.vect <- NULL
  }
  
  if(!is.null(sigma.num.par)){
    sigma.beta.vect <- rep(0, sigma.num.par)
  }else{
    sigma.beta.vect <- NULL
  }
  
  if(!is.null(nu.num.par)){
    nu.beta.vect <- rep(0, nu.num.par)
  }else{
    nu.beta.vect <- NULL
  }
  
  if(!is.null(tau.num.par)){
    tau.beta.vect <- rep(0, tau.num.par)
  }else{
    tau.beta.vect <- NULL
  }
  
  # vectors in legal transmission format ("0","0",...)
  mu.beta.vect.trans <- paste0(as.character(mu.beta.vect), collapse=",")
  sigma.beta.vect.trans <- paste0(as.character(sigma.beta.vect), collapse=",")
  nu.beta.vect.trans <- paste0(as.character(nu.beta.vect), collapse=",")
  tau.beta.vect.trans <- paste0(as.character(tau.beta.vect), collapse=",")
  
  ## Vectors of gamma values with correct dimensions
  if(!is.null(mu.num.par.gamma)){
    mu.gamma.vect <- rep(0, sum(mu.num.par.gamma))
  }else{
    mu.gamma.vect <- NULL
  }
  
  if(!is.null(sigma.num.par.gamma)){
    sigma.gamma.vect <- rep(0, sum(sigma.num.par.gamma))
  }else{
    sigma.gamma.vect <- NULL
  }
  
  if(!is.null(nu.num.par.gamma)){
    nu.gamma.vect <- rep(0, sum(nu.num.par.gamma))
  }else{
    nu.gamma.vect <- NULL
  }
  
  if(!is.null(tau.num.par.gamma)){
    tau.gamma.vect <- rep(0, sum(tau.num.par.gamma))
  }else{
    tau.gamma.vect <- NULL
  }
  
  # vectors in legal transmission format ("0","0",...)
  mu.gamma.vect.trans <- paste0(as.character(mu.gamma.vect), collapse=",")
  sigma.gamma.vect.trans <- paste0(as.character(sigma.gamma.vect), collapse=",")
  nu.gamma.vect.trans <- paste0(as.character(nu.gamma.vect), collapse=",")
  tau.gamma.vect.trans <- paste0(as.character(tau.gamma.vect), collapse=",")
  
  # Left & right boundary for the knots for the pb.smoother variables
  if (!is.null(min.max.names)){
    if (length(grep("$", min.max.names, fixed=TRUE))==0){
      min.max.names <- paste(data, min.max.names, sep="$")
    }
    position <- match(smoother.names, min.max.names)
  }
  if (!is.null(min.values)){
    smoother.xmin[position] <- min.values
  }
  if (!is.null(max.values)){
    smoother.xmax[position] <- max.values
  }
  smoother.xl <- smoother.xmin - 0.01 * (smoother.xmax - smoother.xmin)
  smoother.xr <- smoother.xmax + 0.01 * (smoother.xmax - smoother.xmin)
  
  # convert them to character vector
  smoother.xl.trans <- paste0(as.character(smoother.xl), collapse=",")
  smoother.xr.trans <- paste0(as.character(smoother.xr), collapse=",")
  
  ## Initialize deviance
  G.dev <- sum(G.dev)
  G.dev.old <- G.dev+1
  
  # Convergence state needs to be monitored.
  outer.converge.state <- FALSE
  
  # get the control parameters
  c.crit <- control[1]  # convergence criterion for outer iteration
  n.cyc <- control[2]  # maximum number of cycles for outer iteration
  gd.tol <- control[3]  # global deviance tolerance level
  mu.step <- control[4]  # step length for parameter mu
  sigma.step <- control[5]  # step length for parameter sigma
  nu.step <- control[6]  # step length for parameter nu
  tau.step <- control[7]  # step length for parameter tau
  
  cc <- i.control[1]  # convergence criterion for inner iteration
  cyc <- i.control[2]  # maximum number of cycles for inner iteration
  bf.cyc <- i.control[3]  # maximum number of cycles for backfitting
  bf.tol <- i.control[4]  # convergence criterion (tolerance level) for backfitting
  
  ## Control parameters for smoothers 
  mu.smoother.coef <- mu.coef.names[grep(pattern="pb(", x=tolower(mu.coef.names), fixed=TRUE)]
  sigma.smoother.coef <- sigma.coef.names[grep(pattern="pb(", x=tolower(sigma.coef.names), fixed=TRUE)]
  nu.smoother.coef <- nu.coef.names[grep(pattern="pb(", x=tolower(nu.coef.names), fixed=TRUE)]
  tau.smoother.coef <- tau.coef.names[grep(pattern="pb(", x=tolower(tau.coef.names), fixed=TRUE)]
  # check whether lambda is fixed
  mu.fixed.lambda <- grepl(pattern="lambda", x=mu.smoother.coef, fixed=TRUE)
  sigma.fixed.lambda <- grepl(pattern="lambda", x=sigma.smoother.coef, fixed=TRUE)
  nu.fixed.lambda <- grepl(pattern="lambda", x=nu.smoother.coef, fixed=TRUE)
  tau.fixed.lambda <- grepl(pattern="lambda", x=tau.smoother.coef, fixed=TRUE)
  # only keep the arguments for the pb() function
  mu.pb.args <- substr(mu.smoother.coef, start=4, stop=nchar(mu.smoother.coef)-1)
  mu.pb.args <- strsplit(mu.pb.args, split=",", fixed=TRUE)
  sigma.pb.args <- substr(sigma.smoother.coef, start=4, stop=nchar(sigma.smoother.coef)-1)
  sigma.pb.args <- strsplit(sigma.pb.args, split=",", fixed=TRUE)
  nu.pb.args <- substr(nu.smoother.coef, start=4, stop=nchar(nu.smoother.coef)-1)
  nu.pb.args <- strsplit(nu.pb.args, split=",", fixed=TRUE)
  tau.pb.args <- substr(tau.smoother.coef, start=4, stop=nchar(tau.smoother.coef)-1)
  tau.pb.args <- strsplit(tau.pb.args, split=",", fixed=TRUE)
  # get control parameters for the pb smoother
  mu.pb.control <- lapply(mu.pb.args, FUN=getpbcontrol)
  sigma.pb.control <- lapply(sigma.pb.args, FUN=getpbcontrol)
  nu.pb.control <- lapply(nu.pb.args, FUN=getpbcontrol)
  tau.pb.control <- lapply(tau.pb.args, FUN=getpbcontrol)
  # set starting value for lambda if lambda is not fixed & update knots for smoother
  if (length(mod.gamlss.ds$mu.coefSmo)>0){
    mod.gamlss.ds$mu.coefSmo[[which(!mu.fixed.lambda)]]$lambda <- mu.pb.control[[which(!mu.fixed.lambda)]]$start
    for (i in 1:length(mod.gamlss.ds$mu.coefSmo)){
      name <- mod.gamlss.ds$mu.coefSmo[[i]]$name
      xl <- smoother.xl[which(smoother.names==name)]
      xr <- smoother.xr[which(smoother.names==name)]
      dx <- (xr-xl)/mu.pb.control[[i]]$inter # increment to ensure the desired number of intervals ndx 
      deg <- mu.pb.control[[i]]$degree
      knots <- seq(xl-deg*dx, xr+deg*dx, by=dx)
      n <- length(knots)
      mod.gamlss.ds$mu.coefSmo[[i]]$knots <- knots[-c(1:(deg-1), (n-(deg-2)):n)]
    }
  }
  if (length(mod.gamlss.ds$sigma.coefSmo)>0){
    mod.gamlss.ds$sigma.coefSmo[[which(!sigma.fixed.lambda)]]$lambda <- sigma.pb.control[[which(!sigma.fixed.lambda)]]$start
    for (i in 1:length(mod.gamlss.ds$sigma.coefSmo)){
      name <- mod.gamlss.ds$sigma.coefSmo[[i]]$name
      xl <- smoother.xl[which(smoother.names==name)]
      xr <- smoother.xr[which(smoother.names==name)]
      dx <- (xr-xl)/sigma.pb.control[[i]]$inter # increment to ensure the desired number of intervals ndx 
      deg <- sigma.pb.control[[i]]$degree
      knots <- seq(xl-deg*dx, xr+deg*dx, by=dx)
      n <- length(knots)
      mod.gamlss.ds$sigma.coefSmo[[i]]$knots <- knots[-c(1:(deg-1), (n-(deg-2)):n)]
    }
  }
  if (length(mod.gamlss.ds$nu.coefSmo)>0){
    mod.gamlss.ds$nu.coefSmo[[which(!nu.fixed.lambda)]]$lambda <- nu.pb.control[[which(!nu.fixed.lambda)]]$start
    for (i in 1:length(mod.gamlss.ds$nu.coefSmo)){
      name <- mod.gamlss.ds$nu.coefSmo[[i]]$name
      xl <- smoother.xl[which(smoother.names==name)]
      xr <- smoother.xr[which(smoother.names==name)]
      dx <- (xr-xl)/nu.pb.control[[i]]$inter # increment to ensure the desired number of intervals ndx 
      deg <- nu.pb.control[[i]]$degree
      knots <- seq(xl-deg*dx, xr+deg*dx, by=dx)
      n <- length(knots)
      mod.gamlss.ds$nu.coefSmo[[i]]$knots <- knots[-c(1:(deg-1), (n-(deg-2)):n)]
    }
  }
  if (length(mod.gamlss.ds$tau.coefSmo)>0){
    mod.gamlss.ds$tau.coefSmo[[which(!tau.fixed.lambda)]]$lambda <- tau.pb.control[[which(!tau.fixed.lambda)]]$start
    for (i in 1:length(mod.gamlss.ds$tau.coefSmo)){
      name <- mod.gamlss.ds$tau.coefSmo[[i]]$name
      xl <- smoother.xl[which(smoother.names==name)]
      xr <- smoother.xr[which(smoother.names==name)]
      dx <- (xr-xl)/tau.pb.control[[i]]$inter # increment to ensure the desired number of intervals ndx 
      deg <- tau.pb.control[[i]]$degree
      knots <- seq(xl-deg*dx, xr+deg*dx, by=dx)
      n <- length(knots)
      mod.gamlss.ds$tau.coefSmo[[i]]$knots <- knots[-c(1:(deg-1), (n-(deg-2)):n)]
    }
  }
  
  #**************************************************************************
  # IV) Iteration ----
  #**************************************************************************
  
  #*A) Outer iteration ----
  while(abs(G.dev.old-G.dev) > c.crit && outer.iteration.count < n.cyc){
    outer.iteration.count <- outer.iteration.count+1
    message("Outer iteration ", outer.iteration.count, "...")
    
    for (p in 1:length(parameters)){
      parameter <- parameters[p]
      fixed.lambda <- eval(parse(text=paste(parameter, ".fixed.lambda", sep="")), envir=environment())
      pb.control <- eval(parse(text=paste(parameter, ".pb.control", sep="")), envir=environment())
      num.par.gamma <- eval(parse(text=paste(parameter, ".num.par.gamma", sep="")), envir=environment())
      gamma.vect <- eval(parse(text=paste(parameter, ".gamma.vect", sep="")), envir=environment())
      
      #*A.1) Inner iteration ----
      inner.converge.state <- FALSE
      inner.iteration.count <- 0
      inner.dev.increase <- FALSE  # indicator whether deviance has increased in an inner iteration
      
      # for the first iteration the dv & olddv are set to arbitrary values with (abs(olddv-dv)) < cc to ensure that 
      # the inner iteration starts
      # to ensure that the correct deviance is used for future iterations the values are corrected later
      dv.total <- 100000
      olddv.total <- dv.total+1
      
      while (abs(olddv.total - dv.total) > cc && inner.iteration.count < cyc){
        inner.iteration.count <- inner.iteration.count+1
        
        bf.iteration.count <- 0
        
        # for the first backfitting iteration when the ratio does not exist to ensure that
        # the algorithm continues (ratio > tol)
        ratio <- bf.tol + 1
        
        #*A.1.i) Backfitting ----
        while (ratio > bf.tol & bf.iteration.count < bf.cyc){
          bf.iteration.count <- bf.iteration.count+1
          
          #*A.1.i.a) WLS ----
          ## call second component of gamlssDS to generate matrices and vectors for WLS to estimate beta
          cally2 <- call('gamlssDS2', parameter=parameter, formula=formula.trans, sigma.formula=sigma.formula.trans, 
                         nu.formula=nu.formula.trans, tau.formula=tau.formula.trans, family=family.trans, data=data, 
                         mu.fix=mu.fix, sigma.fix=sigma.fix, nu.fix=nu.fix, tau.fix=tau.fix,
                         mu.beta.vect=mu.beta.vect.trans, sigma.beta.vect=sigma.beta.vect.trans,
                         nu.beta.vect=nu.beta.vect.trans, tau.beta.vect=tau.beta.vect.trans,
                         control=control.trans, i.control=i.control.trans)
          study.summary <- DSI::datashield.aggregate(datasources, cally2)
          
          ## combine the aggregated results from all servers to obtain beta estimates
          disclosure.risk.total <- Reduce(f="+", .select(study.summary, 'disclosure.risk'))
          disclosure.risk <- NULL
          errorMessage2 <- NULL
          
          for(ss2 in 1:numstudies){
            disclosure.risk <- c(disclosure.risk, study.summary[[ss]]$disclosure.risk)
            errorMessage2 <- c(errorMessage2, study.summary[[ss]]$errorMessage2)
          }
          
          disclosure.risk <- as.matrix(disclosure.risk)
          dimnames(disclosure.risk) <- list(names(datasources), "Risk of disclosure")
          
          errorMessage2 <- as.matrix(errorMessage2)
          dimnames(errorMessage2) <- list(names(datasources), "Error messages")
          
          if(disclosure.risk.total>0){
            message("Potential disclosure risk in y, mu.x, sigma.x, nu.x or tau.x \n",
                    "or model overparameterized in at least one study.\n",
                    "In addition clientside function appears to have been modified \n",
                    "to avoid traps in first serverside function.\n",
                    "Vectors and matrices therefore destroyed in all invalid studies \n",
                    "and model fitting terminated. This error is recorded in the log file but \n",
                    "please report it to the DataSHIELD team as we need to understand how\n",
                    "the controlled shutdown traps in gamlssDS1 have been circumvented\n\n")
            
            output.blocked.information.1 <- "Potential disclosure risk in y, mu.x, sigma.x, nu.x or tau.x"
            output.blocked.information.2 <- "or model overparameterized in at least one study."
            output.blocked.information.3 <- "In addition clientside function appears to have been modified"
            output.blocked.information.4 <- "to avoid disclosure traps in first serverside function."
            output.blocked.information.5 <- "Score vectors and information matrices therefore destroyed in all invalid studies"
            output.blocked.information.6 <- "and model fitting terminated. This error is recorded in the log file but"
            output.blocked.information.7 <- "please also report it to the DataSHIELD team as we need to understand how"
            output.blocked.information.8 <- "the controlled shutdown traps in glmDS1 have been circumvented."
            
            return(list(output.blocked.information.1,
                        output.blocked.information.2,
                        output.blocked.information.3,
                        output.blocked.information.4,
                        output.blocked.information.5,
                        output.blocked.information.6,
                        output.blocked.information.7,
                        output.blocked.information.8))
          }
          
          matrix.total <- Reduce(f="+", .select(study.summary, 'matrix'))
          vector.total <- Reduce(f="+", .select(study.summary, 'vector'))
          
          ## get the initial deviance for the inner iteration
          if (bf.iteration.count == 1){
            dv.total <- Reduce(f="+", .select(study.summary, 'dv'))
          }
          
          #*A.1.i.a.i) Calculate beta ----
          inverse.matrix.total <- solve(matrix.total)
          beta.vect.next <- as.vector(inverse.matrix.total %*% vector.total)
          # convert it to legal transmission format
          beta.vect.trans <- paste0(as.character(beta.vect.next), collapse=",")
          # assign the new beta vectors to the current parameter
          eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefficients <- beta.vect.next", sep="")), envir=environment())
          base::assign(paste(parameter, ".beta.vect.trans", sep=""), beta.vect.trans, envir=environment())
          
          ## Reset ratio 
          # to ensure that inner iteration works without backfitting if no smoothers are specifed
          # in the model
          ratio <- bf.tol  # this means that the while loop for the backfitting stops
          
          #*A.1.i.b) PWLS ----
          # if the parameter includes smoothers call the third component of gamlssDS to generate matrices and vectors
          # for PWLS to estimate gamma
          
          deltaf <- 0  # initialize change in function
          
          # loop over all smoothers
          smoothers <- NULL
          if (length(num.par.gamma)>0){
            smoothers <- seq(from=1, to=length(num.par.gamma))
          }
          
          for (s in smoothers){
            # the difference matrix of order
            if(pb.control[[s]]$order==0){
              D.mat <- diag(pb.control[[s]]$inter+pb.control[[s]]$degree)
            } else {
              D.mat <- diff(diag(pb.control[[s]]$inter+pb.control[[s]]$degree), diff=pb.control[[s]]$order)
            }
            
            ## call third component of gamlssDS to generate matrices and vectors for PWLS to estimate gamma
            cally3 <- call('gamlssDS3', parameter=parameter, smoother=s, formula=formula.trans, sigma.formula=sigma.formula.trans, 
                           nu.formula=nu.formula.trans, tau.formula=tau.formula.trans, family=family.trans, data=data, 
                           mu.fix=mu.fix, sigma.fix=sigma.fix, nu.fix=nu.fix, tau.fix=tau.fix,
                           mu.beta.vect=mu.beta.vect.trans, sigma.beta.vect=sigma.beta.vect.trans,
                           nu.beta.vect=nu.beta.vect.trans, tau.beta.vect=tau.beta.vect.trans,
                           mu.gamma.vect=mu.gamma.vect.trans, sigma.gamma.vect=sigma.gamma.vect.trans,
                           nu.gamma.vect=nu.gamma.vect.trans, tau.gamma.vect=tau.gamma.vect.trans,
                           smoother.names=smoother.names, smoother.xl=smoother.xl.trans, smoother.xr=smoother.xr.trans,
                           control=control.trans, i.control=i.control.trans)
            study.summary <- DSI::datashield.aggregate(datasources, cally3)
            
            ## combine the aggregated results from all servers to obtain gamma estimates
            matrix.total <- Reduce(f="+", .select(study.summary, 'matrix'))
            vector.total <- Reduce(f="+", .select(study.summary, 'vector'))
            sumofsquares.total <- Reduce(f="+", .select(study.summary, 'sumofsquares'))
            sumofweights.total <- Reduce(f="+", .select(study.summary, 'sumofweights'))
            
            ## update change in smoothing fitted values (as part of stopping criterion for backfitting)
            deltaf <- deltaf + (sumofsquares.total/sumofweights.total)  # change in smoothing fitted values
            
            #*A.1.i.b.i) Calculate gamma ----
            # figure out which part of the gamma vector should be changed
            if (s==1){
              gamma.start <- 1
            } else{
              gamma.start <- sum(num.par.gamma[1:(s-1)])+1
            }
            gamma.end <- gamma.start+num.par.gamma[s]-1
            
            # update gamma
            if (fixed.lambda[s]==TRUE){
              #*A.1.i.b.i.a) case 1: lambda is known ----
              lambda <- eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", s, "]]$lambda", sep="")), envir=environment())
              G.mat <- t(D.mat) %*% D.mat
              inverse.matrix.total <- solve(matrix.total + lambda * G.mat)
              gamma.vect.update <- as.vector(inverse.matrix.total %*% vector.total)
              gamma.vect[gamma.start:gamma.end] <- gamma.vect.update
              # convert it to legal transmission format
              gamma.vect.trans <- paste0(as.character(gamma.vect), collapse=",")
              # assign the new gamma vectors to the current parameter (use eval since assign does not work on attributes)
              eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[",s, "]]$coef <- gamma.vect.update", sep="")), envir=environment())
              base::assign(paste(parameter, ".gamma.vect.trans", sep=""), gamma.vect.trans, env=parent.frame())
              
              # update edf
              edf <- sum(diag(inverse.matrix.total %*% matrix.total))
              eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", s, "]]$edf <- edf", sep="")), envir=environment())
              
            } else {
              #*A.1.i.b.i.b) case 2: lambda is estimated ----
              if (pb.control[[s]]$method=="ML"){
                for (lambda.it in 1:50){
                  lambda <- eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", s, "]]$lambda", sep="")), envir=environment())
                  # estimate gamma (for current lambda)
                  G.mat <- t(D.mat) %*% D.mat
                  inverse.matrix.total <- solve(matrix.total + lambda*G.mat)
                  gamma.vect.update <- as.vector(inverse.matrix.total %*% vector.total)
                  gamma.vect[gamma.start:gamma.end] <- gamma.vect.update
                  # convert it to legal transmission format
                  gamma.vect.trans <- paste0(as.character(gamma.vect), collapse=",")
                  # assign the new gamma vectors to the current parameter
                  eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[",s, "]]$coef <- gamma.vect.update", sep="")), envir=environment())
                  base::assign(paste(parameter, ".gamma.vect.trans", sep=""), gamma.vect.trans, envir=environment())
                  
                  # estimate lambda
                  # call fourth component of gamlssDS to generate matrices and vectors for PWLS to estimate gamma
                  cally4 <- call('gamlssDS4', parameter=parameter, smoother=s, formula=formula.trans, sigma.formula=sigma.formula.trans, 
                                 nu.formula=nu.formula.trans, tau.formula=tau.formula.trans, family=family.trans, data=data, 
                                 mu.fix=mu.fix, sigma.fix=sigma.fix, nu.fix=nu.fix, tau.fix=tau.fix,
                                 mu.beta.vect=mu.beta.vect.trans, sigma.beta.vect=sigma.beta.vect.trans,
                                 nu.beta.vect=nu.beta.vect.trans, tau.beta.vect=tau.beta.vect.trans,
                                 mu.gamma.vect=mu.gamma.vect.trans, sigma.gamma.vect=sigma.gamma.vect.trans,
                                 nu.gamma.vect=nu.gamma.vect.trans, tau.gamma.vect=tau.gamma.vect.trans,
                                 smoother.names=smoother.names,
                                 smoother.xl=smoother.xl.trans, smoother.xr=smoother.xr.trans,
                                 control=control.trans, i.control=i.control.trans)
                  study.summary <- DSI::datashield.aggregate(datasources, cally4)
                  
                  ## combine the aggregated results from all servers to obtain estimates for lambda
                  nobs.total <- Reduce(f="+", .select(study.summary, 'nobs'))
                  product.total <- Reduce(f="+", .select(study.summary, 'inner.product'))
                  # obtain the appropriate degrees of freedom (according to formula in Rigby et al. 2013,
                  # approximation in appendix 2)
                  gamma.diff <- D.mat %*% as.vector(gamma.vect.update)
                  edf <- sum(diag(inverse.matrix.total %*% matrix.total))
                  sig2 <-  product.total / (nobs.total-edf)  
                  tau2 <- sum(gamma.diff^2) / (edf-pb.control[[s]]$order)
                  if(tau2 < 1e-7){
                    tau2 <- 1.0e-7
                  } 
                  lambda.old <- lambda
                  lambda <- sig2/tau2
                  if (lambda < 1.0e-7){
                    lambda <- 1.0e-7
                  } 
                  if (lambda > 1.0e+7){
                    lambda <- 1.0e+7
                  }
                  # assign the values to the gamlss object
                  eval(parse(text=paste("mod.gamlss.ds$", parameter, ".lambda[", s, "] <- lambda", sep="")), envir=environment())
                  eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", s, "]]$lambda <- lambda", sep="")), envir=environment())
                  eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", s, "]]$edf <- edf", sep="")), envir=environment())
                  eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", s, "]]$sigb2 <- tau2", sep="")), envir=environment())
                  eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", s, "]]$sigb <- sqrt(tau2)", sep="")), envir=environment())
                  eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", s, "]]$sige2 <- sig2", sep="")), envir=environment())
                  eval(parse(text=paste("mod.gamlss.ds$", parameter, ".coefSmo[[", s, "]]$sige <- sqrt(sig2)", sep="")), envir=environment())
                  
                  if (abs(lambda-lambda.old) < 1.0e-7 || lambda > 1.0e10){
                    # change in lambda
                    break
                  } 
                  # update lambda
                  # pb.lambda.start[lambda.start.pos+s] <- lambda
                } # end of for loop for ML method
              } 
              
            } # End of case 2: lambda is estimated
            
            
          } # End for loop over the smoothers
          
          #*A.1.i.c) Stopping criterion backfitting ----
          if(length(num.par.gamma)==0){
            # Ensure backfitting stops if no smoothing terms for the parameter 
            ratio <- bf.tol
          } else{
            
            ## call fifth component of gamlssDS to check convergence of backfitting
            cally5 <- call('gamlssDS5', parameter=parameter, formula=formula.trans, sigma.formula=sigma.formula.trans, 
                           nu.formula=nu.formula.trans, tau.formula=tau.formula.trans, family=family.trans, data=data,
                           mu.fix=mu.fix, sigma.fix=sigma.fix, nu.fix=nu.fix, tau.fix=tau.fix,
                           mu.beta.vect=mu.beta.vect.trans, sigma.beta.vect=sigma.beta.vect.trans,
                           nu.beta.vect=nu.beta.vect.trans, tau.beta.vect=tau.beta.vect.trans,
                           mu.gamma.vect=mu.gamma.vect.trans, sigma.gamma.vect=sigma.gamma.vect.trans,
                           nu.gamma.vect=nu.gamma.vect.trans, tau.gamma.vect=tau.gamma.vect.trans,
                           smoother.names=smoother.names,
                           smoother.xl=smoother.xl.trans, smoother.xr=smoother.xr.trans,
                           control=control.trans, i.control=i.control.trans)
            study.summary <- DSI::datashield.aggregate(datasources, cally5)
            
            ## combine the aggregated results from all servers to obtain gamma estimates
            sumofsquares.total <- Reduce(f="+", .select(study.summary, 'sumofsquares'))
            sumofweights.total <- Reduce(f="+", .select(study.summary, 'sumofweights'))
            sumofsmoothers.total <- Reduce(f="+", .select(study.summary, 'sumofsmoothers'))
            
            ## update change in smoothing fitted values (as part of stopping criterion for backfitting)
            deltaf <- deltaf + (sumofsquares.total/sumofweights.total)  # change in smoothing fitted values
            ratio <- sqrt(deltaf/sumofsmoothers.total)
          }
          
        } # End of backfitting
        
        #*A.1.ii) Stopping criterion inner iteration ----
        cally6 <- call('gamlssDS6', parameter=parameter, formula=formula.trans, sigma.formula=sigma.formula.trans, 
                       nu.formula=nu.formula.trans, tau.formula=tau.formula.trans, family=family.trans, data=data, 
                       mu.fix=mu.fix, sigma.fix=sigma.fix, nu.fix=nu.fix, tau.fix=tau.fix,
                       mu.beta.vect=mu.beta.vect.trans, sigma.beta.vect=sigma.beta.vect.trans,
                       nu.beta.vect=nu.beta.vect.trans, tau.beta.vect=tau.beta.vect.trans,
                       mu.gamma.vect=mu.gamma.vect.trans, sigma.gamma.vect=sigma.gamma.vect.trans,
                       nu.gamma.vect=nu.gamma.vect.trans, tau.gamma.vect=tau.gamma.vect.trans,
                       smoother.names=smoother.names,
                       smoother.xl=smoother.xl.trans, smoother.xr=smoother.xr.trans,
                       control=control.trans, i.control=i.control.trans, autostep=autostep,
                       inner.iteration.count=inner.iteration.count, autostep.count=0)
        study.summary <- DSI::datashield.aggregate(datasources, cally6)
        
        olddv.total <- dv.total
        dv.total <- Reduce(f="+", .select(study.summary, 'dv'))
        
        ## autostep
        # to avoid overjumping (method 2 as described in Stasinopolous et al. 2020, p.66f)
        if (autostep==TRUE){
          
          if (dv.total > olddv.total && inner.iteration.count >= 2){
            # autostep iteration (if deviance increased)
            for(autostep.count in 1:5){
              cally6 <- call('gamlssDS6', parameter=parameter, formula=formula.trans, sigma.formula=sigma.formula.trans, 
                             nu.formula=nu.formula.trans, tau.formula=tau.formula.trans, family=family.trans, data=data, 
                             mu.fix=mu.fix, sigma.fix=sigma.fix, nu.fix=nu.fix, tau.fix=tau.fix,
                             mu.beta.vect=mu.beta.vect.trans, sigma.beta.vect=sigma.beta.vect.trans,
                             nu.beta.vect=nu.beta.vect.trans, tau.beta.vect=tau.beta.vect.trans,
                             mu.gamma.vect=mu.gamma.vect.trans, sigma.gamma.vect=sigma.gamma.vect.trans,
                             nu.gamma.vect=nu.gamma.vect.trans, tau.gamma.vect=tau.gamma.vect.trans,
                             smoother.names=smoother.names,
                             smoother.xl=smoother.xl.trans, smoother.xr=smoother.xr.trans,
                             control=control.trans, i.control=i.control.trans, autostep=autostep,
                             inner.iteration.count=inner.iteration.count, autostep.count=autostep.count)
              study.summary <- DSI::datashield.aggregate(datasources, cally6)
              dv.total <- Reduce(f="+", .select(study.summary, 'dv'))
              if ((olddv.total-dv.total) > cc){
                break # new deviance smaller than old one
              }
            }
          } else {
            # only save parameter estimates (for 1st inner iteration or if deviance did not increase)
            cally6 <- call('gamlssDS6', parameter=parameter, formula=formula.trans, sigma.formula=sigma.formula.trans, 
                           nu.formula=nu.formula.trans, tau.formula=tau.formula.trans, family=family.trans, data=data,
                           mu.fix=mu.fix, sigma.fix=sigma.fix, nu.fix=nu.fix, tau.fix=tau.fix,
                           mu.beta.vect=mu.beta.vect.trans, sigma.beta.vect=sigma.beta.vect.trans,
                           nu.beta.vect=nu.beta.vect.trans, tau.beta.vect=tau.beta.vect.trans,
                           mu.gamma.vect=mu.gamma.vect.trans, sigma.gamma.vect=sigma.gamma.vect.trans,
                           nu.gamma.vect=nu.gamma.vect.trans, tau.gamma.vect=tau.gamma.vect.trans,
                           smoother.names=smoother.names,
                           smoother.xl=smoother.xl.trans, smoother.xr=smoother.xr.trans,
                           control=control.trans, i.control=i.control.trans, autostep=FALSE,
                           inner.iteration.count=inner.iteration.count, autostep.count=0)
            study.summary <- DSI::datashield.aggregate(datasources, cally6)
          }
          
        } # end autostep=TRUE
        
        ## check the deviance
        if ((dv.total > olddv.total + gd.tol ) && inner.iteration.count >= 2 && inner.dev.increase==FALSE){
          warning("The deviance has increased in an inner iteration for ",
                  parameter, "\n","Increase gd.tol and if persist, try different steps",  "\n", "or model maybe inappropriate")
          inner.dev.increase <- TRUE
        }
        
        ## save penalty & check whether updated parameters valid
        errorMessage3 <- NULL
        pen <- 0
        for(ss in 1:numstudies){
          errorMessage3 <- c(errorMessage3, study.summary[[ss]]$errorMessage3)
          pen <- sum(pen, study.summary[[ss]]$pen)
        }
        errorMessage3 <- unique(errorMessage3)
        # assign the penalty to the gamlss object if the parameter includes smoothers
        if (length(smoothers)>0){
          eval(parse(text=paste("mod.gamlss.ds$", parameter, ".pen <- pen", sep="")), envir=environment())
        }
        
        if (!is.null(errorMessage3)){
          stop(paste(errorMessage3))
        } 
        
      } # End of inner iteration
      
    } # End of for loop over the distribution parameters
    
    #*A.2) Stopping criterion outer iteration ----
    G.dev.old <- G.dev
    G.dev <- dv.total  # this should be weighted if weights are considered
    
    if (G.dev > (G.dev.old+gd.tol) && outer.iteration.count>1){
      stop(paste("The global deviance is increasing", "\n", 
                 "Try different steps for the parameters or the model maybe inappropriate"))
    }
    
  } # End of outer iteration
  
  
  #*B) Check convergence outer iteration ----
  # check whether the outer iteration algorithm converged
  if (abs(G.dev.old-G.dev) < c.crit){ 
    outer.converge.state <- TRUE
  } else {
    outer.converge.state <- FALSE
  }
  if (!outer.converge.state){
    warning("Algorithm RS has not yet converged")
  }
  
  #**************************************************************************
  # V) Output ----
  # Update GAMLSS Object
  #**************************************************************************
  mod.gamlss.ds$call <- "Call to ds.gamlss"
  
  mod.gamlss.ds$G.deviance <- G.dev
  mod.gamlss.ds$N <- N
  mod.gamlss.ds$iter <- outer.iteration.count
  mod.gamlss.ds$converged <- outer.converge.state
  mod.gamlss.ds$noObs <- noObs
  
  if("mu" %in% parameters){
    names(mod.gamlss.ds$mu.coefficients) <- mu.coef.names
    if(length(mod.gamlss.ds$mu.coefSmo)>0){
      mod.gamlss.ds$mu.nl.df <- 0
      pb.control <- lapply(mu.pb.args, FUN=getpbcontrol)
      for (i in 1:length(mod.gamlss.ds$mu.coefSmo)){
        mod.gamlss.ds$mu.nl.df <- sum(mod.gamlss.ds$mu.nl.df, mod.gamlss.ds$mu.coefSmo[[i]]$edf-2)
      }
      mod.gamlss.ds$mu.df <- mod.gamlss.ds$mu.nl.df+length(mod.gamlss.ds$mu.coefficients)
    }
  }
  if("sigma" %in% parameters){
    names(mod.gamlss.ds$sigma.coefficients) <- sigma.coef.names
    if(length(mod.gamlss.ds$sigma.coefSmo)>0){
      mod.gamlss.ds$sigma.nl.df <- 0
      for (i in 1:length(mod.gamlss.ds$sigma.coefSmo)){
        mod.gamlss.ds$sigma.nl.df <- sum(mod.gamlss.ds$sigma.nl.df, mod.gamlss.ds$sigma.coefSmo[[i]]$edf-2)
      }
      mod.gamlss.ds$sigma.df <- mod.gamlss.ds$sigma.nl.df+length(mod.gamlss.ds$sigma.coefficients)
    }
  }
  if("nu" %in% parameters){
    names(mod.gamlss.ds$nu.coefficients) <- nu.coef.names
    if(length(mod.gamlss.ds$nu.coefSmo)>0){
      mod.gamlss.ds$nu.nl.df <- 0
      for (i in 1:length(mod.gamlss.ds$nu.coefSmo)){
        mod.gamlss.ds$nu.nl.df <- sum(mod.gamlss.ds$nu.nl.df, mod.gamlss.ds$nu.coefSmo[[i]]$edf-2)
      }
      mod.gamlss.ds$nu.df <- mod.gamlss.ds$nu.nl.df+length(mod.gamlss.ds$nu.coefficients)
    }
  }
  if("tau" %in% parameters){
    names(mod.gamlss.ds$tau.coefficients) <- tau.coef.names
    if(length(mod.gamlss.ds$tau.coefSmo)>0){
      mod.gamlss.ds$tau.nl.df <- 0
      for (i in 1:length(mod.gamlss.ds$tau.coefSmo)){
        mod.gamlss.ds$tau.nl.df <- sum(mod.gamlss.ds$tau.nl.df, mod.gamlss.ds$tau.coefSmo[[i]]$edf-2)
      }
      mod.gamlss.ds$tau.df <- mod.gamlss.ds$tau.nl.df+length(mod.gamlss.ds$tau.coefficients)
    }
  }
  mod.gamlss.ds$df.fit <- sum(mod.gamlss.ds$mu.df, mod.gamlss.ds$sigma.df, mod.gamlss.ds$nu.df, mod.gamlss.ds$tau.df)
  mod.gamlss.ds$df.residual <- noObs - mod.gamlss.ds$df.fit
  mod.gamlss.ds$pen <- sum(mod.gamlss.ds$mu.pen, mod.gamlss.ds$sigma.pen, mod.gamlss.ds$nu.pen, mod.gamlss.ds$tau.pen)
  mod.gamlss.ds$P.deviance <- G.dev + mod.gamlss.ds$pen
  mod.gamlss.ds$aic <- G.dev + 2*mod.gamlss.ds$df.fit
  mod.gamlss.ds$sbc <- G.dev + log(noObs)*mod.gamlss.ds$df.fit
  
  # Residuals
  cally7 <- call('gamlssDS7', formula=formula.trans, sigma.formula=sigma.formula.trans, 
                 nu.formula=nu.formula.trans, tau.formula=tau.formula.trans,
                 family=family.trans, data=data, mu.fix=mu.fix, sigma.fix=sigma.fix, 
                 nu.fix=nu.fix, tau.fix=tau.fix,
                 control=control.trans, i.control=i.control.trans, k=k)
  study.summary.0 <- DSI::datashield.aggregate(datasources, cally7)
  residuals <- NULL
  
  for(s in 1:numstudies){
    residuals <- c(residuals, study.summary.0[[s]])
  }
  mod.gamlss.ds$residuals <- residuals
  
  return(mod.gamlss.ds)
  
}