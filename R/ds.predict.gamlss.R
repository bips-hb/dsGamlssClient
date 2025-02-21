#'
#' @title Derive predictions for \code{ds.gamlss} objects
#' @description This function predicts the specified distribution parameter for a new data set from a
#' \code{ds.gamlss} object that is output by the \code{\link[dsGamlssClient]{ds.gamlss}} function.
#' @details The \code{ds.predict.gamlss} function assumes that the object given in \code{newdata} is a
#' data frame containing the explanatory variables that are used in the model.
#' @param object A \code{ds.gamlss} object output by \code{\link[dsGamlssClient]{ds.gamlss}}.
#' @param newdata A data frame containing new values for the explanatory variables that are used in the model.
#' @param what A string, specifying the distribution parameter that should be predicted. Default \code{what="mu"}.
#' @param type A string, specifying the kind of predictions that should be derived. The default value \code{type="link"}
#' predicts the linear predictor for the specified distribution parameter, whereas \code{type="response"} predicts the
#' fitted values for the distribution parameter.
#' @return A vector with the predictions for \code{newdata}.
#' @author Annika Swenne
#' @export
#' @examples
#' library(DSLite)
#' data(mtcars)
#'
#' ## Create newdata for predictions
#' newdata <- data.frame(wt = seq(2, 5, by = 0.01))
#'
#' ## Fit ds.gamlss model
#' # Set up DSLite server
#' dslite.server1 <- newDSLiteServer(
#'   tables = list(data = mtcars[c(1:15), ]),
#'   config = defaultDSConfiguration(include = c("dsBase", "dsGamlss", "gamlss", "gamlss.dist"))
#' )
#' dslite.server2 <- newDSLiteServer(
#'   tables = list(data = mtcars[c(16:nrow(mtcars)), ]),
#'   config = defaultDSConfiguration(include = c("dsBase", "dsGamlss", "gamlss", "gamlss.dist"))
#' )
#' builder <- DSI::newDSLoginBuilder()
#' builder$append(server = "study1", url = "dslite.server1", table = "data", driver = "DSLiteDriver")
#' builder$append(server = "study2", url = "dslite.server2", table = "data", driver = "DSLiteDriver")
#' logindata.dslite <- builder$build()
#' # Login to the virtualized server
#' conns <- DSI::datashield.login(logindata.dslite, assign = TRUE)
#' DSI::datashield.assign.table(conns = conns, symbol = "D", table = c("data", "data"))
#' # Fit model
#' model <- ds.gamlss(formula = mpg ~ pb(wt), sigma.formula = ~wt, data = "D", family = "NO()")
#' # Logout
#' DSI::datashield.logout(conns)
#'
#' ## Examples
#' # Example 1: Predict mu
#' mu.response <- ds.predict.gamlss(model, newdata, what = "mu", type = "response")
#'
#' # Example 2: Predict linear predictor for sigma
#' sigma.link <- ds.predict.gamlss(model, newdata, what = "sigma", type = "link")
ds.predict.gamlss <- function(object, newdata, what = "mu",
                              type = "link") {
  create.bbase <- function(xvalues, knots, ndx = 20, deg = 3) {
    # Function to create basis for p-splines for new data
    # ***************************************************************************
    # xvalues : vector of x-values
    # knots : knots for the basis matrix (as included in the ds.gamlss object)
    # ndx : number of equal space intervals in x
    # deg : degree of the polynomial

    tpower <- function(x, t, p) {
      # Truncated p-th power function (defined for single values)
      # used to construct B-splines
      # *************************************************************************
      # x : single x-value
      # t : single knot
      # p : power
      (x - t)^p * (x > t) # x>t indicator function TRUE=1, FALSE=0
    } # end tpower function

    xl <- knots[2]
    xr <- knots[length(knots) - 1]
    dx <- (xr - xl) / ndx

    # Add the additional knots at left and right border
    knots.ext <- seq(xl - deg * dx, xr + deg * dx, by = dx)

    # calculate the power in the knots (evaluate function tpower for all combinations of x-values and knots)
    P <- outer(xvalues, knots.ext, tpower, deg) # matrix of dimension length(x)*length(knots)
    n <- dim(P)[2]
    D <- diff(diag(n), diff = deg + 1) / (gamma(deg + 1) * dx^deg)
    # diff(diag(n), diff=deg+1) difference matrix of order deg+1
    # gamma(deg+1) returns gamma function of deg+1 (single value)
    B <- (-1)^(deg + 1) * P %*% t(D)

    return(B)
  } # end create.bbase function

  ## Ensure the input parameters are properly specified
  # object
  if (!(gamlss::is.gamlss(object))) {
    stop("Please provide a valid ds.gamlss.object to derive the predictions.")
  }
  if (is.null(object$dataname)) {
    stop("The function ds.predict.gamlss is not compatible with the output from gamlss(). Please provide the output from ds.gamlss() instead.")
  }

  # newdata
  if (!is.data.frame(newdata)) {
    stop("Please provide a valid data.frame as newdata that can be used to derive the predictions from the ds.gamlss.object.")
  }

  # what
  if (!what %in% object$parameters) {
    stop(paste("The parameter what=", what, " is not among the distribution parameters from the family of the ds.gamlss.object ", paste(object$parameters, collapse = ", "), ".", sep = ""))
  }

  # type
  if (!type %in% c("link", "response")) {
    stop(paste("The specified type ", type, " should either be 'link' or 'response'.", sep = ""))
  }

  pb <- utils::getFromNamespace("pb", "gamlss")
  pb.control <- utils::getFromNamespace("pb.control", "gamlss")

  ## Get estimated model parameters for distribution parameter what
  par.coef <- eval(parse(text = paste("object$", what, ".coefficients", sep = "")), envir = environment())
  par.coef.names <- names(par.coef)
  par.coefSmo <- eval(parse(text = paste("object$", what, ".coefSmo", sep = "")), envir = environment())

  ## Calculate linear predictor for parametric effects
  # Get design matrix
  par.formula <- eval(parse(text = paste("object$", what, ".formula", sep = "")), envir = environment())
  par.formula <- stats::as.formula(gsub(object$dataname, "newdata", Reduce(paste, deparse(par.formula))))
  if (length(par.formula) == 3) {
    par.formula[2] <- NULL
  }
  par.terms <- stats::terms(par.formula)
  par.x <- stats::model.matrix(par.terms, data = newdata, contrasts = object$contrasts)
  # Calculate linear predictor
  eta <- par.x %*% par.coef

  ## Add smoothing fitted values (if included in the model)
  if (length(par.coefSmo) > 0) {
    ## get the control parameters for the smoothers
    par.smoother.coef <- par.coef.names[grep(pattern = "pb(", x = tolower(par.coef.names), fixed = TRUE)]
    # only keep the arguments for the pb() function
    par.pb.args <- substr(par.smoother.coef, start = 4, stop = nchar(par.smoother.coef) - 1)
    par.pb.args <- strsplit(par.pb.args, split = ",", fixed = TRUE)

    for (i in 1:length(par.coefSmo)) {
      ## Get the required variables
      name <- eval(parse(text = paste("object$", what, ".coefSmo[[", i, "]]$name", sep = "")), envir = environment())
      # remove explicit reference to old data frame if necessary
      if (grepl("$", name, fixed = TRUE)) {
        name <- strsplit(name, split = "$", fixed = TRUE)[[1]][2]
      }
      knots <- eval(parse(text = paste("object$", what, ".coefSmo[[", i, "]]$knots", sep = "")), envir = environment())
      args <- par.pb.args[[i]]
      if (length(grep(pattern = "pb.control", x = args, fixed = TRUE)) > 0) {
        # control parameters specified
        text <- args[grep(pattern = "pb.control", x = args, fixed = TRUE)]
        text <- sub(".*(pb\\.control\\(.*\\))", "\\1", text)
        pb.control <- eval(parse(text = paste(text)), envir = asNamespace("gamlss"))
      } else {
        # no control parameters specified - use default
        pb.control <- eval(parse(text = "pb.control()"), envir = asNamespace("gamlss"))
      }
      x <- eval(parse(text = paste("newdata$", name, sep = "")), envir = parent.frame())

      ## Create the basis matrix
      par.B <- create.bbase(xvalues = x, knots = knots, ndx = pb.control$inter, deg = pb.control$degree)

      ## Add to linear predictor
      eta <- eta + par.B %*% par.coefSmo[[i]]$coef
    }
  }

  ## Predict the parameter
  if (type == "response") {
    # get the link function
    family <- gamlss.dist::as.family(eval(parse(text = paste(object$family[1], "()", sep = "")), envir = asNamespace("gamlss.dist")))
    par <- eval(parse(text = paste("family$", what, ".linkinv(eta)", sep = "")), envir = environment())
    return(par)
  } else {
    return(eta)
  }
}
