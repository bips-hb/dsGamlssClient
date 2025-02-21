#'
#' @title Check existence of variables in \code{ds.gamlss} model on each server
#' @description This is an internal function required by the client function \code{{ds.gamlss}}
#' to verify that all variables specified in the formulas for the \code{ds.gamlss} model exist on the server and
#' to ensure the process does not halt inadvertently.
#' @details The variables in the formulas are checked to ensure they exist and are not empty, i.e., not missing
#' completely.
#'
#' @param formula A string, specifying the model for the mu distribution parameter. The response
#' is on the left of an ~ operator, and the terms, separated by + operators, are on the right. Currently, only
#' penalized beta splines, indicated by \code{pb()}, are supported for nonparametric smoothing,
#' e.g. \code{'y~pb(x1)+x2+x2*x3'}.
#' @param sigma.formula A string, specifying the model for the sigma distribution parameter, as in \code{formula}.
#' The only difference is, that it is not necessary to specify the response variable, e.g. \code{sigma.formula='~pb(x)'}.
#' @param nu.formula A string, specifying the model for the nu distribution parameter, as in \code{formula}.
#' The only difference is, that it is not necessary to specify the response variable, e.g. \code{nu.formula='~pb(x)'}.
#' @param tau.formula A string, specifying the model for the tau distribution parameter, as in \code{formula}.
#' The only difference is, that it is not necessary to specify the response variable, e.g. \code{tau.formula='~pb(x)'}.
#' @param data A string, specifying the name of an (optional) data frame on the server-side containing the variables occurring in the formulas.
#' If this is missing, the variables should be on the parent environment on the server-side or referenced explicitly as \code{dataname$varname}.
#' @param datasources A list of \code{\link[DSI]{DSConnection-class}}
#' objects obtained after login. If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link[DSI]{datashield.connections_default}}.
#' @keywords internal
#' @return An integer, 0 if check was passed and 1 if the tests failed
#' @author Annika Swenne
#'
gamlssChecks <- function(formula, sigma.formula, nu.formula, tau.formula, data, datasources) {
  # turn the formulas into a character
  formula <- paste0(Reduce(paste, deparse(formula)))
  sigma.formula <- paste0(Reduce(paste, deparse(sigma.formula)))
  nu.formula <- paste0(Reduce(paste, deparse(nu.formula)))
  tau.formula <- paste0(Reduce(paste, deparse(tau.formula)))
  formulas <- paste(formula, sigma.formula, nu.formula, tau.formula, sep = "|")

  # replace the symbols '~', '+' and '*' by a separator
  formulas <- gsub(" ", "", formulas, fixed = TRUE)
  formulas <- gsub("pb(", "", formulas, fixed = TRUE)
  formulas <- gsub("(1", "", formulas, fixed = TRUE)
  formulas <- gsub("(0", "", formulas, fixed = TRUE)
  formulas <- gsub("(", "", formulas, fixed = TRUE)
  formulas <- gsub(")", "", formulas, fixed = TRUE)
  formulas <- gsub("~1", "|", formulas, fixed = TRUE)
  formulas <- gsub("~0", "|", formulas, fixed = TRUE)
  formulas <- gsub("~", "|", formulas, fixed = TRUE)
  formulas <- gsub("+", "|", formulas, fixed = TRUE)
  formulas <- gsub("*", "|", formulas, fixed = TRUE)
  formulas <- gsub("/", "|", formulas, fixed = TRUE)
  formulas <- gsub(":", "|", formulas, fixed = TRUE)
  formulas <- gsub("||", "|", formulas, fixed = TRUE)

  # split the input formulas by "|" to obtain the names of the variables
  elts <- unlist(strsplit(formulas, split = "|", fixed = TRUE))
  elts <- elts[which(nchar(elts) > 0)]
  elts <- unique(elts)

  # check that each variable is defined and not empty and each study.
  # Stop the process if any check fails
  # the check for the dataframe was already included in ds.gamlss
  stdnames <- names(datasources)
  if (length(elts) > 0) {
    extractobj <- extract(elts)
    for (i in 1:length(extractobj)) {
      holder <- extractobj$holders[i]
      element <- extractobj$elements[i]
      # check that the holder exists on each server
      if (is.na(holder)) {
        holder <- data
      }
      if (is.null(holder)) {
        stop(paste0("No data.frame for the column ", element, " given. Specify it explicitly as dataname$", element, " or provide a valid data argument."))
      }
      isDefined(datasources, holder)
      for (j in 1:length(datasources)) {
        # check that the holder has the respective element as a column on each server
        cally <- call("colnamesDS", holder)
        colnames <- unlist(DSI::datashield.aggregate(datasources[j], cally))
        if (!(element %in% colnames)) {
          stop(paste0("'", element, "' is not defined in ", stdnames[j], "!"), call. = FALSE)
        } else {
          # check that the element is not missing completely (only has NA values)
          call0 <- paste0("isNaDS(", holder, "$", element, ")")
          out1 <- DSI::datashield.aggregate(datasources[j], as.symbol(call0))
          if (out1[[1]]) {
            stop("The variable ", elts[i], " in ", stdnames[j], " is missing at complete (all values are 'NA').", call. = FALSE)
          }
        }
      }
    }
  }
}
