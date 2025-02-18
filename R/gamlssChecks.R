#'
#' @title Checks if the elements in the gamlss model have the right characteristics
#' @description This is an internal function required by the client function \code{\link{ds.gamlss}}
#' to verify all the variables and ensure the process does not halt inadvertently.
#' @details the variables are checked to ensure they are defined, not empty (i.e. are not missing
#' at complete).
#' @param formula a character, a regression formula given as a string character
#' @param sigma.formula a character, a regression formula given as a string character
#' @param nu.formula a character, a regression formula given as a string character
#' @param tau.formula a character, a regression formula given as a string character
#' @param data a character, the name of an optional data frame containing the variables in
#' in the \code{formula}.
#' @param datasources a list of \code{\link[DSI]{DSConnection-class}} objects obtained after login. If the <datasources>
#' the default set of connections will be used: see \link[DSI]{datashield.connections_default}.
#' @keywords internal
#' @return an integer 0 if check was passed and 1 if failed
#' @author Annika Swenne
#'
gamlssChecks <- function(formula, sigma.formula, nu.formula, tau.formula, data, datasources){
  
  # turn the formulas into a character
  formula <- paste0(Reduce(paste, deparse(formula)))
  sigma.formula <- paste0(Reduce(paste, deparse(sigma.formula)))
  nu.formula <- paste0(Reduce(paste, deparse(nu.formula)))
  tau.formula <- paste0(Reduce(paste, deparse(tau.formula)))
  formulas <- paste(formula, sigma.formula, nu.formula, tau.formula, sep="|")
  
  # replace the symbols '~', '+' and '*' by a separator
  formulas <- gsub(" ", "", formulas, fixed=TRUE)
  formulas <- gsub("pb(", "", formulas, fixed=TRUE)
  formulas <- gsub("(1", "", formulas, fixed=TRUE)
  formulas <- gsub("(0", "", formulas, fixed=TRUE)
  formulas <- gsub("(", "", formulas, fixed=TRUE)
  formulas <- gsub(")", "", formulas, fixed=TRUE)
  formulas <- gsub("~1", "|", formulas, fixed=TRUE)
  formulas <- gsub("~0", "|", formulas, fixed=TRUE)
  formulas <- gsub("~", "|", formulas, fixed=TRUE)
  formulas <- gsub("+", "|", formulas, fixed=TRUE)
  formulas <- gsub("*", "|", formulas, fixed=TRUE)
  formulas <- gsub("/", "|", formulas, fixed=TRUE)
  formulas <- gsub(":", "|", formulas, fixed=TRUE)
  formulas <- gsub("||", "|", formulas, fixed=TRUE)

  # split the input formulas by "|" to obtain the names of the variables
  elts <- unlist(strsplit(formulas, split="|", fixed=TRUE))
  elts <- elts[which(nchar(elts)>0)]
  elts <- unique(elts)
  
  # check that each variable is defined and not empty and each study. 
  # Stop the process if any check fails
  # the check for the dataframe was already included in ds.gamlss
  stdnames <- names(datasources)
  if (length(elts)>0){
    extractobj <- extract(elts)
    for (i in 1:length(extractobj)){
      holder <- extractobj$holders[i]
      element <- extractobj$elements[i]
      # check that the holder exists on each server
      if (is.na(holder)){
        holder <- data
      }
      if (is.null(holder)){
        stop(paste0("No data.frame for the column ", element, " given. Specify it explicitly as dataname$", element, " or provide a valid data argument."))
      }
      isDefined(datasources, holder)
      for (j in 1:length(datasources)){
        # check that the holder has the respective element as a column on each server
        cally <- call("colnamesDS", holder)
        colnames <- unlist(DSI::datashield.aggregate(datasources[j], cally))
        if(!(element %in% colnames)){
          stop(paste0("'", element, "' is not defined in ", stdnames[j], "!"), call.=FALSE)
        } else {
          # check that the element is not missing completely (only has NA values)
          call0 <- paste0("isNaDS(", holder, "$", element, ")")
          out1 <- DSI::datashield.aggregate(datasources[j], as.symbol(call0))
          if(out1[[1]]){
            stop("The variable ", elts[i], " in ", stdnames[j], " is missing at complete (all values are 'NA').", call.=FALSE)
          }
        }
      }
    }
  }
}