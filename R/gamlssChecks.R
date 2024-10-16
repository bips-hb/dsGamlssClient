#'
#' @title Checks if the elements in the gamlss model have the right characteristics
#' @description This is an internal function required by the client function \code{\link{ds.gamlss}}
#' to verify all the variables and ensure the process does not halt inadvertanly.
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
  elts <- unique(elts)
  
  # check that each variable is defined and not empty and each study. Stop the process if any check fails
  stdnames <- names(datasources)
  for(i in 1:length(elts)){
    if(is.na(as.numeric(elts[i], options(warn=-1)))){ # making sure an eventual intercept term is not included in the checks
      message(paste0("    ", elts[i], "..."))
      for(j in 1: length(datasources)){
        # check if the variable is defined on the server site
        myterms <- unlist(strsplit(elts[i], split='$', fixed=TRUE))
        if(length(myterms) == 1){
          cally <- call("exists", myterms[1])
          out <- DSI::datashield.aggregate(datasources[j], cally)
          if(!(out[[1]])){
            stop(paste0("'", myterms[1], "' is not defined in ", stdnames[j], "!"), call.=FALSE)
          }else{
            cally <- call("colnamesDS", myterms[1])
            clnames <- unlist(DSI::datashield.aggregate(datasources[j], cally))
            if(!(myterms[2] %in% clnames)){
              stop(paste0("'", myterms[2], "' is not defined in ", stdnames[j], "!"), call.=FALSE)
            }else{
              call0 <- paste0("isNaDS(", elts[i], ")")
            }
          }
        }else{
          if(!(is.null(data))){
            cally <- call("colnamesDS", data)
            clnames <- unlist(DSI::datashield.aggregate(datasources[j], cally))
            if(!(elts[i] %in% clnames)){
              dd <- isDefined(datasources, elts[i])
              call0 <- paste0("isNaDS(", elts[i], ")")
            }else{
              call0 <- paste0("isNaDS(", paste0(data, "$", elts[i]), ")")
            }
          }else{
            defined <- isDefined(datasources, elts[i])
            call0 <- paste0("isNaDS(", elts[i], ")")
          }
        }
        # check if variable is not missing at complete
        out1 <- DSI::datashield.aggregate(datasources[j], as.symbol(call0))
        if(out1[[1]]){
          stop("The variable ", elts[i], " in ", stdnames[j], " is missing at complete (all values are 'NA').", call.=FALSE)
        }
      }
    }
  }
}