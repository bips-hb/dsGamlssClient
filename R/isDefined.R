#'
#' @title Check if an object is defined on all servers
#' @description This is an internal function based on the internal \code{\link[dsBaseClient]{isDefined}} function from \code{dsBaseClient} (version 6.3.0).
#' @details In DataSHIELD an object included in an analysis must be defined (i.e. exists) on all servers. If not the process should halt.
#' @param datasources A list of \code{\link[DSI]{DSConnection-class}} 
#' objects obtained after login. If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link[DSI]{datashield.connections_default}}.
#' @param obj A string with the name of the object(s) to look for.
#' @param error.message A Boolean which specifies if the function should stop and return
#' an error message when the input object is not defined on one or more servers (\code{error.message=TRUE}) 
#' or if it should return a list of \code{TRUE/FALSE} indicating on which server the object is defined
#' (\code{error.message=FALSE}). Default \code{error.message=TRUE}.
#' @keywords internal
#' @return An error message if \code{error.message} argument is set to \code{TRUE} (default)
#' and if the input object is not defined on one or more servers, or a Boolean value if
#' \code{error.message=FALSE}.
#' @author Demetris Avraam for DataSHIELD Development Team
#'
isDefined <- function(datasources=NULL, obj=NULL, error.message=TRUE){
  
  inputobj <- unlist(obj)
  
  for(i in 1:length(inputobj)){
    
    extractObj <- extract(inputobj[i])
    
    if(is.na(extractObj$holders)){
      cally <- call('exists', extractObj$elements)
      out <- DSI::datashield.aggregate(datasources, cally)
    }else{
      dfname <- as.name(extractObj$holders)
      cally <- call('exists', extractObj$elements, dfname)
      out <- DSI::datashield.aggregate(datasources, cally)
    }
    
    if(error.message==TRUE & any(out==FALSE)){
      stop("The input object ", inputobj[i], " is not defined in ", paste(names(which(out==FALSE)), collapse=", "), "!" , call.=FALSE)
    }else{
      return(out)
    }
  }
}
