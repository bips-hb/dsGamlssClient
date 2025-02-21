#'
#' @title Gets a pooled variance from the servers
#' @description This is an internal function based on the internal \code{\link[dsBaseClient]{getPooledVar}} function from \code{dsBaseClient} (version 6.3.0).
#' @details This function is called to avoid calling the client function \code{\link[dsBaseClient]{ds.var}}
#' which may stop the process due to some checks not required when computing a mean inside a function.
#' @param datasources A list of \code{\link[DSI]{DSConnection-class}}
#' objects obtained after login. If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link[DSI]{datashield.connections_default}}.
#' @param x A string with the name of a numeric vector for which the variance should be computed.
#' @keywords internal
#' @return A numeric value, giving the pooled variance.
#' @author DataSHIELD Development Team
#'
getPooledVar <- function(datasources, x) {
  num.sources <- length(datasources)

  cally <- paste0("varDS(", x, ")")
  out.var <- DSI::datashield.aggregate(datasources, as.symbol(cally))

  length.total <- 0
  sum.weighted <- 0
  var.global <- NA

  for (i in 1:num.sources) {
    if ((!is.null(out.var[[i]][[5]])) & (out.var[[i]][[5]] != 0)) {
      var.local <- out.var[[i]][[2]] / (out.var[[i]][[4]] - 1) - (out.var[[i]][[1]])^2 / (out.var[[i]][[4]] * (out.var[[i]][[4]] - 1))
      completeLength <- out.var[[i]][[5]] - out.var[[i]][[3]]
      length.total <- length.total + completeLength
      sum.weighted <- sum.weighted + completeLength * var.local
    }
  }

  var.global <- sum.weighted / length.total
  return(var.global)
}
