#'
#' @title Gets a pooled statistical mean from the servers
#' @description This is an internal function based on the internal \code{\link[dsBaseClient]{getPooledMean}} function from \code{dsBaseClient} (version 6.3.0).
#' @details This function is called to avoid calling the client function \code{\link[dsBaseClient]{ds.mean}}
#' which may stop the process due to some checks not required when computing a mean inside
#' a function.
#' @param datasources A list of \code{\link[DSI]{DSConnection-class}}
#' objects obtained after login. If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link[DSI]{datashield.connections_default}}.
#' @param x A string with the name of a numeric vector for which the mean should be computed
#' @keywords internal
#' @return A numeric value, giving the pooled mean.
#'
getPooledMean <- function(datasources, x) {
  num.sources <- length(datasources)

  cally <- paste0("meanDS(", x, ")")
  out.mean <- DSI::datashield.aggregate(datasources, as.symbol(cally))

  length.total <- 0
  sum.weighted <- 0
  mean.global <- NA

  for (i in 1:num.sources) {
    if ((!is.null(out.mean[[i]][[4]])) & (out.mean[[i]][[4]] != 0)) {
      completeLength <- out.mean[[i]][[4]] - out.mean[[i]][[2]]
      length.total <- length.total + completeLength
      sum.weighted <- sum.weighted + completeLength * out.mean[[i]][[1]]
    }
  }

  mean.global <- sum.weighted / length.total
  return(mean.global)
}
