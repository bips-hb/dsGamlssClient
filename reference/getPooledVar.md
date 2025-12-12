# Gets a pooled variance from the servers

This is an internal function based on the internal `getPooledVar`
function from `dsBaseClient` (version 6.3.0).

## Usage

``` r
getPooledVar(datasources, x)
```

## Arguments

- datasources:

  A list of
  [`DSConnection-class`](https://datashield.github.io/DSI/reference/DSConnection-class.html)
  objects obtained after login. If the `datasources` argument is not
  specified the default set of connections will be used: see
  [`datashield.connections_default`](https://datashield.github.io/DSI/reference/datashield.connections_default.html).

- x:

  A string with the name of a numeric vector for which the variance
  should be computed.

## Value

A numeric value, giving the pooled variance.

## Details

This function is called to avoid calling the client function `ds.var`
which may stop the process due to some checks not required when
computing a mean inside a function.

## Author

DataSHIELD Development Team
