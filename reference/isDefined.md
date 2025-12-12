# Check if an object is defined on all servers

This is an internal function based on the internal `isDefined` function
from `dsBaseClient` (version 6.3.0).

## Usage

``` r
isDefined(datasources = NULL, obj = NULL, error.message = TRUE)
```

## Arguments

- datasources:

  A list of
  [`DSConnection-class`](https://datashield.github.io/DSI/reference/DSConnection-class.html)
  objects obtained after login. If the `datasources` argument is not
  specified the default set of connections will be used: see
  [`datashield.connections_default`](https://datashield.github.io/DSI/reference/datashield.connections_default.html).

- obj:

  A string with the name of the object(s) to look for.

- error.message:

  A Boolean which specifies if the function should stop and return an
  error message when the input object is not defined on one or more
  servers (`error.message=TRUE`) or if it should return a list of
  `TRUE/FALSE` indicating on which server the object is defined
  (`error.message=FALSE`). Default `error.message=TRUE`.

## Value

An error message if `error.message` argument is set to `TRUE` (default)
and if the input object is not defined on one or more servers, or a
Boolean value if `error.message=FALSE`.

## Details

In DataSHIELD an object included in an analysis must be defined (i.e.
exists) on all servers. If not the process should halt.

## Author

Demetris Avraam for DataSHIELD Development Team
