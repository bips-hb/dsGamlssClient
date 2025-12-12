# Check existence of variables in `ds.gamlss` model on each server

This is an internal function required by the client function
`{ds.gamlss}` to verify that all variables specified in the formulas for
the `ds.gamlss` model exist on the server and to ensure the process does
not halt inadvertently.

## Usage

``` r
gamlssChecks(
  formula,
  sigma.formula,
  nu.formula,
  tau.formula,
  data,
  datasources
)
```

## Arguments

- formula:

  A string, specifying the model for the mu distribution parameter. The
  response is on the left of an ~ operator, and the terms, separated
  by + operators, are on the right. Currently, only penalized beta
  splines, indicated by `pb()`, are supported for nonparametric
  smoothing, e.g. `'y~pb(x1)+x2+x2*x3'`.

- sigma.formula:

  A string, specifying the model for the sigma distribution parameter,
  as in `formula`. The only difference is, that it is not necessary to
  specify the response variable, e.g. `sigma.formula='~pb(x)'`.

- nu.formula:

  A string, specifying the model for the nu distribution parameter, as
  in `formula`. The only difference is, that it is not necessary to
  specify the response variable, e.g. `nu.formula='~pb(x)'`.

- tau.formula:

  A string, specifying the model for the tau distribution parameter, as
  in `formula`. The only difference is, that it is not necessary to
  specify the response variable, e.g. `tau.formula='~pb(x)'`.

- data:

  A string, specifying the name of an (optional) data frame on the
  server-side containing the variables occurring in the formulas. If
  this is missing, the variables should be on the parent environment on
  the server-side or referenced explicitly as `dataname$varname`.

- datasources:

  A list of
  [`DSConnection-class`](https://datashield.github.io/DSI/reference/DSConnection-class.html)
  objects obtained after login. If the `datasources` argument is not
  specified the default set of connections will be used: see
  [`datashield.connections_default`](https://datashield.github.io/DSI/reference/datashield.connections_default.html).

## Value

An integer, 0 if check was passed and 1 if the tests failed

## Details

The variables in the formulas are checked to ensure they exist and are
not empty, i.e., not missing completely.

## Author

Annika Swenne
