# Generalized Additive Models for Location Scale and Shape using DataSHIELD

Fits a Generalized Additive Model for Location, Scale and shape (GAMLSS)
using DataSHIELD on data from a single source or multiple sources on the
server side.

## Usage

``` r
ds.gamlss(
  formula = NULL,
  sigma.formula = ~1,
  nu.formula = ~1,
  tau.formula = ~1,
  family = "NO()",
  data = NULL,
  min.values = NULL,
  max.values = NULL,
  min.max.names = NULL,
  checks = FALSE,
  mu.fix = FALSE,
  sigma.fix = FALSE,
  nu.fix = FALSE,
  tau.fix = FALSE,
  control = c(0.001, 20, 1, 1, 1, 1, Inf),
  i.control = c(0.001, 50, 30, 0.001),
  autostep = TRUE,
  datasources = NULL
)
```

## Arguments

- formula:

  A formula object, specifying the model for the mu distribution
  parameter. The response is on the left of an ~ operator, and the
  terms, separated by + operators, are on the right. Currently, only
  penalized beta splines, indicated by `pb()`, are supported for
  nonparametric smoothing, e.g. `y~pb(x1)+x2+x2*x3`.

- sigma.formula:

  A formula object, specifying the model for the sigma distribution
  parameter, as in `formula`. The only difference is, that it is not
  necessary to specify the response variable, e.g.
  `sigma.formula=~pb(x)`.

- nu.formula:

  A formula object, specifying the model for the nu distribution
  parameter, as in `formula`. The only difference is, that it is not
  necessary to specify the response variable, e.g. `nu.formula=~pb(x)`.

- tau.formula:

  A formula object, specifying the model for the tau distribution
  parameter, as in `formula`. The only difference is, that it is not
  necessary to specify the response variable, e.g. `tau.formula=~pb(x)`.

- family:

  A string, specifying the distribution of the response variable and the
  link functions of the distribution parameters. Currently, the
  following families are supported:
  `family=c('NO()', 'NO2()', 'BCCG()', 'BCPE()')`. Details on the
  distributions can be found in
  [`gamlss.family`](https://rdrr.io/pkg/gamlss.dist/man/gamlss.family.html).
  Default `family='NO()'`.

- data:

  A string, specifying the name of an (optional) data frame on the
  server-side containing the variables occurring in the formulas. If
  this is missing, the variables should be on the parent environment on
  the server-side or referenced explicitly as `dataname$varname`.

- min.values:

  A numeric vector specifying minimum values for the covariates, which
  are used to determine the knots for `pb()`. If `min.values=NULL` an
  anonymized (noisy) minimum is used instead to determine the knots on
  all servers. Default `min.values=NULL`.

- max.values:

  A numeric vector specifying maximum values for the covariates, which
  are used to determine the knots for `pb()`. If `max.values=NULL` an
  anonymized (noisy) maximum is used instead to determine the knots on
  all servers. Default `min.values=NULL`.

- min.max.names:

  A string vector specifying the names for the minimum (`min.values`)
  and maximum values (`max.values`). Only required if `min.values` or
  `max.values` are given. Default `min.max.names=NULL`.

- checks:

  Logical, if `checks=TRUE` `ds.gamlss` checks whether the required
  variables for the model exist on each server and are not completely
  missing. Default `checks=FALSE`.

- mu.fix:

  Logical, indicating whether the mu distribution parameter should be
  kept fixed during the fitting processes. Default `mu.fix=FALSE`.

- sigma.fix:

  Logical, indicating whether the sigma distribution parameter should be
  kept fixed during the fitting processes. Default `sigma.fix=FALSE`.

- nu.fix:

  Logical, indicating whether the nu distribution parameter should be
  kept fixed during the fitting processes. Default `nu.fix=FALSE`.

- tau.fix:

  Logical, indicating whether the tau distribution parameter should be
  kept fixed during the fitting processes. Default `tau.fix=FALSE`.

- control:

  Numeric vector with seven elements that sets the control parameters of
  the outer iterations algorithm using the
  [`gamlss.control`](https://rdrr.io/pkg/gamlss/man/gamlss.control.html)
  function: (i) c.crit (the convergence criterion for the
  algorithm), (ii) n.cyc (the number of cycles of the algorithm), (iii)
  mu.step (the step length for the distribution parameter mu), (iv)
  sigma.step (the step length for the distribution parameter sigma), (v)
  nu.step (the step length for the distribution parameter nu), (vi)
  tau.step (the step length for the distribution parameter tau), (vii)
  gd.tol (global deviance tolerance level). Default
  `control=c(0.001, 20, 1, 1, 1, 1, Inf)`.

- i.control:

  Numeric vector with four elements that sets the control parameters of
  the inner iterations of the RS algorithm using the
  [`glim.control`](https://rdrr.io/pkg/gamlss/man/glim.control.html)
  function: (i) cc (the convergence criterion for the algorithm), (ii)
  cyc (the number of cycles of the algorithm), (iii) bf.cyc (the number
  of cycles of the backfitting algorithm), (iv) bf.tol (the convergence
  criterion (tolerance level) for the backfitting algorithm). Default
  `i.control=c(0.001, 50, 30, 0.001)`.

- autostep:

  Logical, indicating whether the steps should be halved automatically
  if the new global deviance is greater than the old one. Default
  `autostep=TRUE`.

- datasources:

  A list of
  [`DSConnection-class`](https://datashield.github.io/DSI/reference/DSConnection-class.html)
  objects obtained after login. If the `datasources` argument is not
  specified the default set of connections will be used: see
  [`datashield.connections_default`](https://datashield.github.io/DSI/reference/datashield.connections_default.html).

## Value

A ds.gamlss object with all components as in the
[`gamlss`](https://rdrr.io/pkg/gamlss/man/gamlss.html) function.
Individual-level information like the components `y` (the response) and
`residuals` (the normalised quantile residuals of the model) are not
disclosed to the client-side.

## Details

Fits a Generalized Additive model for Location, scale and shape (GAMLSS)
using DataSHIELD on data from a single source or multiple sources on the
server side. In the latter case, the data are co-analysed (when using
`ds.gamlss`) by using an approach that is mathematically equivalent to
placing all individual-level data from all sources in one central
warehouse and analysing those data using the conventional
[`gamlss`](https://rdrr.io/pkg/gamlss/man/gamlss.html) function in R.
For additional details please see the header of the
[`gamlss`](https://rdrr.io/pkg/gamlss/man/gamlss.html) function.

Server functions called: `gamlssDS1`, `gamlssDS2`, `gamlssDS3`,
`gamlssDS4`, `gamlssDS5`, `gamlssDS6`

## Author

Annika Swenne

## Examples

``` r
library(DSLite)
#> Loading required package: DSI
#> Loading required package: progress
#> Loading required package: R6
#> Loading required package: rly
data(mtcars)

## Set up DSLite server
dslite.server1 <- newDSLiteServer(
  tables = list(data = mtcars[c(1:15), ]),
  config = defaultDSConfiguration(include = c("dsBase", "dsGamlss"))
)
dslite.server2 <- newDSLiteServer(
  tables = list(data = mtcars[c(16:nrow(mtcars)), ]),
  config = defaultDSConfiguration(include = c("dsBase", "dsGamlss"))
)
builder <- DSI::newDSLoginBuilder()
builder$append(server = "study1", url = "dslite.server1", table = "data", driver = "DSLiteDriver")
builder$append(server = "study2", url = "dslite.server2", table = "data", driver = "DSLiteDriver")
logindata.dslite <- builder$build()
# Login to the virtualized server
conns <- DSI::datashield.login(logindata.dslite, assign = TRUE)
#> 
#> Logging into the collaborating servers
#> Error in base::get(url, envir = getOption("datashield.env", parent.frame())): object 'dslite.server1' not found
DSI::datashield.assign.table(conns = conns, symbol = "D", table = c("data", "data"))
#> Error: object 'conns' not found

## Examples
# Example 1: parametric model
model1 <- ds.gamlss(formula = mpg ~ wt, sigma.formula = ~disp, data = "D", family = "NO()")
#> Error:  Are you logged in to any server? Please provide a valid DSConnection object! 

# Example 2: penalized beta splines
model2 <- ds.gamlss(formula = mpg ~ pb(wt), sigma.formula = ~disp, data = "D", family = "NO()")
#> Error:  Are you logged in to any server? Please provide a valid DSConnection object! 

# Example 3: penalized beta splines with known minimum and maximum
model3 <- ds.gamlss(
  formula = mpg ~ pb(wt), sigma.formula = ~disp,
  min.values = min(mtcars$wt),
  max.values = max(mtcars$wt),
  min.max.names = "wt",
  data = "D", family = "NO()"
)
#> Error:  Are you logged in to any server? Please provide a valid DSConnection object! 

## Logout
DSI::datashield.logout(conns)
#> Error: object 'conns' not found
```
