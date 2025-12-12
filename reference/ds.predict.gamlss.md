# Derive predictions for `ds.gamlss` objects

This function predicts the specified distribution parameter for a new
data set from a `ds.gamlss` object that is output by the
[`ds.gamlss`](https://bips-hb.github.io/dsGamlssClient/reference/ds.gamlss.md)
function.

## Usage

``` r
ds.predict.gamlss(object, newdata, what = "mu", type = "link")
```

## Arguments

- object:

  A `ds.gamlss` object output by
  [`ds.gamlss`](https://bips-hb.github.io/dsGamlssClient/reference/ds.gamlss.md).

- newdata:

  A data frame containing new values for the explanatory variables that
  are used in the model.

- what:

  A string, specifying the distribution parameter that should be
  predicted. Default `what="mu"`.

- type:

  A string, specifying the kind of predictions that should be derived.
  The default value `type="link"` predicts the linear predictor for the
  specified distribution parameter, whereas `type="response"` predicts
  the fitted values for the distribution parameter.

## Value

A vector with the predictions for `newdata`.

## Details

The `ds.predict.gamlss` function assumes that the object given in
`newdata` is a data frame containing the explanatory variables that are
used in the model.

## Author

Annika Swenne

## Examples

``` r
library(DSLite)
data(mtcars)

## Create newdata for predictions
newdata <- data.frame(wt = seq(2, 5, by = 0.01))

## Fit ds.gamlss model
# Set up DSLite server
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
# Fit model
model <- ds.gamlss(formula = mpg ~ pb(wt), sigma.formula = ~wt, data = "D", family = "NO()")
#> Error:  Are you logged in to any server? Please provide a valid DSConnection object! 
# Logout
DSI::datashield.logout(conns)
#> Error: object 'conns' not found

## Examples
# Example 1: Predict mu
mu.response <- ds.predict.gamlss(model, newdata, what = "mu", type = "response")
#> Error: object 'model' not found

# Example 2: Predict linear predictor for sigma
sigma.link <- ds.predict.gamlss(model, newdata, what = "sigma", type = "link")
#> Error: object 'model' not found
```
