
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dsGamlssClient

<!-- badges: start -->

[![R-CMD-check](https://github.com/bips-hb/dsGamlssClient/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bips-hb/dsGamlssClient/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The `dsGamlssClient` package is a
[DataSHIELD](https://www.datashield.org) client-side package that
includes the server-side functions to fit Generalized Additive Models
for Location, Scale and Shape (GAMLSS) \[1\] using DataSHIELD. It is
based on the original
[gamlss](https://cran.r-project.org/package=gamlss) implementation \[1\]
and the [dsBaseClient](https://github.com/datashield/dsBaseClient)
package \[2\].

### DataSHIELD

DataSHIELD is a software infrastructure which allows you to do
non-disclosive federated analysis on sensitive data. The [DataSHIELD
website](https://www.datashield.org) has in depth descriptions of what
it is, how it works and how to install it. A key point to highlight is
that DataSHIELD has a client-server infrastructure, so the
`dsGamlssClient` package needs to be used in conjunction with the
[dsGamlss](https://github.com/bips-hb/dsGamlss) package - trying to use
one without the other makes no sense. Detailed instructions on how to
install DataSHIELD can be found at the [DataSHIELD
Wiki](https://www.datashield.org/wiki). Discussion and help with using
DataSHIELD can be obtained from the [DataSHIELD
Forum](https://datashield.discourse.group/).

## Installation

You can install the `dsGamlssClient` package from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bips-hb/dsGamlssClient")
```

To successfully run the package, the
[dsBase](https://github.com/datashield/dsBase) and
[dsGamlss](https://github.com/bips-hb/dsGamlss) server-side packages
must be installed on each DataSHIELD server. Instructions on how to
install a server-side package on a DataSHIELD server can be found in the
[Data Manager
Section](https://wiki.datashield.org/en/getting-started/data-manager/overview)
at the DataSHIELD Wiki.

## Example

The example uses the server-less DataSHIELD implementation
[DSLite](https://cran.r-project.org/package=DSLite) \[3\] to illustrate
the use of the package’s main functions `ds.gamlss` and
`ds.predict.gamlss`. Thus, access to a DataSHIELD server is not required
to follow the example.

First, you need to install the `DSLite` package, which is available on
[CRAN](https://cran.r-project.org/).

``` r
install.packages("DSLite", repos = "https://cloud.r-project.org/")
```

Note that the `DSLite` package is only required if you want to use the
server-less DataSHIELD implementation. If instead you are using
[Armadillo](https://molgenis.github.io/molgenis-service-armadillo/) or
[Opal](https://opaldoc.obiba.org/en/latest/) DataSHIELD servers you must
to install the
[DSMolgenisArmadillo](https://cran.r-project.org/package=DSMolgenisArmadillo)
or [DSOpal](https://cran.r-project.org/package=DSOpal) package from CRAN
with `install.packages("DSMolgenisArmadillo")` or
`install.packages("DSOpal")`.

Furthermore, to follow the example, the server-side packages `dsBase`
and `dsGamlss` must be installed from GitHub in the global environment.
Again, this is only necessary for the server-less DataSHIELD
implementation. If you are using Armadillo or Opal DataSHIELD servers,
you must install the `dsBase` and `dsGamlss` packages on the server, as
described e.g. in the [Data Manager
Section](https://wiki.datashield.org/en/getting-started/data-manager/overview)
at the DataSHIELD Wiki.

``` r
# install.packages("devtools")
devtools::install_github("datashield/dsBase")
devtools::install_github("bips-hb/dsGamlss")
library(dsGamlssClient)
library(DSLite)
library(dsBase)
library(dsGamlss)
data(mtcars)
```

For illustrative purposes, the `mtcars` example data is split across two
servers. Therefore, two `DSLite` servers are set up, with the required
server-side packages `dsBase` and `dsGamlss`. Furthermore, each server
holds a subset of the `mtcars` data. Analyzing the data on the two
servers jointly with `ds.gamlss`, is mathematically equivalent to
fitting a GAMLSS model to the whole `mtcars` data. Note however, that if
nonparametric terms are included in the model, there might be slight
numerical differences between `ds.gamlss` and `gamlss::gamlss`, since
the matrix equation to obtain the regression coefficients is solved
differently.

``` r
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
#> 
#>   No variables have been specified. 
#>   All the variables in the table 
#>   (the whole dataset) will be assigned to R!
#> 
#> Assigning table data...
DSI::datashield.assign.table(conns = conns, symbol = "D", table = c("data", "data"))
```

Then, one can apply the `ds.gamlss` function to fit a GAMLSS model. In
this case, it is assumed that the minimum and maximum of `wt`, which are
used to determine the knots for the penalized beta spline, are known and
non-disclosive.

``` r
model <- ds.gamlss(
  formula = mpg ~ pb(wt), sigma.formula = ~wt,
  min.values = min(mtcars$wt),
  max.values = max(mtcars$wt),
  min.max.names = "wt",
  data = "D", family = "NO()"
)
#> Outer iteration 1...
#> Outer iteration 2...
#> Outer iteration 3...
```

However, this might not always be the case and instead, an anonymized
minimum and maximum can be used to determine the knots.

``` r
model <- ds.gamlss(formula = mpg ~ pb(wt), sigma.formula = ~wt, data = "D", family = "NO()")
#> Outer iteration 1...
#> Outer iteration 2...
#> Outer iteration 3...
```

Note that this implies that the knots are different from a GAMLSS model
that is fit to the whole `mctars` dataset, and hence the results from
`ds.gamlss` are slightly different from the GAMLSS model that is fit to
the whole data. However, this should have no impact on the
interpretation of the model.

After the `ds.gamlss` model has been fit it can be used to predict the
different distribution parameters, e.g. the mu distribution parameter.
Therefore a new data frame with the explanatory variables included in
the model must be created.

``` r
newdata <- data.frame(wt = seq(2, 5, by = 0.01))
mu.response <- ds.predict.gamlss(model, newdata, what = "mu", type = "response")
head(mu.response)
#>       [,1]
#> 1 27.71961
#> 2 27.64003
#> 3 27.56050
#> 4 27.48102
#> 5 27.40159
#> 6 27.32221
```

## References

1.  Rigby RA, Stasinopoulos DM. Generalized additive models for
    location, scale and shape. Journal of the Royal Statistical Society:
    Series C (Applied Statistics). 2005;54(3):507-54.
2.  DataSHIELD Developers (2023). *dsBaseClient: DataSHIELD Client
    Functions*. R package version 6.3.0.
3.  Marcon Y (2022). *DSLite: ‘DataSHIELD’ Implementation on Local
    Datasets*. R package version 1.4.0,
    <https://CRAN.R-project.org/package=DSLite>.
