
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dsGamlssClient

<!-- badges: start -->
<!-- badges: end -->

The dsGamlssClient package is a DataSHIELD client-side package that
includes the server-side functions to fit Generalized Additive Models
for Location, Scale and Shape (GAMLSS) \[1\] using DataSHIELD.

DataSHIELD is a software package which allows you to do non-disclosive
federated analysis on sensitive data. The DataSHIELD website
(<https://www.datashield.org>) has in depth descriptions of what it is,
how it works and how to install it. A key point to highlight is that
DataSHIELD has a client-server infrastructure, so the dsGamlssClient
package needs to be used in conjunction with the dsGamlss package
(<https://github.com/bips-hb/dsGamlss>) - trying to use one without the
other makes no sense.

Detailed instructions on how to install DataSHIELD are at
<https://www.datashield.org/wiki>.

Discussion and help with using DataSHIELD can be obtained from the
DataSHIELD Forum <https://datashield.discourse.group/>

## Installation

You can install the development version of dsGamlssClient from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bips-hb/dsGamlssClient")
```

## References

1.  Rigby RA, Stasinopoulos DM. Generalized additive models for
    location, scale and shape. Journal of the Royal Statistical Society:
    Series C (Applied Statistics). 2005;54(3):507-54.
