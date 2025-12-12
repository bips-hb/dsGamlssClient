# Split input string by '\$'

This is an internal function based on the internal `extract` function
from `dsBaseClient` (version 6.3.0). It splits the input by the '\$'
symbol and returns the single elements as a list.

## Usage

``` r
extract(input)
```

## Arguments

- input:

  A vector or a list of strings.

## Value

a list with the following elements

- `holders`:

  The strings before the '\$' symbol.

- `elements`:

  The strings after the '\$' symbol.

## Author

DataSHIELD Development Team
