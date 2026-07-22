# Getting the CNA length

This method returns the CNA length.

## Value

The CNA length.

## See also

[`CNA`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA_class.md)

## Examples

``` r
# create an amplification CNA
cna <- CNA("A", "X", 20002, 100)

# get the length of `cna` (i.e., 100)
cna$get_length()
#> [1] 100
```
