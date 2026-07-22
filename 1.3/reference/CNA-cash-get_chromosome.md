# Getting the CNA chromosome

This method returns the identifier of the chromosome where the CNA
occurred.

## Value

The identifier of the chromosome in which the CNA occurred.

## See also

[`CNA`](https://caravagnalab.github.io/ProCESS/1.3/reference/CNA_class.md)

## Examples

``` r
# create an amplification CNA
cna <- CNA("A", "X", 20002, 100)

# get the chromosome in which `cna` occurs (i.e., "X")
cna$get_chromosome()
#> [1] "X"
```
