# Getting the CNA chromosome position

This method returns the position in chromosome where the CNA occurred.

## Value

The position in chromosome where the CNA occurred.

## Examples

``` r
# create an amplification CNA
cna <- Amplification("X", 20002, 100, 1, 0)

# get the position in chromosome where `cna` occurs (i.e., 20002)
cna$get_position_in_chromosome()
#> [1] 20002
```
