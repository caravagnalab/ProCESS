# Getting the mutation chromosome position

This method returns the position in the chromosome where the mutation
occurred.

## Value

The position in chromosome where the mutation occurred.

## Examples

``` r
snv <- SNV("X", 20002, "T", "A")

# get the position in chromosome where `snv` occurs (i.e., 20002)
snv$get_position_in_chromosome()
#> [1] 20002
```
