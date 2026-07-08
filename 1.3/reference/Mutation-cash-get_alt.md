# Getting the mutation altered sequence

This method returns the sequence after the mutation occurred.

## Value

The sequence after the mutation occurred.

## Examples

``` r
snv <- SNV("X", 20002, "T", "A")

# get the sequence after `snv` occurs (i.e., "T")
snv$get_alt()
#> [1] "T"
```
