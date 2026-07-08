# Getting the mutation reference sequence

This method returns the reference sequence that is altered by the
mutation.

## Value

The reference sequence before the mutation.

## Examples

``` r
snv <- SNV("X", 20002, "T", "A")

# get the reference base in which `snv` occurs (i.e., "A")
snv$get_ref()
#> [1] "A"
```
