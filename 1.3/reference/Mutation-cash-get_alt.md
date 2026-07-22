# Getting the mutation altered sequence

This method returns the sequence after the mutation occurred.

## Value

The sequence after the mutation occurred.

## See also

[`Mutation`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation_class.md)

## Examples

``` r
snv <- SNV("X", 20002, "T", "A")

# get the sequence after `snv` occurs (i.e., "T")
snv$get_alt()
#> [1] "T"
```
