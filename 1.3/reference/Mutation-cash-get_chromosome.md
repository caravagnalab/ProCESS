# Getting the mutation chromosome

This method identify the chromosome where the mutation occurred.

## Value

The identifier of the chromosome in which the mutation occurred.

## See also

[`Mutation`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation_class.md)

## Examples

``` r
snv <- SNV("X", 20002, "T", "A")

# get the chromosome in which `snv` occurs (i.e., "X")
snv$get_chromosome()
#> [1] "X"
```
