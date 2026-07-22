# Getting the mutation chromosome position

This method returns the position in the chromosome where the mutation
occurred.

## Value

The position in chromosome where the mutation occurred.

## See also

[`Mutation`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation_class.md)

## Examples

``` r
snv <- SNV("X", 20002, "T", "A")

# get the position in chromosome where `snv` occurs (i.e., 20002)
snv$get_position_in_chromosome()
#> [1] 20002
```
