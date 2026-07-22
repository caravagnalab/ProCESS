# Getting the mutation data frame

This method builds a data frame representing the mutation.

## Details

The data frame has the columns `chr`, `from`, `ref`, `alt`, `type`
(i.e., `SNV` and `indel`), and `cause`.

## See also

[`Mutation`](https://caravagnalab.github.io/ProCESS/1.3/reference/Mutation_class.md)

## Examples

``` r
snv <- SNV("X", 20002, "T", "A")

snv$get_dataframe()
#>   chr  from allele ref alt type cause
#> 1   X 20002 random   A   T  SNV  <NA>
```
