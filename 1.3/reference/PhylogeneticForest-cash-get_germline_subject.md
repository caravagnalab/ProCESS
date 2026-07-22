# Getting the germline subject

This method returns a data frame reporting the germline subject name
(column "sample"), population (column "pop"), super-population (column
"super_pop"), and gender (column "gender").

## Value

The name of the subject whose germline is used.

## See also

[`PhylogeneticForest$get_sampled_cell_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sampled_cell_mutations.md),
[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the germline subject
forest$get_germline_subject()
#>    sample pop super_pop gender
#> 1 NA20513 TSI       EUR   male
```
