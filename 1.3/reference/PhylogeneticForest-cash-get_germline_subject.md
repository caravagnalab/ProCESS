# Getting the germline subject

This method returns a data frame reporting the germline subject name
(column "sample"), population (column "pop"), super-population (column
"super_pop"), and gender (column "gender").

## Value

The name of the subject whose germline is used.

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the germline subject
forest$get_germline_subject()
#>    sample pop super_pop gender
#> 1 NA20513 TSI       EUR   male
```
