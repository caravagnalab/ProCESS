# Building sub-forests

This method builds a sub-forest using as leaves some of the original
samples.

## Arguments

- sample_names:

  The names of the samples whose cells will be used as leaves of the new
  forest

## Value

A sample forest built on the samples mentioned in `sample_names`

## See also

[`SampleForest$get_subforest_for()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_subforest_for.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the subforest for sample "S_1_2"
forest$get_subforest_for("S_1_2")
#> PhylogeneticForest
#>   # of trees: 1
#>   # of nodes: 880
#>   # of leaves: 119
#>   samples: {"S_1_2"}
#> 
#>   # of emerged SNVs and indels: 2050
#>   # of emerged CNAs: 5
#> 
```
