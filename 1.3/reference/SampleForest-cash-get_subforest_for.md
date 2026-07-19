# Building sub-forests

This method builds a sub-forest using as leaves some of the original
samples.

## Arguments

- sample_names:

  The names of the samples whose cells will be used as leaves of the new
  forest.

## Value

A sample forest built on the samples mentioned in `sample_names`

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# get the subforest for sample "S_1_2"
forest$get_subforest_for("S_1_2")
#> SampleForest
#>   # of trees: 1
#>   # of nodes: 880
#>   # of leaves: 119
#>   samples: {"S_1_2"}
#> 
```
