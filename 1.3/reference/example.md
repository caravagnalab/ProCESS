# Get examples

This function loads data structure examples.

## Usage

``` r
example(name)
```

## Arguments

- name:

  The name of the data structure example that is aimed. The supported
  examples are "SampleForest" and "PhylogeneticForest".

## Value

An example object of the specified data type.

## See also

[`available_examples()`](https://caravagnalab.github.io/ProCESS/1.3/reference/available_examples.md)

## Examples

``` r
# get an example of `SampleForest` object
forest <- example("SampleForest")

# see the forest
forest
#> SampleForest
#>   # of trees: 1
#>   # of nodes: 6233
#>   # of leaves: 1424
#>   samples: {"S_1_1", "S_1_2", "S_2_1", "S_2_2"}
#> 

# see the first nodes of the forest
head(forest$get_nodes())
#>   cell_id ancestor depth mutant epistate sample birth_time
#> 1       1       NA     0      A       E1   <NA>   0.000000
#> 2       2        1     1      A       E1   <NA>   5.741436
#> 3       3        1     1      A       E1   <NA>   5.741436
#> 4       4        3     2      A       E1   <NA>   9.232361
#> 5       5        3     2      A       E1   <NA>   9.232361
#> 6       6        2     2      A       E2   <NA>  15.850958

# load an example of `PhylogeneticForest` object
forest <- example("PhylogeneticForest")
```
