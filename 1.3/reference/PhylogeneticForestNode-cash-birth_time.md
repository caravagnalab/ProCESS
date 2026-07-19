# Getting the corresponding cell birth time

This property is the simulated time at which the corresponding time was
born.

## Value

The simulated time at which the corresponding time was born.

## Examples

``` r
# use a phylogenetic forest example
forest <- ProCESS::example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the birth time
node$birth_time
#> [1] 5.741436
```
