# Getting the corresponding cell birth time

This field is the simulated time at which the corresponding time was
born.

## Value

The simulated time at which the corresponding time was born.

## See also

[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the birth time
node$birth_time
#> [1] 5.741436
```
