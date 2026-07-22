# Getting the corresponding cell death time

This field is the simulated time at which the corresponding time died.

## Value

The simulated time at which the corresponding time died or the sampling
time.

## See also

[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the death time
node$death_time
#> [1] 15.85096
```
