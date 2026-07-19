# Getting the corresponding cell death time

This property is the simulated time at which the corresponding time
died.

## Value

The simulated time at which the corresponding time died or the sampling
time.

## Examples

``` r
# use a phylogenetic forest example
forest <- ProCESS::example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the death time
node$death_time
#> [1] 15.85096
```
