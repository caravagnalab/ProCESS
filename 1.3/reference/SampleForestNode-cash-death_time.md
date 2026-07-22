# Getting the corresponding cell death time

This property is the simulated time at which the corresponding time
died.

This field is the simulated time at which the corresponding time died.

## Value

The simulated time at which the corresponding time died or the sampling
time.

The simulated time at which the corresponding time died or the sampling
time.

## See also

[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the death time
node$death_time
#> [1] 15.85096
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the death time
node$death_time
#> [1] 15.85096
```
