# Getting the corresponding cell birth time

This property is the simulated time at which the corresponding time was
born.

This field is the simulated time at which the corresponding time was
born.

## Value

The simulated time at which the corresponding time was born.

The simulated time at which the corresponding time was born.

## See also

[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the birth time
node$birth_time
#> [1] 5.741436
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the birth time
node$birth_time
#> [1] 5.741436
```
