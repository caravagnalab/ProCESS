# Checking whether the node is a root

This property is `TRUE` if and only if a root of the forest.

This field is `TRUE` if and only if a root of the forest.

## Value

`TRUE` if and only if the node is a root of the forest.

`TRUE` if and only if the node is a root of the forest.

## See also

[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# check whether the node is a root
node$is_root
#> [1] FALSE

# get the node corresponding to the cell whose identifier is 1, i.e,
# the root
node <- forest$get_node(1)

# check whether the node is a root
node$is_root
#> [1] TRUE
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# check whether the node is a root
node$is_root
#> [1] FALSE

# get the node corresponding to the cell whose identifier is 1, i.e,
# the root
node <- forest$get_node(1)

# check whether the node is a root
node$is_root
#> [1] TRUE
```
