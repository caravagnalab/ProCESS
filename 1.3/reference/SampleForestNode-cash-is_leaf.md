# Checking whether the node is a leaf

This property is `TRUE` if and only if a leaf of the forest.

This field is `TRUE` if and only if a leaf of the forest.

## Value

`TRUE` if and only if the node is a leaf of the forest.

`TRUE` if and only if the node is a leaf of the forest.

## See also

[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# check whether the node is a leaf
node$is_leaf
#> [1] FALSE

# get a tour over forest leaves
node_tour <- get_node_tour(forest, only_leaves = TRUE)

# check whether the first node in the tour is a leaf
node_tour$node$is_leaf
#> [1] TRUE
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# check whether the node is a leaf
node$is_leaf
#> [1] FALSE

# get a tour over forest leaves
node_tour <- get_node_tour(forest, only_leaves = TRUE)

# check whether the first node in the tour is a leaf
node_tour$node$is_leaf
#> [1] TRUE
```
