# Getting the parent node

This property stores the parent of the node.

This field stores the parent of the node.

## Value

The parent of the node.

The parent of the node.

## See also

[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the node's parent
node$parent
#> SampleForestNode(cell_id = 1, species = "A[E1]")
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the node's parent
node$parent
#> SampleForestNode(cell_id = 1, species = "A[E1]")
```
