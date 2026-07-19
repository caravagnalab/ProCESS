# Getting the parent node

This property stores the parent of the node.

## Value

The parent of the node.

## Examples

``` r
# use a sample forest example
forest <- ProCESS::example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the node's parent
node$parent
#> SampleForestNode(cell_id = 1, species = "A[E1]")
```
