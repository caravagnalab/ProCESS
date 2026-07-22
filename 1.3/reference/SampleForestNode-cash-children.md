# Getting the children of the node

This property stores the children of the node.

This field stores the children of the node.

## Value

The children of the node.

The children of the node.

## See also

[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the node's children
node$children
#> [[1]]
#> SampleForestNode(cell_id = 6, species = "A[E2]")
#> 
#> [[2]]
#> SampleForestNode(cell_id = 7, species = "A[E1]")
#> 
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the node's children
node$children
#> [[1]]
#> SampleForestNode(cell_id = 6, species = "A[E2]")
#> 
#> [[2]]
#> SampleForestNode(cell_id = 7, species = "A[E1]")
#> 
```
