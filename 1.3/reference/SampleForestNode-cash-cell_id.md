# Getting the identifier of the associated cell

This property stores the identifier of the cell associated to the node.

This field stores the identifier of the cell associated to the node.

## Value

The identifier of the cell associated to the node.

The identifier of the cell associated to the node.

## See also

[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the identifier of the cell associated to the node
node$cell_id
#> [1] 2
# use a sample forest example
forest <- example("SampleForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the identifier of the cell associated to the node
node$cell_id
#> [1] 2
```
