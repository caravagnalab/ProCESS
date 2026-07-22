# Getting the identifier of the associated cell

This field stores the identifier of the cell associated to the node.

## Value

The identifier of the cell associated to the node.

## See also

[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the identifier of the cell associated to the node
node$cell_id
#> [1] 2
```
