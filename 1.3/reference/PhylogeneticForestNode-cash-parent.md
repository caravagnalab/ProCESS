# Getting the parent node

This field stores the parent of the node.

## Value

The parent of the node.

## See also

[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the node's parent
node$parent
#> PhylogeneticForestNode(cell_id = 1, species = "A[E1]")
```
