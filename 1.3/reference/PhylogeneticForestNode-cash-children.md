# Getting the children of the node

This field stores the children of the node.

## Value

The children of the node.

## See also

[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the node's children
node$children
#> [[1]]
#> PhylogeneticForestNode(cell_id = 6, species = "A[E2]")
#> 
#> [[2]]
#> PhylogeneticForestNode(cell_id = 7, species = "A[E1]")
#> 
```
