# Getting the parent node

This property stores the parent of the node.

## Value

The parent of the node.

## Examples

``` r
# use a phylogenetic forest example
forest <- ProCESS::example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the node's parent
node$parent
#> PhylogeneticForestNode(cell_id = 1, species = "A[E1]")
```
