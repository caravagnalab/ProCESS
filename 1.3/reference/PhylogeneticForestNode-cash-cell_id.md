# Getting the identifier of the associated cell

This property stores the identifier of the cell associated to the node.

## Value

The identifier of the cell associated to the node.

## Examples

``` r
# use a phylogenetic forest example
forest <- ProCESS::example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the identifier of the cell associated to the node
node$cell_id
#> [1] 2
```
