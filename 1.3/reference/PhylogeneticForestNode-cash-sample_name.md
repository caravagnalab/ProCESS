# Getting the corresponding cell sample

This property is the name of the sample that collected the corresponding
cell.

## Value

The name of of the sample that collected the corresponding cell if the
corresponding cell was collected. Otherwise, `NA`.

## Details

This property is the name of the sample that collected the corresponding
cell. If the node is not a leaf, it was not collected by any sample and
this property is `NA`.

## Examples

``` r
# use a phylogenetic forest example
forest <- ProCESS::example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# the corresponding cell was not collected and the property is `NA`
node$sample_name
#> [1] NA

# get a tour over forest leaves
node_tour <- get_node_tour(forest, only_leaves = TRUE)

# in this case, the node is a leaf and then was collected
node_tour$node$sample_name
#> [1] "S_1_2"
```
