# Checking whether the node is a root

This field is `TRUE` if and only if a root of the forest.

## Value

`TRUE` if and only if the node is a root of the forest.

## See also

[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# check whether the node is a root
node$is_root
#> [1] FALSE

# get the node corresponding to the cell whose identifier is 1, i.e,
# the root
node <- forest$get_node(1)

# check whether the node is a root
node$is_root
#> [1] TRUE
```
