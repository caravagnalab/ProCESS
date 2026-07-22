# Getting the name of the corresponding cell species

This field is the name of the corresponding cell species.

## Value

The name of the corresponding cell species.

## See also

[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# get the node corresponding to the cell whose identifier is 2
node <- forest$get_node(2)

# get the name of the corresponding cell species
node$species_name
#> [1] "A[E1]"
```
