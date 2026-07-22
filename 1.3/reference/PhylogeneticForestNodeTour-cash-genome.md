# Getting the genome of the associated tour

This property stores the genome of the cell associated to the current
node in the tour.

## Value

The genome of the cell associated to the current node in the tour.

## Details

This property is optional and is available only if the
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)'s
optional parameter `with_genomes` was set to `TRUE`.

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

#' # build a tour for the forest nodes
node_tour <- get_node_tour(forest, )

# show the first node in the tour
node_tour$node
#> PhylogeneticForestNode(cell_id = 1, species = "A[E1]")
```
