# Moving to the next node in the tour

This method moves to the next node in the tour.

## See also

[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md),
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`vignette("node_tour")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_tour.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# build a tour for the forest nodes
node_tour <- get_node_tour(forest)

# show the first node in the tour
node_tour$node
#> PhylogeneticForestNode(cell_id = 1, species = "A[E1]")

# move to the next node in the tour
node_tour$step()

# show the second node in the tour
node_tour$node
#> PhylogeneticForestNode(cell_id = 2, species = "A[E1]")
```
