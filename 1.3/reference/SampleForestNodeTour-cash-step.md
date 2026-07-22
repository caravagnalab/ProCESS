# Moving to the next node in the tour

This method moves to the next node in the tour.

## See also

[`SampleForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md),
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`vignette("node_tour")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_tour.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# build a tour for the forest nodes
node_tour <- get_node_tour(forest)

# show the first node in the tour
node_tour$node
#> SampleForestNode(cell_id = 1, species = "A[E1]")

# move to the next node in the tour
node_tour$step()

# show the second node in the tour
node_tour$node
#> SampleForestNode(cell_id = 2, species = "A[E1]")
```
