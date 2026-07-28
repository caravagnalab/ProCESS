# Getting node labelling

The current node labelling.

## Details

This property represents the labelling of the node currently pointed by
the tour and is optional. It is available only if the
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)'s
optional parameter `labelling_functor` was set.

## See also

[`SampleForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md),
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`vignette("node_tour")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_tour.md)

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

l_funct <- function(parent_label, node) {
  node$cell_id
}

# build a tour for the forest nodes
node_tour <- get_node_tour(forest, labelling_functor = l_funct)

# get the label of the first node in the tour
node_tour$label
#> [1] 1

# move to the next node
node_tour$step()

# get the label of the second node in the tour node
node_tour$label
#> [1] 2
```
