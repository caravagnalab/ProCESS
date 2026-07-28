# Testing whether the tour end has been reached

A Boolean flag representing the tour end.

## Details

This Boolean property is set to `TRUE` if and only if the end of the
tour has been reached.

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

# set the counter to 0
counter <- 0

# repeat until the end has been reached
while(!node_tour$done) {
  # increase the counter
  counter <- counter + 1

  # move to the next node
  node_tour$step()
}

# show the node count
counter
#> [1] 6233

# show the forest (and the number of its nodes)
forest
#> SampleForest
#>   # of trees: 1
#>   # of nodes: 6233
#>   # of leaves: 1424
#>   samples: {"S_1_1", "S_1_2", "S_2_1", "S_2_2"}
#> 
```
