# Plot a sample forest

Plot a sample forest. This plot is carried out using `ggraph` and for
simplicity of visualisation the forest is plot as a set of trees
connected to a generic wildtype cell.

## Usage

``` r
plot_forest(forest, highlight_sample = NULL, color_map = NULL)
```

## Arguments

- forest:

  The sample forest to be plot.

- highlight_sample:

  If a sample name, the path from root to the sampled cells in the
  sample is highlighted. If `NULL` (default), nothing is highlighted.

- color_map:

  A named vector representing the simulation species color map
  (optional).

## Value

A `ggraph` tree plot.

## Examples

``` r
# use a sample forest example
forest <- example("SampleForest")

# plot the forest
plot_forest(forest)


# define a custom color map for the forest species
color_map <- c("#7FC97F", "#BEAED4", "#FDC086", "#FFFF99")

plot_forest(forest, color_map = color_map)
```
