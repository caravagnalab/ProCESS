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
sim <- TissueSimulation()
sim$add_mutant("A", c(duplication = 0.08))
sim$place_cell("A", 500, 500)
sim$run_up_to_time(60)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    

sim$sample_cells("MySample", c(500, 500), c(510, 510))
forest <- sim$get_sample_forest()

plot_forest(forest)


# define a custom color map
color_map <- c("#B2DF8A")
names(color_map) <- c("A")

plot_forest(forest, color_map = color_map)
```
