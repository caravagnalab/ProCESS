# Plot a sample forest

Plot a sample forest. This plot is carried out using `ggraph` and for
simplicity of visualisation the forest is plot as a set of trees
connected to a generic wild-type cell.

## Usage

``` r
plot_forest(
  forest,
  highlight_sample = NULL,
  color_map = NULL,
  alpha_function = NULL
)
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

- alpha_function:

  A function whose input is the data frame returned by
  `SampleForest$get_cells()` or `PhylogeneticForest$get_cells()` and
  returns a real vector whose values are in the interval \\\[0,1\]\\ and
  whose length is the number of rows in the input data frame. Each value
  in the output is used as alpha level of the corresponding cell. When
  the parameter is set to `NULL`, all tumour simulation cells have alpha
  level `1` (default: `NULL`).

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


# define an alpha function hiding the nodes representing cells
# born after the simulated time 400
library(dplyr)

alpha_f <- function(nodes) {
   nodes %>%
     dplyr::mutate(alpha = dplyr::case_when(birth_time <= 400 ~ 1,
                                            TRUE ~ 0)) %>%
     dplyr::pull(alpha)
}

plot_forest(forest, alpha_function = alpha_f)
```
