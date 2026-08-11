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
  alpha_function = NULL,
  shape_label_function = sample_shape,
  color_label_function = NULL
)
```

## Arguments

- forest:

  The sample forest to be plot.

- highlight_sample:

  If a sample name, the path from root to the sampled cells in the
  sample is highlighted. If `NULL` (default), nothing is highlighted.

- color_map:

  A named vector representing the simulation species color map (optional
  when `color_label_function` is `NULL`; mandatory, otherwise).

- alpha_function:

  A function whose input is the data frame returned by
  `SampleForest$get_cells()` or `PhylogeneticForest$get_cells()` and
  returns a real vector whose values are in the interval \\\[0,1\]\\ and
  whose length is the number of rows in the input data frame. Each value
  in the output is used as alpha level of the corresponding cell. When
  the parameter is set to `NULL`, all nodes have alpha level `1`
  (default: `NULL`).

- shape_label_function:

  A function whose input is the data frame returned by
  `SampleForest$get_cells()` or `PhylogeneticForest$get_cells()` and
  returns a string vector whose length is the number of rows in the
  input data frame. Each value in the output is the label of the
  associated node and determine the node shape. When the parameter is
  set to `NULL`, all nodes have the same shape. By default, the node
  shapes correspond to the sample which collected the corresponding
  cell.

- color_label_function:

  A function whose input is the data frame returned by
  `SampleForest$get_cells()` or `PhylogeneticForest$get_cells()` and
  returns a string vector whose length is the number of rows in the
  input data frame. Each value in the output is the label of the
  associated node and determine the node color. When the parameter is
  set to `NULL`, the nodes are colored according their species (default:
  `NULL`).

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


# define an alpha function assigning to each node an alpha
# value according to the corresponding cell's birth time
alpha_f <- function(nodes) {
  library(dplyr)

  T <- max(nodes$birth_time)
  nodes %>%
    mutate(alpha = (T-.data$birth_time)/T) %>%
    pull(alpha)
}

plot_forest(forest, alpha_function = alpha_f)


# plot the forest avoiding to denote samples by using node shapes
plot_forest(forest, shape_label_function = NULL)


# a shape function that labels each node as the sample that
# collected the corresponding cells
shape_label_function <- function(nodes) {
  library(dplyr)

  nodes %>%
    mutate(label = case_when(birth_time <= 300 ~ "Old",
                             TRUE ~ "Young")) %>%
    pull(label)
}

# plot the forest using the node shapes to denote cell birth time
plot_forest(forest, shape_label_function = shape_label_function)


color_label_function <- function(nodes) {
  library(dplyr)

  nodes %>%
    mutate(age = case_when(birth_time <= 300 ~ "Old",
                           TRUE ~ "Young")) %>%
    mutate(label = paste(mutant, age)) %>%
    pull(label)
}

# get the plot labels (i.e., paste(mutants, "Old") +
# paste(mutants, "Young"))
mutants <- unique(forest$get_nodes() %>% dplyr::pull(mutant))
labels <- c(paste(mutants, "Old"), paste(mutants, "Young"))

# create a color map for the labels
color_map <- RColorBrewer::brewer.pal(n = length(labels), name = "Set1")
names(color_map) <- labels

# plot the tissue labelling the cells according to
# `label_function`. The parameter `color_map` is mandatory.
plot_forest(forest, color_label_function = color_label_function,
            color_map = color_map, shape_label_function = NULL)
```
