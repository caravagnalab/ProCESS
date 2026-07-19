# Annotate a plot of cell divisions

It annotates a plot of cell divisions with information from sampling
times and MRCAs for all available samples

## Usage

``` r
annotate_forest(
  forest_plot,
  forest,
  samples = TRUE,
  MRCAs = TRUE,
  exposures = FALSE,
  facet_signatures = TRUE,
  drivers = TRUE,
  add_driver_label = TRUE
)
```

## Arguments

- forest_plot:

  The output of `plot_forest`.

- forest:

  The original forest object from which the input to `plot_forest`has
  been derived.

- samples:

  If `TRUE` it annotates samples.

- MRCAs:

  If `TRUE` it annotates MRCAs.

- exposures:

  If `TRUE` it annotates exposures to mutational signatures.

- facet_signatures:

  If `TRUE` and if `exposures` is `TRUE` it creates a faceted forest
  plot where the exposure to each signature is annotated on a separated
  plot.

- drivers:

  If `TRUE` it annotates drivers on the node they originated.

- add_driver_label:

  If `TRUE` and if `drivers` is `TRUE` it annotates the driver name.

## Value

A `ggraph` tree plot.

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# generate a plot for the forest
plot <- plot_forest(forest)

# annotate the forest plot
annotate_forest(plot, forest)
```
