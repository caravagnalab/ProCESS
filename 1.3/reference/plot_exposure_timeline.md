# Plot the signature exposure timeline of a phylogenetic forest

Plots the signatures exposure changes along a phylogenetic forest.

## Usage

``` r
plot_exposure_timeline(
  phylogenetic_forest,
  linewidth = 0.8,
  emphasize_switches = FALSE,
  pal_name = "Set3"
)
```

## Arguments

- phylogenetic_forest:

  A phylogenetic forest.

- linewidth:

  The width of the lines in the plot.

- emphasize_switches:

  A Boolean flag to emphasize the exposure switches.

- pal_name:

  The name of a `RColorBrewer` palette.

## Value

A `ggplot2` plot.

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# plotting the phylogenetic forest
plot_exposure_timeline(forest)


# plotting the phylogenetic forest emphasizing the exposure switches
plot_exposure_timeline(forest, emphasize_switches = TRUE)
```
