# Getting forest cell mutations

This method generates a
[`PhylogeneticForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestLabelTour.md)

## Usage

``` r
get_genome_tour(forest, only_leaves)
```

## Arguments

- forest:

  A
  [`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)
  object.

- only_leaves:

  A Boolean value (default: `FALSE`).

## Value

A
[`PhylogeneticForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestLabelTour.md)
iterates over the genome mutations of the cells associated to the
`forest`'s nodes. The returned object exclusively iterates over
`forest`'s leaves if and only if `only_leaves` is set to `TRUE`.

## See also

[`get_label_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_label_tour.md),
[`PhylogeneticForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestLabelTour.md),
[`vignette("node_labelling")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_labelling.md)
