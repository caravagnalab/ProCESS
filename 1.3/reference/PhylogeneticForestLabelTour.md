# An iterator class over sample forest labels

This class represents iterators over sample forest labels. The objects
of this class are built by
[`get_label_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_label_tour.md).

## Fields

- `value`:

  A list consisting of the corresponding cell identifier, the label,
  and, whenever requested, the corresponding cell genome of the current
  node in the tour.

- `step`:

  A method that moves to the next node in the tour.

- `done`:

  A Boolean flag that is set to TRUE only when the tour ended.

## See also

[`get_label_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_label_tour.md),
[`get_genome_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_genome_tour.md),
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md),
[`SampleForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestLabelTour.md),
[`vignette("node_labelling")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_labelling.md)
