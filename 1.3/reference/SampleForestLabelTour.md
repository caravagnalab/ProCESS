# An iterator class over sample forest labels

This class represents iterators over sample forest labels. The objects
of this class are built by
[`get_label_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_label_tour.md).

## Fields

- `value`:

  A list `cell id`-`label` for the current node in the tour.

- `step`:

  A method that moves to the next node in the tour.

- `done`:

  A Boolean flag that is set to TRUE only when the tour ended.

## See also

[`get_label_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_label_tour.md),
[`PhylogeneticForestLabelTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestLabelTour.md)
