# An iterator class over sample forest nodes

This class represents iterators over sample forest nodes. The objects of
this class are built by
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md).

## Fields

- `node`:

  An object of the class
  [`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)
  representing the node pointed by the iterator.

- `label`:

  (OPTIONAL) The label of the of the node pointed by the iterator. The
  presence of this field depends on the
  [`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)'s
  parameters used to create the tour object.

- `step`:

  A method that moves to the next node in the tour.

- `done`:

  A Boolean flag that is set to TRUE only when the tour ended.

## See also

[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md),
[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md),
`vignette("node_labelling")`
