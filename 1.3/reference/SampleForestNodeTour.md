# An iterator class over sample forest nodes

Iterators over sample forest nodes.

## Details

This class represents iterators over sample forest nodes. The objects of
this class are built by
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)
and provide the following methods and properties:

- [`node`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour-cash-node.md)
  is an object of the class
  [`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)
  and represents the node pointed by the iterator.

- [`label`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour-cash-label.md)
  (OPTIONAL) represents the label of the of the node pointed by the
  iterator. The presence of this field depends on the type of the
  [`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)'s
  parameters used to create the tour object.

- [`step()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour-cash-step.md)
  moves the iterator to the next node in the tour.

- [`done`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour-cash-done.md)
  is a Boolean flag that is set to `TRUE` only when the tour ended.

## See also

[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md),
[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md),
[`vignette("node_tour")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_tour.md)
