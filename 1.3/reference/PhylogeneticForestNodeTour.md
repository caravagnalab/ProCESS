# An iterator class over phylogenetic forest nodes

This class represents iterators over phylogenetic forest nodes. The
objects of this class are built by
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md).

## Fields

- `node`:

  An object of the class
  [`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md)
  representing the node pointed by the iterator.

- `label`:

  (OPTIONAL) The label of the of the node pointed by the iterator. The
  presence of this field depends on the
  [`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)'s
  parameters used to create the tour object.

- `genome`:

  (OPTIONAL) An object of the class
  [`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md)
  that represent the genome of the node pointed by the iterator.

- `step`:

  A method that moves to the next node in the tour.

- `done`:

  A Boolean flag that is set to TRUE only when the tour ended.

## See also

[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md),
[`SampleForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md),
[`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md),
`vignette("node_labelling")`
