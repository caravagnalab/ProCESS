# An iterator class over phylogenetic forest nodes

Iterators over phylogenetic forest nodes.

## Details

This class represents iterators over phylogenetic forest nodes. The
objects of this class are built by
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)
and provide the following methods and properties:

- [`node`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour-cash-node.md)
  is an object of the class
  [`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md)
  and represents the node pointed by the iterator.

- [`label`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour-cash-label.md)
  represents the label of the of the node pointed by the iterator. The
  presence of this field depends on the type of the
  [`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)'s
  parameters used to create the tour object (OPTIONAL).

- [`genome`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour-cash-genome.md)
  is an object of the class
  [`GenomeMutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/GenomeMutations.md)
  that represent the genome of the node pointed by the iterator
  (OPTIONAL).

- [`step()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour-cash-step.md)
  moves the iterator to the next node in the tour.

- [`done`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour-cash-done.md)
  is a Boolean flag that is set to `TRUE` only when the tour ended.

## See also

[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md),
`PhylogeneticForestNodeTour`,
[`vignette("node_tour")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_tour.md)
