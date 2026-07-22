# The node of a sample forest

A class representing the nodes of a sample forest.

## Details

This class represents the nodes of a sample forest. It does not have a
user constructor. Its objects are produced by
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)
and
[`SampleForest$get_node()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_node.md).

The objects of this class provide the following methods and properties:

- [`cell_id`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-cell_id.md)
  represents the identifier of the associated cell.

- [`parent`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-parent.md)
  represents the node's parent.

- [`children`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-children.md)
  represents a list of the node's children.

- [`is_root`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-is_root.md)
  is a Boolean flag that is `TRUE` if and only if the node is a forest
  root.

- [`is_leaf`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-is_leaf.md)
  is a Boolean flag that is `TRUE` if and only if the node is a forest
  leaf.

- [`sample_name`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-sample_name.md)
  is the name of the sample that collected the associated cell.

- [`birth_time`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-birth_time.md)
  is the birth time of the cell associated to the node.

- [`death_time`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-death_time.md)
  is the death time of the cell associated to the node.

- [`life_span`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-life_span.md)
  is the life span of the cell associated to the node.

- [`species_name`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-species_name.md)
  is the name of the associated cell's species.

- [`epistate_name`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-epistate_name.md)
  is the name of the associated cell's epigenetic state.

- [`mutant_name`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode-cash-mutant_name.md)
  is the name of the associated cell's mutant.

## See also

[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`SampleForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNodeTour.md),
[`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md),
[`vignette("node_tour")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_tour.md)
