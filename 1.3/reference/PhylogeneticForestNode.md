# The node of a phylogenetic forest

This class represents the nodes of a phylogenetic forest.

## Details

This class represents the nodes of a phylogenetic forest. It does not
have a user constructor. Its objects are produced by
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md)
and
[`PhylogeneticForest$get_node()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_node.md).

The objects of this class provide the following methods and properties:

- [`cell_id`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-cell_id.md)
  represents the identifier of the associated cell.

- [`parent`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-parent.md)
  represents the node's parent.

- [`children`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-children.md)
  represents a list of the node's children.

- [`is_root`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-is_root.md)
  is a Boolean flag that is `TRUE` if and only if the node is a forest
  root.

- [`is_leaf`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-is_leaf.md)
  is a Boolean flag that is `TRUE` if and only if the node is a forest
  leaf.

- [`sample_name`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-sample_name.md)
  is the name of the sample that collected the associated cell.

- [`birth_time`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-birth_time.md)
  is the birth time of the cell associated to the node.

- [`death_time`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-death_time.md)
  is the death time of the cell associated to the node.

- [`life_span`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-life_span.md)
  is the life span of the cell associated to the node.

- [`species_name`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-species_name.md)
  is the name of the associated cell's species.

- [`epistate_name`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-epistate_name.md)
  is the name of the associated cell's epigenetic state.

- [`mutant_name`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-mutant_name.md)
  is the name of the associated cell's mutant.

- [`arising_mutations`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-arising_mutations.md)
  stores the mutations arising in the associated cell.

- [`get_genome()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_genome.md)
  returns the genome of the associated cell.

## See also

[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md),
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md),
[`vignette("node_tour")`](https://caravagnalab.github.io/ProCESS/1.3/articles/node_tour.md)
