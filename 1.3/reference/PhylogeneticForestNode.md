# The node of a phylogenetic forest

This class represents the nodes of a phylogenetic forest. It does not
have a user constructor and its objects are produced by ProCESS and
passed to the labelling function of
[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md).

## Fields

- `cell_id`:

  The identifier of the associated cell.

- `parent`:

  The node's parent.

- `children`:

  A list of the node's children.

- `is_root`:

  A flag that is set to TRUE if and only if the node is a root.

- `is_leaf`:

  A flag that is set to TRUE if and only if the node is a leaf.

- `sample_name`:

  The name of the sample that collected the associated cell.

- `birth_time`:

  The birth time of the cell associated to the node.

- `death_time`:

  The death time of the cell associated to the node.

- `life_span`:

  The life span of the cell associated to the node.

- `species_id`:

  The identifier of the associated cell's species.

- `species_name`:

  The name of the associated cell's species.

- `epistate_name`:

  The name of the associated cell's epigenetic state.

- `mutant_id`:

  The identifier of the associated cell's mutant.

- `mutant_name`:

  The name of the associated cell's mutant.

- `arising_mutations`:

  The mutations arising for the first time in the associated cell.

## See also

[`get_node_tour()`](https://caravagnalab.github.io/ProCESS/1.3/reference/get_node_tour.md),
[`PhylogeneticForestNodeTour`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNodeTour.md),
[`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md),
`vignette("node_labelling")`
