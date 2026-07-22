# The sample cell ancestor forest

The ancestor forest of sampled cells

## Details

This class represents the forest of the ancestors of the cells sampled
during the computation. The leaves of this forest are the sampled cells.

The objects of this class provide the following methods and properties:

- [`get_coalescent_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_coalescent_cells.md)
  returns the most recent common ancestors of the sampled cells.

- [`get_node()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_node.md)
  returns an object of type
  [`SampleForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForestNode.md).

- [`get_nodes()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_nodes.md)
  returns information about the nodes in the forest.

- [`get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_samples_info.md)
  returns information about the samples generating the forest.

- [`get_species_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_species_info.md)
  returns information about the simulated species.

- [`get_sticks()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_sticks.md)
  computes the forest sticks.

- [`get_subforest_for()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_subforest_for.md)
  extracts a sub-forest for some of the samples.

- [`save()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-save.md)
  saves the forest.

## See also

[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)
