# The phylogenetic forest of cells in samples

This class represents the phylogenetic forest of the cells sampled
during the computation.

## Details

The leaves of his forest are the sampled cells. This class is analogous
to the class
[`SampleForest()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest.md),
but each node is labelled with the mutations occurring for the first
time on the cell represented by the node itself. Moreover each leaf is
also associated with the genome mutations occurring in the corresponding
cell.

The objects of this class are built by using
[MutationEngine](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine.md)
objects (see `[vignette("mutations")]`) and provide the following
methods and properties:

- [`get_absolute_chromosome_positions()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_absolute_chromosome_positions.md)
  returns the absolute chromosome positions.

- [`get_bulk_allelic_fragmentation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_bulk_allelic_fragmentation.md)
  returns the genome bulk allelic fragmentation.

- [`get_cell_allelic_fragmentation()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_cell_allelic_fragmentation.md)
  returns the cell allelic fragmentation.

- [`get_coalescent_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_coalescent_cells.md)
  returns the most recent common ancestors of the sampled cells.

- [`get_driver_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_driver_mutations.md)
  returns the data frame of the driver mutations.

- [`get_exposures()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_exposures.md)
  returns the timed exposure data frame.

- [`get_first_occurrences()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_first_occurrences.md)
  returns the identifier of the first cell in which a mutation emerged.

- [`get_germline_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_germline_mutations.md)
  returns the data frame of the forest germinal mutations.

- [`get_germline_subject()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_germline_subject.md)
  returns the germline subject of the phylogenetic forest.

- [`get_mutation_statistics()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_mutation_statistics.md)
  returns the statistics about mutations on each node.

- [`get_node()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_node.md)
  returns an object of type
  [`PhylogeneticForestNode`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForestNode.md).

- [`get_nodes()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_nodes.md)
  returns information about the nodes in the forest.

- [`get_reference_path()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_reference_path.md)
  returns the path to the reference genome FASTA file.

- [`get_sampled_cell_CNAs()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sampled_cell_CNAs.md)
  returns the data frame of the CNAs in the sampled cells.

- [`get_sampled_cell_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sampled_cell_mutations.md)
  returns the data frame of the SNVs and indels in the sampled cells.

- [`get_samples_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_samples_info.md)
  returns information about the samples generating the forest.

- [`get_species_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_species_info.md)
  returns information about the simulated species.

- [`get_sticks()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_sticks.md)
  computes the forest sticks.

- [`get_subforest_for()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_subforest_for.md)
  extracts a sub-forest for some of the samples.

- [`partition_samples()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-partition_samples.md)
  partitions the samples according to a labelling.

- [`represents_cell()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-represents_cell.md)
  tests whether a cell having a given identifier is represented by the
  forest.

- [`save()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-save.md)
  saves the forest.

- [`set_reference_path()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-set_reference_path.md)
  sets the path to the reference genome FASTA file.

## See also

[`SampleForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest.md),
[`MutationEngine$place_mutations()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-place_mutations.md)
