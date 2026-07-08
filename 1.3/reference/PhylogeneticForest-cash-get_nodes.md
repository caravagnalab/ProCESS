# Getting the forest nodes

This method returns the nodes of the forest.

## Value

A data frame representing, for each node in the forest, the identified
(column `cell_id`), the ancestor identifier (column `ancestor`), the
node's depth (column `depth`), the name of the sample containing the
node (column `sample`), the mutant (column `mutant`), the birth time
(column `birth_time`), and, whenever the simulation has epigenetic
states, the epigenetic state (column `epistate`).

## See also

[`SampleForest$get_nodes()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_nodes.md)
for usage examples.
