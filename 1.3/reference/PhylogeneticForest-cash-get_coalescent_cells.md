# Retrieving the most recent common ancestors

This method retrieves the most recent common ancestors of a set of
cells.

## Arguments

- cell_ids:

  The list of the identifiers of the cells whose most recent common
  ancestors are aimed (optional).

## Value

A data frame representing, for each of the identified cells, the
identified (column `cell_id`), the ancestor identifier (column
`ancestor`), the node's depth (column `depth`), the name of the sample
containing the node (column `sample`), the mutant (column `mutant`), the
birth time (column `birth_time`), and, whenever the simulation has
epigenetic states, the epigenetic state (column `epistate`).

## Details

If the optional parameter `cell_ids` is used, this method find the most
recent common ancestors of the cells having an identifier among those in
`cell_ids`. If, otherwise, the optional parameter is not used, this
method find the most recent common ancestors of the forest leaves.

## See also

[`SampleForest$get_coalescent_cells()`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest-cash-get_coalescent_cells.md)
for usage examples
