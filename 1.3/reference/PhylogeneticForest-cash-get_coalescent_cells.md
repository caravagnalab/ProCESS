# Retrieving the most recent common ancestors

This method retrieves the most recent common ancestors of a set of
cells.

## Arguments

- cell_ids:

  The list of the identifiers of the cells whose most recent common
  ancestors are aimed (optional).

## Value

A data frame reporting the identified (column `cell_id`), the ancestor
identifier (column `ancestor`), the name of the sample containing the
node (column `sample`), the mutant (column `mutant`), and the birth time
(column `birth_time`). Whenever, the simulation has epigenetic states,
the data frame also contains the column `epistate`.

## Details

If the optional parameter `cell_ids` is used, this method find the most
recent common ancestors of the cells having an identifier among those in
`cell_ids`. If, otherwise, the optional parameter is not used, this
method find the most recent common ancestors of the forest leaves.

## See also

[`PhylogeneticForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest.md)

## Examples

``` r
# use a forest example
forest <- example("PhylogeneticForest")

# get the most recent common ancestor of all the leaves in the forest
forest$get_coalescent_cells()
#>   cell_id ancestor depth mutant epistate sample birth_time
#> 1       1       NA     0      A       E1   <NA>          0
```
