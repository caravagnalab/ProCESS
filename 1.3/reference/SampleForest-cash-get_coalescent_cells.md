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

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add the mutant "A"
sim$add_mutant("A", c(duplication = 0.2, death = 0.01))

# place a cell in the tissue
sim$place_cell("A", 500, 500)

# run the simulation until "A" has less than 50000 cells
sim$run_up_to_size("A", 50000)
#> 
 [█████████████████████████████████████---] 92% [00m:00s] Cells: 46431                                        

 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# sample the region [450,500]x[475,550]
sim$sample_cells("S1", lower_corner = c(450, 475),
                 upper_corner = c(500, 550))

# build the sample forest
forest <- sim$get_sample_forest()

# get the most recent common ancestor of all the leaves in the forest
forest$get_coalescent_cells()
#>   cell_id ancestor depth mutant sample birth_time
#> 1       1       NA     0      A   <NA>          0
```
