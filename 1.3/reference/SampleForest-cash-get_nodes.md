# Getting forest nodes

This method builds a data frame containing forest nodes.

## Value

A data frame reporting, for each node in the forest, the identified
(column `cell_id`), the ancestor identifier (column `ancestor`), the
node's depth (column `depth`), the name of the sample containing the
node, (column `sample`), and the mutant (column `mutant`), the birth
time (column `birth_time`). Whenever, the simulation has epigenetic
states, the data frame also contains the column `epistate`.

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
 [█████████████████████████████████████---] 91% [00m:00s] Cells: 45601                                        

 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                    


# sample the region [450,500]x[475,550]
sim$sample_cells("S1", lower_corner = c(450, 475),
                 upper_corner = c(500, 550))

# build the sample forest
forest <- sim$get_sample_forest()

# print the first five nodes
forest$get_nodes() %>% head(5)
#>   cell_id ancestor depth mutant sample birth_time
#> 1       1       NA     0      A   <NA>   0.000000
#> 2       2        1     1      A   <NA>   2.870718
#> 3       3        1     1      A   <NA>   2.870718
#> 4       4        3     2      A   <NA>   4.616181
#> 5       5        3     2      A   <NA>   4.616181
```
