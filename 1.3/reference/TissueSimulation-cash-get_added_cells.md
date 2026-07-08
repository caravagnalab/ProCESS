# Getting the cells manually added to the simulation

This method returns the cells manually added to the simulation.

## Value

A data frame having the columns `mutant`, `position_x`, `position_y`,
`time`, and, when the simulation has epigenetic states, `epistate`. The
data frame contains a row for each cell manually added to the
simulation.

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# add mutant "A" and set its rates
sim$add_mutant("A", c(duplication = 0.2, death = 0.1))

# place an "A" cell in the tissue
sim$place_cell("A", 500, 500)

# run the simulation up to time 70
sim$run_up_to_time(70)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                            


# get the cells
sim$get_added_cells()
#>   mutant position_x position_y time
#> 1      A        500        500    0
```
