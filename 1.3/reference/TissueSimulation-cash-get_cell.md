# Getting one of the tissue cells

This method collects some data of the aimed cell without altering the
tissue.

## Arguments

- x:

  The position of the aimed cell on the x axis.

- y:

  The position of the aimed cell on the y axis.

## Value

A data frame reporting `cell_id`, `mutant`, `position_x`, and
`position_y` of the aimed cell. If the simulation has epigenetic states,
the data frame also contains the column `epistate`.

## See also

[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

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


# collect a cell in the tissue
sim$get_cell(501, 502)
#>   cell_id mutant position_x position_y
#> 1    1491      A        501        502
```
