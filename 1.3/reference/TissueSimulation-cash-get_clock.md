# Getting the simulated time

This method returns the current simulation time.

## Value

The time simulated by the simulation.

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

# run the simulation up to time 40
sim$run_up_to_time(40)
#> 
 [████████████████████████████████████████] 100% [00m:00s] Saving snapshot                                                                                                       


# get the simulated time
sim$get_clock()
#> [1] 40.06212
```
