# Getting the epigenetic states in the simulation

This method returns the epigenetic states in the simulations.

## Value

A data frame having a single column `epistate`. The column contains the
names of the epigenetic states added to the simulation.

## See also

[`TissueSimulation$add_epigenetic_state()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_epigenetic_state.md),
[`TissueSimulation$add_epigenetic_states()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_epigenetic_states.md)
[`TissueSimulation$get_mutants()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_mutants.md)

## Examples

``` r
# set the seed of the random number generator
set.seed(0)

# create a simulation
sim <- TissueSimulation()

# get the data frame of the mutants
sim$get_epigenetic_states()
#> [1] epistate
#> <0 rows> (or 0-length row.names)

# add epigenetic states
sim$add_epigenetic_states(c("E1","E2"))

# get the data frame of the mutants
sim$get_epigenetic_states()
#>   epistate
#> 1       E1
#> 2       E2
```
