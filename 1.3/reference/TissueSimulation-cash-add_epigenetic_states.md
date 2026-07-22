# Adding epigenetic states and their species

This method adds epigenetic states and their species to the simulation.

## Arguments

- epigenetic_states:

  A list of epigenetic state names

## Details

This method introduces novel epigenetic states into the simulation.
Additionally, the method adds to each known mutant as many species as
the number of newly introduced epigenetic states. The default rate of
these new species is set to zero.

## See also

[`TissueSimulation$add_epigenetic_state()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_epigenetic_state.md),
[`TissueSimulation$get_epigenetic_states()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_epigenetic_states.md),
[`TissueSimulation$add_mutant()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_mutant.md),
[`TissueSimulation$add_mutants()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_mutants.md),
[`TissueSimulation$set_rate()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rate.md),
[`TissueSimulation$set_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rates.md),
[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

sim$get_epigenetic_states()
#> [1] epistate
#> <0 rows> (or 0-length row.names)

# add the epigenetic state "E1", "E2", and "E3" to the simulation.
sim$add_epigenetic_states(c("E1", "E2", "E3"))

sim$get_epigenetic_states()
#>   epistate
#> 1       E1
#> 2       E2
#> 3       E3
```
