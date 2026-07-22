# Adding an epigenetic state and its species

This method adds an epigenetic state and its species to the simulation.

## Arguments

- epigenetic_state:

  The name of the epigenetic state to add.

## Details

This method introduces a new epigenetic state into the simulation.
Additionally, the method adds to each known mutant a new species. The
default rate of the new species is set to zero.

## See also

[`TissueSimulation$add_epigenetic_states()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_epigenetic_states.md),
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

# add the epigenetic state "E1" to the simulation.
sim$add_epigenetic_state("E1")

sim$get_epigenetic_states()
#>   epistate
#> 1       E1
```
