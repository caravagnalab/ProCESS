# Adding mutants and their species

This method adds mutants and their species to the simulation.

## Arguments

- mutants:

  A list of mutant names

## Details

This method adds mutants to the simulation. The method also creates the
species associated to the new mutants according to the known epigenetic
states. The default rate of the new species is set to zero.

## See also

[`TissueSimulation$add_mutant()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_mutant.md),
[`TissueSimulation$add_epigenetic_state()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_epigenetic_state.md),
[`TissueSimulation$add_epigenetic_states()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_epigenetic_states.md),
[`TissueSimulation$get_mutants()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_mutants.md),
[`TissueSimulation$set_rate()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rate.md),
[`TissueSimulation$set_rates()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-set_rates.md)

## Examples

``` r
# create a simulation
sim <- TissueSimulation()

sim$get_mutants()
#> [1] mutant
#> <0 rows> (or 0-length row.names)

# add the mutants "A", "B", and "C" to the simulation.
sim$add_mutants(c("A", "B", "C"))

sim$get_mutants()
#>   mutant
#> 1      A
#> 2      B
#> 3      C
```
