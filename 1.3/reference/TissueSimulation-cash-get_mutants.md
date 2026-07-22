# Getting the names of the simulated mutants

This method returns the names of the simulated mutants.

## Value

A data frame containing the column `mutant`. Each row of the data frame
reports the name of one of the simulated mutants.

## See also

[`TissueSimulation$add_mutant()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_mutant.md),
[`TissueSimulation$add_mutants()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-add_mutants.md),
[`TissueSimulation$get_epigenetic_states()`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation-cash-get_epigenetic_states.md),
[`TissueSimulation`](https://caravagnalab.github.io/ProCESS/1.3/reference/TissueSimulation_class.md)

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
