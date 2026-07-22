# Getting forest species

This method builds a data frame containing information about the
simulated species.

## Value

A data frame reporting `mutant` and, if the simulation has epigenetic
states, `epistate` for each registered species.

## See also

[`MutationEngine$get_species_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_species_info.md),
[`PhylogeneticForest$get_species_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/PhylogeneticForest-cash-get_species_info.md),
[`SampleForest`](https://caravagnalab.github.io/ProCESS/1.3/reference/SampleForest.md)

## Examples

``` r
# use a forest example without epistates
forest <- example("SampleForest - no epistates")

# get species information. Since the simulation has no epigenetic
# state, the species correspond to the mutants
forest$get_species_info()
#>   mutant
#> 1      A
#> 2      B

# use a forest example
forest <- example("SampleForest")

# get species information
forest$get_species_info()
#>   mutant epistate
#> 1      A       E1
#> 2      A       E2
#> 3      B       E1
#> 4      B       E2
```
