# Getting the species and their rates

This method returns the species and their rates.

## Value

A data frame reporting `species`, `time`, `SNV_rate`, `indel_rate`, and
`CNA_rate` for each species.

## Details

This method returns the species and their rates during the simulation
are returned in a data frame. The column `species` contains the species
names; the columns `time`, `SNV_rate`, `indel_rate`, and `CNA_rate`
store the time from which rates hold, and the corresponding the SNV,
indel, and CNA rates, respectively.

## See also

[`MutationEngine$get_species_info()`](https://caravagnalab.github.io/ProCESS/1.3/reference/MutationEngine-cash-get_species_info.md)

## Examples

``` r
# use a phylogenetic forest example
forest <- example("PhylogeneticForest")

# show information about samples
forest$get_species_info()
#>   mutant epistate time SNV_rate indel_rate CNA_rate
#> 1      A       E1    0    1e-09      1e-09    0e+00
#> 2      A       E2    0    3e-09      0e+00    1e-11
#> 3      B       E1    0    8e-09      0e+00    0e+00
#> 4      B       E2    0    2e-08      0e+00    0e+00
```
